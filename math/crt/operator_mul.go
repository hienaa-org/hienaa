package crt

import (
	"slices"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type mulType uint64

const (
	typeTrivialMul mulType = iota
	typeAnyCyclotomicMul
	typeReduceMul
)

func mulTypeOf(params dft.RingParameters) mulType {
	if params.RingType() == dft.TypeCyclotomic && !num.IsPowerOfTwo(params.CycloIndex()) {
		return typeAnyCyclotomicMul
	}

	if params.RingType() == dft.TypeOther {
		return typeReduceMul
	}

	return typeTrivialMul
}

// trivialMulOperator is a MulOperator for [typeTrivialMul].
type trivialMulOperator struct {
	params dft.RingParameters
	mod    []*num.Modulus

	expFac    uint64
	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []*dft.Transformer
	embedder  []*VecEmbedder

	pool *pool.Pool[*Poly]
}

// newTrivialMulOperator creates a new [baseMulOperator].
func newTrivialMulOperator(params dft.RingParameters, mod []*num.Modulus) *trivialMulOperator {
	var expFac uint64
	switch params.RingType() {
	case dft.TypeCyclotomic, dft.TypeCyclic:
		expFac = uint64(params.Rank())
	case dft.TypeAutFixed:
		if num.IsPowerOfTwo(params.Rank()) {
			expFac = uint64(params.CycloIndex()) >> 1
		} else {
			expFac = 2*uint64(params.CycloIndex()) - 2
		}
	}

	ambMod := ambientNTTPrimes(params, num.Log2(expFac)+2*num.MaxModulusBits)
	ambNTT := make([]*dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(params, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		ambModLen[i] = ambientModLen(ambMod, num.Log2(expFac)+2*num.Log2(mod[i].Value()))
	}

	embedder := make([]*VecEmbedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewVecEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
	}

	return &trivialMulOperator{
		params: params,
		mod:    mod,

		expFac:    expFac,
		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		pool: pool.NewPool(func() *Poly {
			return NewPoly(params.Rank(), len(ambMod))
		}),
	}
}

func (op *trivialMulOperator) mulTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var p0Amb, p1Amb *Poly
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				p0Amb = op.pool.Get()
				p1Amb = op.pool.Get()
				break
			}
		}
		defer func() {
			if p0Amb != nil && p1Amb != nil {
				op.pool.Put(p0Amb)
				op.pool.Put(p1Amb)
			}
		}()

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(eOutPoly.Coeffs[i:i+1], p0Amb.Coeffs[:op.ambModLen[i]])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		for i := range op.mod {
			vec.MulTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
		eOutPoly.IsNTT = p.IsNTT
	}
}

func (op *trivialMulOperator) mulAddTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Add(eOutScalar.Value[i], num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i]), op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !eOutPoly.IsNTT || !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var p0Amb, p1Amb *Poly
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				p0Amb = op.pool.Get()
				p1Amb = op.pool.Get()
				break
			}
		}
		defer func() {
			if p0Amb != nil && p1Amb != nil {
				op.pool.Put(p0Amb)
				op.pool.Put(p1Amb)
			}
		}()

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulAddTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				vec.AddTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0], op.mod[i])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		if p.IsNTT != eOutPoly.IsNTT {
			panic("inconsistent NTT flags")
		}

		for i := range op.mod {
			vec.MulAddTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
	}
}

func (op *trivialMulOperator) mulSubTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Sub(eOutScalar.Value[i], num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i]), op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !eOutPoly.IsNTT || !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var p0Amb, p1Amb *Poly
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				p0Amb = op.pool.Get()
				p1Amb = op.pool.Get()
				break
			}
		}
		defer func() {
			if p0Amb != nil && p1Amb != nil {
				op.pool.Put(p0Amb)
				op.pool.Put(p1Amb)
			}
		}()

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulSubTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				vec.SubTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0], op.mod[i])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		if p.IsNTT != eOutPoly.IsNTT {
			panic("inconsistent NTT flags")
		}

		for i := range op.mod {
			vec.MulSubTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
	}
}

func (op *trivialMulOperator) withModIdx(idx ...int) *trivialMulOperator {
	return &trivialMulOperator{
		params: op.params,
		mod:    vec.Gather(op.mod, idx...),

		expFac:    op.expFac,
		ambModLen: vec.Gather(op.ambModLen, idx...),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Gather(op.embedder, idx...),

		pool: op.pool,
	}
}

func (op *trivialMulOperator) slice(lo, hi int) *trivialMulOperator {
	return &trivialMulOperator{
		params: op.params,
		mod:    op.mod[lo:hi],

		expFac:    op.expFac,
		ambModLen: op.ambModLen[lo:hi],
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  op.embedder[lo:hi],

		pool: op.pool,
	}
}

func (op *trivialMulOperator) append(op0 *trivialMulOperator) *trivialMulOperator {
	return &trivialMulOperator{
		params: op0.params,
		mod:    slices.Concat(op.mod, op0.mod),

		expFac:    op.expFac,
		ambModLen: slices.Concat(op.ambModLen, op0.ambModLen),
		ambMod:    op0.ambMod,
		ambNTT:    op0.ambNTT,
		embedder:  slices.Concat(op.embedder, op0.embedder),

		pool: op.pool,
	}
}

func (op *trivialMulOperator) appendTmpModulus(mod *num.Modulus) *trivialMulOperator {
	ambModLen := ambientModLen(op.ambMod, num.Log2(op.expFac)+2*num.Log2(mod.Value()))
	embedder := NewVecEmbedder([]*num.Modulus{mod}, op.ambMod[:ambModLen])

	return &trivialMulOperator{
		params: op.params,
		mod:    slices.Concat(op.mod, []*num.Modulus{mod}),

		expFac:    op.expFac,
		ambModLen: slices.Concat(op.ambModLen, []int{ambModLen}),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, []*VecEmbedder{embedder}),

		pool: op.pool,
	}
}

// anyCyclotomicMulOperator is a MulOperator for [typeAnyCyclotomicMul].
type anyCyclotomicMulOperator struct {
	params dft.RingParameters
	mod    []*num.Modulus

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []*dft.Transformer
	embedder  []*VecEmbedder

	reducer *CyclotomicReducer

	pool *pool.Pool[*Poly]
}

// newAnyCyclotomicMulOperator creates a new [anyCyclotomicMulOperator].
func newAnyCyclotomicMulOperator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) *anyCyclotomicMulOperator {
	ambParams := dft.NewCyclicParameters(params.CycloIndex())
	ambExpFac := uint64(params.CycloIndex())
	ambMod := ambientNTTPrimes(ambParams, num.Log2(ambExpFac)+2*num.MaxModulusBits)
	ambNTT := make([]*dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		ambModLen[i] = ambientModLen(ambMod, num.Log2(ambExpFac)+2*num.Log2(mod[i].Value()))
	}

	embedder := make([]*VecEmbedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewVecEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
	}

	return &anyCyclotomicMulOperator{
		params: params,
		mod:    mod,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		reducer: reducer,

		pool: pool.NewPool(func() *Poly {
			return NewPoly(ambParams.Rank(), len(ambMod))
		}),
	}
}

func (op *anyCyclotomicMulOperator) mulTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var p0Amb, p1Amb *Poly
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				p0Amb = op.pool.Get()
				p1Amb = op.pool.Get()
				break
			}
		}
		defer func() {
			if p0Amb != nil && p1Amb != nil {
				op.pool.Put(p0Amb)
				op.pool.Put(p1Amb)
			}
		}()

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					clear(p0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					clear(p1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOutPoly.Coeffs[i], p0Amb.Coeffs[0], i)
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		for i := range op.mod {
			vec.MulTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
		eOutPoly.IsNTT = p.IsNTT
	}
}

func (op *anyCyclotomicMulOperator) mulAddTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Add(eOutScalar.Value[i], num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i]), op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !eOutPoly.IsNTT || !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var p0Amb, p1Amb *Poly
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				p0Amb = op.pool.Get()
				p1Amb = op.pool.Get()
				break
			}
		}
		defer func() {
			if p0Amb != nil && p1Amb != nil {
				op.pool.Put(p0Amb)
				op.pool.Put(p1Amb)
			}
		}()

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulAddTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					clear(p0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					clear(p1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:rank], p0Amb.Coeffs[0], i)
				vec.AddTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0][:rank], op.mod[i])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		if p.IsNTT != eOutPoly.IsNTT {
			panic("inconsistent NTT flags")
		}

		for i := range op.mod {
			vec.MulAddTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
	}
}

func (op *anyCyclotomicMulOperator) mulSubTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Sub(eOutScalar.Value[i], num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i]), op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !eOutPoly.IsNTT || !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var p0Amb, p1Amb *Poly
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				p0Amb = op.pool.Get()
				p1Amb = op.pool.Get()
				break
			}
		}
		defer func() {
			if p0Amb != nil && p1Amb != nil {
				op.pool.Put(p0Amb)
				op.pool.Put(p1Amb)
			}
		}()

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulSubTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					clear(p0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					clear(p1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:rank], p0Amb.Coeffs[0], i)
				vec.SubTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0][:rank], op.mod[i])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		if p.IsNTT != eOutPoly.IsNTT {
			panic("inconsistent NTT flags")
		}

		for i := range op.mod {
			vec.MulSubTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
	}
}

func (op *anyCyclotomicMulOperator) withModIdx(idx ...int) *anyCyclotomicMulOperator {
	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    vec.Gather(op.mod, idx...),

		ambModLen: vec.Gather(op.ambModLen, idx...),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Gather(op.embedder, idx...),

		reducer: op.reducer.WithModIdx(idx...),

		pool: op.pool,
	}
}

func (op *anyCyclotomicMulOperator) slice(lo, hi int) *anyCyclotomicMulOperator {
	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    op.mod[lo:hi:hi],

		ambModLen: op.ambModLen[lo:hi],
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  op.embedder[lo:hi],

		reducer: op.reducer.Slice(lo, hi),

		pool: op.pool,
	}
}

func (op *anyCyclotomicMulOperator) append(op0 *anyCyclotomicMulOperator) *anyCyclotomicMulOperator {
	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    slices.Concat(op.mod, op0.mod),

		ambModLen: slices.Concat(op.ambModLen, op0.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, op0.embedder),

		reducer: op.reducer.Append(op0.reducer),

		pool: op.pool,
	}
}

func (op *anyCyclotomicMulOperator) appendTmpModulus(mod *num.Modulus) *anyCyclotomicMulOperator {
	ambExpFac := uint64(op.params.CycloIndex())
	ambModLen := ambientModLen(op.ambMod, num.Log2(ambExpFac)+2*num.Log2(mod.Value()))

	embedder := NewVecEmbedder([]*num.Modulus{mod}, op.ambMod[:ambModLen])

	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    slices.Concat(op.mod, []*num.Modulus{mod}),

		ambModLen: slices.Concat(op.ambModLen, []int{ambModLen}),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, []*VecEmbedder{embedder}),

		reducer: op.reducer.AppendTmpModulus(mod),

		pool: op.pool,
	}
}

// reduceMulOperator is a MulOperator for [typeReduceMul].
type reduceMulOperator struct {
	rank      int
	ambParams dft.RingParameters
	mod       []*num.Modulus

	ntt []*dft.Transformer

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []*dft.Transformer
	embedder  []*VecEmbedder

	reducer *Reducer

	pool *pool.Pool[*Poly]
}

// newReduceMulOperator creates a new [reduceMulOperator].
func newReduceMulOperator(mod []*num.Modulus, modPoly []int64, reducer *Reducer) *reduceMulOperator {
	ambParams := dft.NewCyclicParameters(num.NextProdPower(2*len(modPoly)-1, []int{2}))
	ambMod := ambientNTTPrimes(ambParams, num.Log2(ambParams.Rank())+2*num.MaxModulusBits)
	ambNTT := make([]*dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	ntt := make([]*dft.Transformer, len(mod))
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(ambParams, mod[i]) {
			ntt[i] = dft.NewTransformer(ambParams, mod[i])
			continue
		}
		ambModLen[i] = ambientModLen(ambMod, num.Log2(ambParams.Rank())+2*num.Log2(mod[i].Value()))
	}

	embedder := make([]*VecEmbedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewVecEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
	}

	return &reduceMulOperator{
		rank:      len(modPoly) - 1,
		ambParams: ambParams,
		mod:       mod,

		ntt: ntt,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		reducer: reducer,

		pool: pool.NewPool(func() *Poly {
			return NewPoly(ambParams.Rank(), len(ambMod))
		}),
	}
}

func (op *reduceMulOperator) mulTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		p0Amb := op.pool.Get()
		p1Amb := op.pool.Get()
		defer op.pool.Put(p0Amb)
		defer op.pool.Put(p1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(p0Amb.Coeffs[0], e0Poly.Coeffs[i])
				clear(p0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])

				copy(p1Amb.Coeffs[0], e1Poly.Coeffs[i])
				clear(p1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p1Amb.Coeffs[0], p1Amb.Coeffs[0])

				vec.MMulTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0], p1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])
				op.reducer.reduceTo(eOutPoly.Coeffs[i], p0Amb.Coeffs[0], i)
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					clear(p0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					clear(p1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOutPoly.Coeffs[i], p0Amb.Coeffs[0], i)
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		for i := range op.mod {
			vec.MulTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
		eOutPoly.IsNTT = p.IsNTT
	}
}

func (op *reduceMulOperator) mulAddTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Add(eOutScalar.Value[i], num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i]), op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !eOutPoly.IsNTT || !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		p0Amb := op.pool.Get()
		p1Amb := op.pool.Get()
		defer op.pool.Put(p0Amb)
		defer op.pool.Put(p1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(p0Amb.Coeffs[0], e0Poly.Coeffs[i])
				clear(p0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])

				copy(p1Amb.Coeffs[0], e1Poly.Coeffs[i])
				clear(p1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p1Amb.Coeffs[0], p1Amb.Coeffs[0])

				vec.MMulTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0], p1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.AddTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					clear(p0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					clear(p1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.AddTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		if p.IsNTT != eOutPoly.IsNTT {
			panic("inconsistent NTT flags")
		}

		for i := range op.mod {
			vec.MulAddTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
	}
}

func (op *reduceMulOperator) mulSubTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Sub(eOutScalar.Value[i], num.Mul(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i]), op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		if !eOutPoly.IsNTT || !e0Poly.IsNTT || !e1Poly.IsNTT {
			panic("input(s) must be in NTT form")
		}

		p0Amb := op.pool.Get()
		p1Amb := op.pool.Get()
		defer op.pool.Put(p0Amb)
		defer op.pool.Put(p1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(p0Amb.Coeffs[0], e0Poly.Coeffs[i])
				clear(p0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])

				copy(p1Amb.Coeffs[0], e1Poly.Coeffs[i])
				clear(p1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p1Amb.Coeffs[0], p1Amb.Coeffs[0])

				vec.MMulTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0], p1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.SubTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0Poly.Coeffs[i])
					clear(p0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1Poly.Coeffs[i])
					clear(p1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.SubTo(eOutPoly.Coeffs[i], eOutPoly.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			}
		}

		eOutPoly.IsNTT = true

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		if p.IsNTT != eOutPoly.IsNTT {
			panic("inconsistent NTT flags")
		}

		for i := range op.mod {
			vec.MulSubTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
		}
	}
}

func (op *reduceMulOperator) withModIdx(idx ...int) *reduceMulOperator {
	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       vec.Gather(op.mod, idx...),

		ntt: vec.Gather(op.ntt, idx...),

		ambModLen: vec.Gather(op.ambModLen, idx...),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Gather(op.embedder, idx...),

		reducer: op.reducer.WithModIdx(idx...),

		pool: op.pool,
	}
}

func (op *reduceMulOperator) slice(lo, hi int) *reduceMulOperator {
	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       op.mod[lo:hi],

		ntt: op.ntt[lo:hi],

		ambModLen: op.ambModLen[lo:hi],
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  op.embedder[lo:hi],

		reducer: op.reducer.Slice(lo, hi),

		pool: op.pool,
	}
}

func (op *reduceMulOperator) append(op0 *reduceMulOperator) *reduceMulOperator {
	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       slices.Concat(op.mod, op0.mod),

		ntt: slices.Concat(op.ntt, op0.ntt),

		ambModLen: slices.Concat(op.ambModLen, op0.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, op0.embedder),

		reducer: op.reducer.Append(op0.reducer),

		pool: op.pool,
	}
}

func (op *reduceMulOperator) appendTmpModulus(mod *num.Modulus) *reduceMulOperator {
	ambModLen := ambientModLen(op.ambMod, num.Log2(op.ambParams.Rank())+2*num.Log2(mod.Value()))

	embedder := NewVecEmbedder([]*num.Modulus{mod}, op.ambMod[:ambModLen])

	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       slices.Concat(op.mod, []*num.Modulus{mod}),

		ntt: slices.Concat(op.ntt, []*dft.Transformer{nil}),

		ambModLen: slices.Concat(op.ambModLen, []int{ambModLen}),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, []*VecEmbedder{embedder}),

		reducer: op.reducer.AppendTmpModulus(mod),

		pool: op.pool,
	}
}
