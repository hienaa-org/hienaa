package crt

import (
	"slices"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type mulType byte

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

	pool *pool.Pool[*Element]
}

// newTrivialMulOperator creates a new [trivialMulOperator].
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

		pool: pool.NewPool(func() *Element {
			return NewPoly(params.Rank(), len(ambMod), dft.FormCoeff)
		}),
	}
}

func (op *trivialMulOperator) mulTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		var p0Amb, p1Amb *Element
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
				vec.MMulTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(eOut.Coeffs[i:i+1], p0Amb.Coeffs[:op.ambModLen[i]])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.Form = p.Form
	}
}

func (op *trivialMulOperator) mulAddTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if eOut.Form != dft.FormNTT || e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		var p0Amb, p1Amb *Element
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
				vec.MMulAddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0], op.mod[i])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		if p.Form != eOut.Form {
			panic("inconsistent forms")
		}

		for i := range op.mod {
			vec.MulAddTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
	}
}

func (op *trivialMulOperator) mulSubTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if eOut.Form != dft.FormNTT || e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		var p0Amb, p1Amb *Element
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
				vec.MMulSubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0], op.mod[i])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		if p.Form != eOut.Form {
			panic("inconsistent forms")
		}

		for i := range op.mod {
			vec.MulSubTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
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

func (op *trivialMulOperator) append(opIn *trivialMulOperator) *trivialMulOperator {
	return &trivialMulOperator{
		params: opIn.params,
		mod:    slices.Concat(op.mod, opIn.mod),

		expFac:    op.expFac,
		ambModLen: slices.Concat(op.ambModLen, opIn.ambModLen),
		ambMod:    opIn.ambMod,
		ambNTT:    opIn.ambNTT,
		embedder:  slices.Concat(op.embedder, opIn.embedder),

		pool: op.pool,
	}
}

func (op *trivialMulOperator) appendTmpMod(mod *num.Modulus) *trivialMulOperator {
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

	pool *pool.Pool[*Element]
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

		pool: pool.NewPool(func() *Element {
			return NewPoly(ambParams.Rank(), len(ambMod), dft.FormCoeff)
		}),
	}
}

func (op *anyCyclotomicMulOperator) mulTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		var p0Amb, p1Amb *Element
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
				vec.MMulTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0.Coeffs[i])
					clear(p0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1.Coeffs[i])
					clear(p1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOut.Coeffs[i], p0Amb.Coeffs[0], i)
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.Form = p.Form
	}
}

func (op *anyCyclotomicMulOperator) mulAddTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if eOut.Form != dft.FormNTT || e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		var p0Amb, p1Amb *Element
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
				vec.MMulAddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0.Coeffs[i])
					clear(p0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1.Coeffs[i])
					clear(p1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:rank], p0Amb.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0][:rank], op.mod[i])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		if p.Form != eOut.Form {
			panic("inconsistent forms")
		}

		for i := range op.mod {
			vec.MulAddTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
	}
}

func (op *anyCyclotomicMulOperator) mulSubTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if eOut.Form != dft.FormNTT || e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		var p0Amb, p1Amb *Element
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
				vec.MMulSubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0.Coeffs[i])
					clear(p0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1.Coeffs[i])
					clear(p1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:rank], p0Amb.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0][:rank], op.mod[i])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		if p.Form != eOut.Form {
			panic("inconsistent forms")
		}

		for i := range op.mod {
			vec.MulSubTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
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

		reducer: op.reducer.WithModSlice(lo, hi),

		pool: op.pool,
	}
}

func (op *anyCyclotomicMulOperator) append(opIn *anyCyclotomicMulOperator) *anyCyclotomicMulOperator {
	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    slices.Concat(op.mod, opIn.mod),

		ambModLen: slices.Concat(op.ambModLen, opIn.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, opIn.embedder),

		reducer: op.reducer.Append(opIn.reducer),

		pool: op.pool,
	}
}

func (op *anyCyclotomicMulOperator) appendTmpMod(mod *num.Modulus) *anyCyclotomicMulOperator {
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

		reducer: op.reducer.AppendTmpMod(mod),

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

	pool *pool.Pool[*Element]
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

		pool: pool.NewPool(func() *Element {
			return NewPoly(ambParams.Rank(), len(ambMod), dft.FormCoeff)
		}),
	}
}

func (op *reduceMulOperator) mulTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		p0Amb := op.pool.Get()
		p1Amb := op.pool.Get()
		defer op.pool.Put(p0Amb)
		defer op.pool.Put(p1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(p0Amb.Coeffs[0], e0.Coeffs[i])
				clear(p0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])

				copy(p1Amb.Coeffs[0], e1.Coeffs[i])
				clear(p1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p1Amb.Coeffs[0], p1Amb.Coeffs[0])

				vec.MMulTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0], p1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])
				op.reducer.reduceTo(eOut.Coeffs[i], p0Amb.Coeffs[0], i)
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0.Coeffs[i])
					clear(p0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1.Coeffs[i])
					clear(p1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOut.Coeffs[i], p0Amb.Coeffs[0], i)
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.Form = p.Form
	}
}

func (op *reduceMulOperator) mulAddTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if eOut.Form != dft.FormNTT || e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		p0Amb := op.pool.Get()
		p1Amb := op.pool.Get()
		defer op.pool.Put(p0Amb)
		defer op.pool.Put(p1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(p0Amb.Coeffs[0], e0.Coeffs[i])
				clear(p0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])

				copy(p1Amb.Coeffs[0], e1.Coeffs[i])
				clear(p1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p1Amb.Coeffs[0], p1Amb.Coeffs[0])

				vec.MMulTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0], p1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0.Coeffs[i])
					clear(p0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1.Coeffs[i])
					clear(p1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		if p.Form != eOut.Form {
			panic("inconsistent forms")
		}

		for i := range op.mod {
			vec.MulAddTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
	}
}

func (op *reduceMulOperator) mulSubTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		if eOut.Form != dft.FormNTT || e0.Form != dft.FormNTT || e1.Form != dft.FormNTT {
			panic("inputs must be in NTT form")
		}

		p0Amb := op.pool.Get()
		p1Amb := op.pool.Get()
		defer op.pool.Put(p0Amb)
		defer op.pool.Put(p1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(p0Amb.Coeffs[0], e0.Coeffs[i])
				clear(p0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])

				copy(p1Amb.Coeffs[0], e1.Coeffs[i])
				clear(p1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(p1Amb.Coeffs[0], p1Amb.Coeffs[0])

				vec.MMulTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0], p1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(p0Amb.Coeffs[0], p0Amb.Coeffs[0])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(p0Amb.Coeffs[j], e0.Coeffs[i])
					clear(p0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])

					copy(p1Amb.Coeffs[j], e1.Coeffs[i])
					clear(p1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(p1Amb.Coeffs[j], p1Amb.Coeffs[j])

					vec.MMulTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j], p1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(p0Amb.Coeffs[j], p0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedTo(p0Amb.Coeffs[:1], p0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(p0Amb.Coeffs[0][:op.rank], p0Amb.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], p0Amb.Coeffs[0][:op.rank], op.mod[i])
			}
		}

		eOut.Form = dft.FormNTT

	default:
		c, p := orderByType(e0, e1)
		if p.Form != eOut.Form {
			panic("inconsistent forms")
		}

		for i := range op.mod {
			vec.MulSubTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
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

		reducer: op.reducer.WithModSlice(lo, hi),

		pool: op.pool,
	}
}

func (op *reduceMulOperator) append(opIn *reduceMulOperator) *reduceMulOperator {
	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       slices.Concat(op.mod, opIn.mod),

		ntt: slices.Concat(op.ntt, opIn.ntt),

		ambModLen: slices.Concat(op.ambModLen, opIn.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  slices.Concat(op.embedder, opIn.embedder),

		reducer: op.reducer.Append(opIn.reducer),

		pool: op.pool,
	}
}

func (op *reduceMulOperator) appendTmpMod(mod *num.Modulus) *reduceMulOperator {
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

		reducer: op.reducer.AppendTmpMod(mod),

		pool: op.pool,
	}
}
