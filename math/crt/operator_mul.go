package crt

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// mulOperator implements multiplication operations.
type mulOperator interface {
	// Mul returns e0 * e1.
	// When e0, e1 are both polynomials, they must be in NTT form.
	Mul(e0, e1 *Element) *Element
	// MulTo computes eOut = e0 * e1.
	// When e0, e1 are both polynomials, they must be in NTT form.
	MulTo(eOut, e0, e1 *Element)
	// MulAddTo computes eOut += e0 * e1.
	// When e0, e1 are both polynomials, they must be in NTT form.
	MulAddTo(eOut, e0, e1 *Element)
	// MulSubTo computes eOut -= e0 * e1.
	// When e0, e1 are both polynomials, they must be in NTT form.
	MulSubTo(eOut, e0, e1 *Element)

	withModIdx(idx ...int) mulOperator
	append(op0 mulOperator) mulOperator
	appendAuxModulus(mod *num.Modulus) mulOperator
}

// baseMulOperator is a [mulOperator] for rings that do not require reduction after
// multiplication.
// This includes all rings except aribtrary cyclotomic and quotient ring.
type baseMulOperator struct {
	params dft.RingParameters
	mod    []*num.Modulus

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	pool *pool.Pool[*Element]
}

// newBaseMulOperator creates a new [baseMulOperator].
func newBaseMulOperator(params dft.RingParameters, mod []*num.Modulus) *baseMulOperator {
	ambMod := dft.MustFindAmbientPrimes(params, num.Log2(params.ExpandFactor())+2*num.MaxModulusBits)
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(params, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}

		maxBits := num.Log2(params.ExpandFactor()) + 2*num.Log2(mod[i].Value())
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits {
				ambModLen[i] = j + 1
				break
			}
		}
	}

	embedder := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
	}

	return &baseMulOperator{
		params: params,
		mod:    mod,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		pool: pool.NewPool(func() *Element {
			return NewPoly(params.Rank(), len(ambMod))
		}),
	}
}

// Mul returns e0 * e1.
// When e0, e1 are both polynomials, they must be in NTT form.
func (op *baseMulOperator) Mul(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
// When e0, e1 are both polynomials, they must be in NTT form.
func (op *baseMulOperator) MulTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var e0Amb, e1Amb *Element
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				e0Amb = op.pool.Get()
				e1Amb = op.pool.Get()
				defer op.pool.Put(e0Amb)
				defer op.pool.Put(e1Amb)
				break
			}
		}

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(eOut.Coeffs[i:i+1], e0Amb.Coeffs[:op.ambModLen[i]])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

// MulAddTo computes eOut += e0 * e1.
func (op *baseMulOperator) MulAddTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var e0Amb, e1Amb *Element
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				e0Amb = op.pool.Get()
				e1Amb = op.pool.Get()
				defer op.pool.Put(e0Amb)
				defer op.pool.Put(e1Amb)
				break
			}
		}

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulAddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0], op.mod[i])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulAddScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

// MulSubTo computes eOut -= e0 * e1.
func (op *baseMulOperator) MulSubTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var e0Amb, e1Amb *Element
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				e0Amb = op.pool.Get()
				e1Amb = op.pool.Get()
				defer op.pool.Put(e0Amb)
				defer op.pool.Put(e1Amb)
				break
			}
		}

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulSubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0], op.mod[i])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulSubScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

func (op *baseMulOperator) withModIdx(idx ...int) mulOperator {
	return &baseMulOperator{
		params: op.params,
		mod:    vec.Gather(op.mod, idx...),

		ambModLen: vec.Gather(op.ambModLen, idx...),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Gather(op.embedder, idx...),

		pool: op.pool,
	}
}

func (op *baseMulOperator) append(op0 mulOperator) mulOperator {
	opOther := op0.(*baseMulOperator)
	return &baseMulOperator{
		params: op.params,
		mod:    vec.Concat(op.mod, opOther.mod),

		ambModLen: vec.Concat(op.ambModLen, opOther.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Concat(op.embedder, opOther.embedder),

		pool: op.pool,
	}
}

func (op *baseMulOperator) appendAuxModulus(mod *num.Modulus) mulOperator {
	ambModLen := 0
	ambBits := num.Log2(op.params.ExpandFactor()) + 2*num.Log2(mod.Value())
	bits := 0.0
	for j := range op.ambMod {
		bits += num.Log2(op.ambMod[j].Value())
		if bits >= ambBits {
			ambModLen = j + 1
			break
		}
	}

	embedder := NewEmbedder([]*num.Modulus{mod}, op.ambMod[:ambModLen])

	return &baseMulOperator{
		params: op.params,
		mod:    vec.Concat(op.mod, []*num.Modulus{mod}),

		ambModLen: vec.Concat(op.ambModLen, []int{ambModLen}),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Concat(op.embedder, []*Embedder{embedder}),

		pool: op.pool,
	}
}

// anyCyclotomicMulOperator is a [mulOperator] for arbitrary cyclotomic rings.
type anyCyclotomicMulOperator struct {
	params dft.RingParameters
	mod    []*num.Modulus

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	reducer *CyclotomicReducer

	pool *pool.Pool[*Element]
}

// newAnyCyclotomicMulOperator creates a new [anyCyclotomicMulOperator].
func newAnyCyclotomicMulOperator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) *anyCyclotomicMulOperator {
	ambParams := dft.NewCyclicParameters(params.CycloOrder())
	ambMod := dft.MustFindAmbientPrimes(ambParams, num.Log2(params.ExpandFactor())+2*num.MaxModulusBits)
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}

		maxBits := num.Log2(params.ExpandFactor()) + 2*num.Log2(mod[i].Value())
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits {
				ambModLen[i] = j + 1
				break
			}
		}
	}

	embedder := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
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
			return NewPoly(params.CycloOrder(), len(ambMod))
		}),
	}
}

// Mul returns e0 * e1.
// When e0, e1 are both polynomials, they must be in NTT form.
func (op *anyCyclotomicMulOperator) Mul(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
// When e0, e1 are both polynomials, they must be in NTT form.
func (op *anyCyclotomicMulOperator) MulTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var e0Amb, e1Amb *Element
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				e0Amb = op.pool.Get()
				e1Amb = op.pool.Get()
				defer op.pool.Put(e0Amb)
				defer op.pool.Put(e1Amb)
				break
			}
		}

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(e0Amb.Coeffs[j], e0.Coeffs[i])
					clear(e0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])

					copy(e1Amb.Coeffs[j], e1.Coeffs[i])
					clear(e1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1Amb.Coeffs[j])

					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOut.Coeffs[i], e0Amb.Coeffs[0], i)
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

// MulAddTo computes eOut += e0 * e1.
func (op *anyCyclotomicMulOperator) MulAddTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var e0Amb, e1Amb *Element
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				e0Amb = op.pool.Get()
				e1Amb = op.pool.Get()
				defer op.pool.Put(e0Amb)
				defer op.pool.Put(e1Amb)
				break
			}
		}

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulAddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(e0Amb.Coeffs[j], e0.Coeffs[i])
					clear(e0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])

					copy(e1Amb.Coeffs[j], e1.Coeffs[i])
					clear(e1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1Amb.Coeffs[j])

					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(e0Amb.Coeffs[0][:rank], e0Amb.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0][:rank], op.mod[i])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulAddScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

// MulSubTo computes eOut -= e0 * e1.
func (op *anyCyclotomicMulOperator) MulSubTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		var e0Amb, e1Amb *Element
		for i := range op.ambModLen {
			if op.ambModLen[i] > 0 {
				e0Amb = op.pool.Get()
				e1Amb = op.pool.Get()
				defer op.pool.Put(e0Amb)
				defer op.pool.Put(e1Amb)
				break
			}
		}

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulSubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(e0Amb.Coeffs[j], e0.Coeffs[i])
					clear(e0Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])

					copy(e1Amb.Coeffs[j], e1.Coeffs[i])
					clear(e1Amb.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1Amb.Coeffs[j])

					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(e0Amb.Coeffs[0][:rank], e0Amb.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0][:rank], op.mod[i])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulSubScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

func (op *anyCyclotomicMulOperator) withModIdx(idx ...int) mulOperator {
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

func (op *anyCyclotomicMulOperator) append(op0 mulOperator) mulOperator {
	opOther := op0.(*anyCyclotomicMulOperator)
	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    vec.Concat(op.mod, opOther.mod),

		ambModLen: vec.Concat(op.ambModLen, opOther.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Concat(op.embedder, opOther.embedder),

		reducer: op.reducer.Append(opOther.reducer),

		pool: op.pool,
	}
}

func (op *anyCyclotomicMulOperator) appendAuxModulus(mod *num.Modulus) mulOperator {
	ambModLen := 0
	ambBits := num.Log2(op.params.ExpandFactor()) + 2*num.Log2(mod.Value())
	bits := 0.0
	for j := range op.ambMod {
		bits += num.Log2(op.ambMod[j].Value())
		if bits >= ambBits {
			ambModLen = j + 1
			break
		}
	}

	embedder := NewEmbedder([]*num.Modulus{mod}, op.ambMod[:ambModLen])

	return &anyCyclotomicMulOperator{
		params: op.params,
		mod:    vec.Concat(op.mod, []*num.Modulus{mod}),

		ambModLen: vec.Concat(op.ambModLen, []int{ambModLen}),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Concat(op.embedder, []*Embedder{embedder}),

		reducer: op.reducer.AppendAuxModulus(mod),

		pool: op.pool,
	}
}

// reduceMulOperator is a [mulOperator] for arbitrary modulo rings.
type reduceMulOperator struct {
	rank      int
	ambParams dft.RingParameters
	mod       []*num.Modulus

	ntt []dft.Transformer

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	reducer *Reducer

	pool *pool.Pool[*Element]
}

// newReduceMulOperator creates a new [reduceMulOperator].
func newReduceMulOperator(mod []*num.Modulus, modPoly []int64, reducer *Reducer) *reduceMulOperator {
	ambParams := dft.NewCyclicParameters(num.NextProdPower(2*len(modPoly)-1, []int{2}))
	ambMod := dft.MustFindAmbientPrimes(ambParams, num.Log2(ambParams.ExpandFactor())+2*num.MaxModulusBits)
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	ntt := make([]dft.Transformer, len(mod))
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(ambParams, mod[i]) {
			ntt[i] = dft.NewTransformer(ambParams, mod[i])
			continue
		}

		maxBits := num.Log2(ambParams.ExpandFactor()) + 2*num.Log2(mod[i].Value())
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits {
				ambModLen[i] = j + 1
				break
			}
		}
	}

	embedder := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
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
			return NewPoly(ambParams.Rank(), len(ambMod))
		}),
	}
}

// Mul returns e0 * e1.
// When e0, e1 are both polynomials, they must be in NTT form.
func (op *reduceMulOperator) Mul(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
// When e0, e1 are both polynomials, they must be in NTT form.
func (op *reduceMulOperator) MulTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		e0Amb := op.pool.Get()
		e1Amb := op.pool.Get()
		defer op.pool.Put(e0Amb)
		defer op.pool.Put(e1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(e0Amb.Coeffs[0], e0.Coeffs[i])
				clear(e0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0])

				copy(e1Amb.Coeffs[0], e1.Coeffs[i])
				clear(e1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(e1Amb.Coeffs[0], e1Amb.Coeffs[0])

				vec.MMulTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0], e1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0])
				op.reducer.reduceTo(eOut.Coeffs[i], e0Amb.Coeffs[0], i)
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(e0Amb.Coeffs[j], e0.Coeffs[i])
					clear(e0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])

					copy(e1Amb.Coeffs[j], e1.Coeffs[i])
					clear(e1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1Amb.Coeffs[j])

					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOut.Coeffs[i], e0Amb.Coeffs[0], i)
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

// MulAddTo computes eOut += e0 * e1.
func (op *reduceMulOperator) MulAddTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		e0Amb := op.pool.Get()
		e1Amb := op.pool.Get()
		defer op.pool.Put(e0Amb)
		defer op.pool.Put(e1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(e0Amb.Coeffs[0], e0.Coeffs[i])
				clear(e0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0])

				copy(e1Amb.Coeffs[0], e1.Coeffs[i])
				clear(e1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(e1Amb.Coeffs[0], e1Amb.Coeffs[0])

				vec.MMulTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0], e1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0])
				op.reducer.reduceTo(e0Amb.Coeffs[0][:op.rank], e0Amb.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(e0Amb.Coeffs[j], e0.Coeffs[i])
					clear(e0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])

					copy(e1Amb.Coeffs[j], e1.Coeffs[i])
					clear(e1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1Amb.Coeffs[j])

					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(e0Amb.Coeffs[0][:op.rank], e0Amb.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0][:op.rank], op.mod[i])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulAddScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

// MulSubTo computes eOut -= e0 * e1.
func (op *reduceMulOperator) MulSubTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		e0Amb := op.pool.Get()
		e1Amb := op.pool.Get()
		defer op.pool.Put(e0Amb)
		defer op.pool.Put(e1Amb)

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(e0Amb.Coeffs[0], e0.Coeffs[i])
				clear(e0Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0])

				copy(e1Amb.Coeffs[0], e1.Coeffs[i])
				clear(e1Amb.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(e1Amb.Coeffs[0], e1Amb.Coeffs[0])

				vec.MMulTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0], e1Amb.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(e0Amb.Coeffs[0], e0Amb.Coeffs[0])
				op.reducer.reduceTo(e0Amb.Coeffs[0][:op.rank], e0Amb.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(e0Amb.Coeffs[j], e0.Coeffs[i])
					clear(e0Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])

					copy(e1Amb.Coeffs[j], e1.Coeffs[i])
					clear(e1Amb.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(e1Amb.Coeffs[j], e1Amb.Coeffs[j])

					vec.MMulTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j], e1Amb.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(e0Amb.Coeffs[j], e0Amb.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(e0Amb.Coeffs[:1], e0Amb.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(e0Amb.Coeffs[0][:op.rank], e0Amb.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], e0Amb.Coeffs[0][:op.rank], op.mod[i])
			}
		}

		eOut.IsNTT = true
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			vec.MulSubScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
		}
		eOut.IsNTT = p.IsNTT
	}
}

func (op *reduceMulOperator) withModIdx(idx ...int) mulOperator {
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

func (op *reduceMulOperator) append(op0 mulOperator) mulOperator {
	opOther := op0.(*reduceMulOperator)
	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       vec.Concat(op.mod, opOther.mod),

		ntt: vec.Concat(op.ntt, opOther.ntt),

		ambModLen: vec.Concat(op.ambModLen, opOther.ambModLen),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Concat(op.embedder, opOther.embedder),

		reducer: op.reducer.Append(opOther.reducer),

		pool: op.pool,
	}
}

func (op *reduceMulOperator) appendAuxModulus(mod *num.Modulus) mulOperator {
	ambModLen := 0
	ambBits := num.Log2(op.ambParams.ExpandFactor()) + 2*num.Log2(mod.Value())
	bits := 0.0
	for j := range op.ambMod {
		bits += num.Log2(op.ambMod[j].Value())
		if bits >= ambBits {
			ambModLen = j + 1
			break
		}
	}

	embedder := NewEmbedder([]*num.Modulus{mod}, op.ambMod[:ambModLen])

	return &reduceMulOperator{
		rank:      op.rank,
		ambParams: op.ambParams,
		mod:       vec.Concat(op.mod, []*num.Modulus{mod}),

		ntt: vec.Concat(op.ntt, []dft.Transformer{nil}),

		ambModLen: vec.Concat(op.ambModLen, []int{ambModLen}),
		ambMod:    op.ambMod,
		ambNTT:    op.ambNTT,
		embedder:  vec.Concat(op.embedder, []*Embedder{embedder}),

		reducer: op.reducer.AppendAuxModulus(mod),

		pool: op.pool,
	}
}
