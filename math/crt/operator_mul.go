package crt

import (
	"sync"

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
}

// baseMulOperator is a [mulOperator] for rings that do not require reduction after
// multiplication.
// This includes all rings except aribtrary cyclotomic and quotient ring.
type baseMulOperator struct {
	rank int
	mod  []*num.Modulus

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	pool *sync.Pool
}

// newBaseMulOperator creates a new [baseMulOperator].
func newBaseMulOperator(params dft.RingParameters, mod []*num.Modulus) baseMulOperator {
	maxBits := make([]float64, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		maxBits[i] = 2 * num.Log2(mod[i].Value())
		switch params.RingType() {
		case dft.TypeCyclic, dft.TypeCyclotomic:
			maxBits[i] += num.Log2(params.Rank())
		case dft.TypeAutFixed:
			maxBits[i] += num.Log2(params.CycloOrder())
		}
	}

	ambMod := dft.MustFindAmbientPrimes(params, vec.Max(maxBits))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(params, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if maxBits[i] == 0 {
			continue
		}
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits[i] {
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

	return baseMulOperator{
		rank: params.Rank(),
		mod:  mod,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		pool: &sync.Pool{
			New: func() any {
				return NewPoly(params.Rank(), max(1, vec.Max(ambModLen)))
			},
		},
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
	checkTernaryOperable(op.rank, len(op.mod), eOut, e0, e1)

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
				e0Amb = op.pool.Get().(*Element)
				e1Amb = op.pool.Get().(*Element)
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
	checkTernaryOperable(op.rank, len(op.mod), eOut, e0, e1)

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
				e0Amb = op.pool.Get().(*Element)
				e1Amb = op.pool.Get().(*Element)
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
	checkTernaryOperable(op.rank, len(op.mod), eOut, e0, e1)

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
				e0Amb = op.pool.Get().(*Element)
				e1Amb = op.pool.Get().(*Element)
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

func (op *baseMulOperator) subOperator(idx ...int) baseMulOperator {
	modCopy := make([]*num.Modulus, len(idx))
	ambModLenCopy := make([]int, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		ambModLenCopy[i] = op.ambModLen[idx[i]]
	}

	embedderCopy := make([]*Embedder, len(idx))
	for i := range idx {
		embedderCopy[i] = op.embedder[idx[i]]
	}

	maxAmbModLen := vec.Max(ambModLenCopy)

	return baseMulOperator{
		rank: op.rank,
		mod:  modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    op.ambMod[:maxAmbModLen],
		ambNTT:    op.ambNTT[:maxAmbModLen],
		embedder:  embedderCopy,

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

	pool *sync.Pool
}

// newAnyCyclotomicMulOperator creates a new [anyCyclotomicMulOperator].
func newAnyCyclotomicMulOperator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) anyCyclotomicMulOperator {
	maxBits := make([]float64, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		maxBits[i] = num.Log2(params.CycloOrder()) + 2*num.Log2(mod[i].Value())
	}

	ambParams := dft.NewCyclicParameters(params.CycloOrder())
	ambMod := dft.MustFindAmbientPrimes(ambParams, vec.Max(maxBits))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if maxBits[i] == 0 {
			continue
		}
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits[i] {
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

	return anyCyclotomicMulOperator{
		params: params,
		mod:    mod,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		reducer: reducer,

		pool: &sync.Pool{
			New: func() any {
				return NewPoly(params.CycloOrder(), max(1, vec.Max(ambModLen)))
			},
		},
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
	checkTernaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

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
				e0Amb = op.pool.Get().(*Element)
				e1Amb = op.pool.Get().(*Element)
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
	checkTernaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

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
				e0Amb = op.pool.Get().(*Element)
				e1Amb = op.pool.Get().(*Element)
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
	checkTernaryOperable(op.params.Rank(), len(op.mod), eOut, e0, e1)

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
				e0Amb = op.pool.Get().(*Element)
				e1Amb = op.pool.Get().(*Element)
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

func (op *anyCyclotomicMulOperator) subOperator(idx ...int) anyCyclotomicMulOperator {
	modCopy := make([]*num.Modulus, len(idx))
	ambModLenCopy := make([]int, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		ambModLenCopy[i] = op.ambModLen[idx[i]]
	}

	maxAmbModLen := vec.Max(ambModLenCopy)

	embedderCopy := make([]*Embedder, len(idx))
	for i := range idx {
		embedderCopy[i] = op.embedder[idx[i]]
	}

	return anyCyclotomicMulOperator{
		params: op.params,
		mod:    modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    op.ambMod[:maxAmbModLen],
		ambNTT:    op.ambNTT[:maxAmbModLen],
		embedder:  embedderCopy,

		reducer: op.reducer.SubReducer(idx...),

		pool: op.pool,
	}
}

// reduceMulOperator is a [mulOperator] for arbitrary modulo rings.
type reduceMulOperator struct {
	rank    int
	ambRank int
	mod     []*num.Modulus

	ntt []dft.Transformer

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	reducer *Reducer

	pool *sync.Pool
}

// newReduceMulOperator creates a new [reduceMulOperator].
func newReduceMulOperator(mod []*num.Modulus, modPoly []int64, reducer *Reducer) reduceMulOperator {
	ambParams := dft.NewCyclicParameters(num.NextProdPower(2*len(modPoly)-1, []int{2}))

	maxBits := make([]float64, len(mod))
	ntt := make([]dft.Transformer, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(ambParams, mod[i]) {
			ntt[i] = dft.NewTransformer(ambParams, mod[i])
			continue
		}
		maxBits[i] = num.Log2(ambParams.Rank()) + 2*num.Log2(mod[i].Value())
	}

	ambMod := dft.MustFindAmbientPrimes(ambParams, vec.Max(maxBits))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if maxBits[i] == 0 {
			continue
		}
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits[i] {
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

	return reduceMulOperator{
		rank:    len(modPoly) - 1,
		ambRank: ambParams.Rank(),
		mod:     mod,

		ntt: ntt,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		reducer: reducer,

		pool: &sync.Pool{
			New: func() any {
				return NewPoly(ambParams.Rank(), max(1, vec.Max(ambModLen)))
			},
		},
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
	checkTernaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		e0Amb := op.pool.Get().(*Element)
		e1Amb := op.pool.Get().(*Element)
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
	checkTernaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		e0Amb := op.pool.Get().(*Element)
		e1Amb := op.pool.Get().(*Element)
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
	checkTernaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], num.Mul(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i]), op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		if !e0.IsNTT || !e1.IsNTT {
			panic("input(s) must be in NTT form")
		}

		e0Amb := op.pool.Get().(*Element)
		e1Amb := op.pool.Get().(*Element)
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

func (op *reduceMulOperator) subOperator(idx ...int) reduceMulOperator {
	modCopy := make([]*num.Modulus, len(idx))
	nttCopy := make([]dft.Transformer, len(idx))
	ambModLenCopy := make([]int, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		nttCopy[i] = op.ntt[idx[i]]
		ambModLenCopy[i] = op.ambModLen[idx[i]]
	}

	maxAmbModLen := vec.Max(ambModLenCopy)

	embedderCopy := make([]*Embedder, len(idx))
	for i := range idx {
		embedderCopy[i] = op.embedder[idx[i]]
	}

	return reduceMulOperator{
		rank:    op.rank,
		ambRank: op.ambRank,
		mod:     modCopy,

		ntt: nttCopy,

		ambModLen: ambModLenCopy,
		ambMod:    op.ambMod[:maxAmbModLen],
		ambNTT:    op.ambNTT[:maxAmbModLen],
		embedder:  embedderCopy,

		reducer: op.reducer.SubReducer(idx...),

		pool: op.pool,
	}
}
