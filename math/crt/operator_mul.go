package crt

import (
	"math"

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

// mulOperatorBuffer is a buffer for [mulOperator].
type mulOperatorBuffer struct {
	p0 *Element
	p1 *Element
}

// newMulOperatorBuffer creates a new [mulOperatorBuffer].
func newMulOperatorBuffer(rank, modLen int) mulOperatorBuffer {
	return mulOperatorBuffer{
		p0: NewPoly(rank, modLen),
		p1: NewPoly(rank, modLen),
	}
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

	buf mulOperatorBuffer
}

// newBaseMulOperator creates a new [baseMulOperator].
func newBaseMulOperator(params dft.RingParameters, mod []*num.Modulus) baseMulOperator {
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		maxBits := 2 * num.Log2(mod[i].Value())
		switch params.RingType() {
		case dft.TypeCyclic, dft.TypeCyclotomic:
			maxBits += num.Log2(params.Rank())
		case dft.TypeAutFixed:
			maxBits += num.Log2(params.CycloOrder())
		}
		ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
	}

	ambMod := dft.MustFindPrevNTTPrimes(params, num.MaxModulusBits, vec.Max(ambModLen))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(params, ambMod[i])
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

		buf: newMulOperatorBuffer(params.Rank(), max(1, vec.Max(ambModLen))),
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(eOut.Coeffs[i:i+1], op.buf.p0.Coeffs[:op.ambModLen[i]])
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulAddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0], op.mod[i])
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulSubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0], op.mod[i])
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

func (e *baseMulOperator) subOperator(idx ...int) baseMulOperator {
	modCopy := make([]*num.Modulus, len(idx))
	ambModLenCopy := make([]int, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		ambModLenCopy[i] = e.ambModLen[idx[i]]
	}

	maxAmbModLen := vec.Max(ambModLenCopy)
	ambNTTCopy := make([]dft.Transformer, maxAmbModLen)
	for i := range ambNTTCopy {
		ambNTTCopy[i] = e.ambNTT[i].SafeCopy()
	}

	embedderCopy := make([]*Embedder, len(idx))
	for i := range idx {
		if e.embedder[idx[i]] != nil {
			embedderCopy[i] = e.embedder[idx[i]].SafeCopy()
		}
	}

	return baseMulOperator{
		rank: e.rank,
		mod:  modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    e.ambMod[:maxAmbModLen],
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		buf: newMulOperatorBuffer(e.rank, max(1, maxAmbModLen)),
	}
}

func (e *baseMulOperator) safeCopy() baseMulOperator {
	ambNTTCopy := make([]dft.Transformer, len(e.ambNTT))
	for i := range e.ambNTT {
		ambNTTCopy[i] = e.ambNTT[i].SafeCopy()
	}

	embedderCopy := make([]*Embedder, len(e.embedder))
	for i := range e.embedder {
		if e.embedder[i] != nil {
			embedderCopy[i] = e.embedder[i].SafeCopy()
		}
	}

	return baseMulOperator{
		rank: e.rank,
		mod:  e.mod,

		ambModLen: e.ambModLen,
		ambMod:    e.ambMod,
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		buf: newMulOperatorBuffer(e.rank, max(1, vec.Max(e.ambModLen))),
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

	buf mulOperatorBuffer
}

// newAnyCyclotomicMulOperator creates a new [anyCyclotomicMulOperator].
func newAnyCyclotomicMulOperator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) anyCyclotomicMulOperator {
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		maxBits := num.Log2(params.CycloOrder()) + 2*num.Log2(mod[i].Value())
		ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
	}

	ambParams := dft.NewCyclicParameters(params.CycloOrder())
	ambMod := dft.MustFindPrevNTTPrimes(ambParams, num.MaxModulusBits, vec.Max(ambModLen))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
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

		reducer: reducer.SafeCopy(),

		buf: newMulOperatorBuffer(params.CycloOrder(), max(1, vec.Max(ambModLen))),
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					clear(op.buf.p0.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])

					copy(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					clear(op.buf.p1.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], op.buf.p1.Coeffs[j])

					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOut.Coeffs[i], op.buf.p0.Coeffs[0], i)
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulAddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					clear(op.buf.p0.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])

					copy(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					clear(op.buf.p1.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], op.buf.p1.Coeffs[j])

					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(op.buf.p0.Coeffs[0][:rank], op.buf.p0.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0][:rank], op.mod[i])
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				vec.MMulSubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
			} else {
				rank := op.params.Rank()
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					clear(op.buf.p0.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])

					copy(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					clear(op.buf.p1.Coeffs[j][rank:])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], op.buf.p1.Coeffs[j])

					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(op.buf.p0.Coeffs[0][:rank], op.buf.p0.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0][:rank], op.mod[i])
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
	ambNTTCopy := make([]dft.Transformer, maxAmbModLen)
	for i := range ambNTTCopy {
		ambNTTCopy[i] = op.ambNTT[i].SafeCopy()
	}

	embedderCopy := make([]*Embedder, len(idx))
	for i := range idx {
		if op.embedder[idx[i]] != nil {
			embedderCopy[i] = op.embedder[idx[i]].SafeCopy()
		}
	}

	return anyCyclotomicMulOperator{
		params: op.params,
		mod:    modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    op.ambMod[:maxAmbModLen],
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: op.reducer.SubReducer(idx...),

		buf: newMulOperatorBuffer(op.params.CycloOrder(), max(1, maxAmbModLen)),
	}
}

func (e *anyCyclotomicMulOperator) safeCopy() anyCyclotomicMulOperator {
	ambNTTCopy := make([]dft.Transformer, len(e.ambNTT))
	for i := range e.ambNTT {
		ambNTTCopy[i] = e.ambNTT[i].SafeCopy()
	}

	embedderCopy := make([]*Embedder, len(e.embedder))
	for i := range e.embedder {
		if e.embedder[i] != nil {
			embedderCopy[i] = e.embedder[i].SafeCopy()
		}
	}

	return anyCyclotomicMulOperator{
		params: e.params,
		mod:    e.mod,

		ambModLen: e.ambModLen,
		ambMod:    e.ambMod,
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: e.reducer.SafeCopy(),

		buf: newMulOperatorBuffer(e.params.CycloOrder(), max(1, vec.Max(e.ambModLen))),
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

	buf mulOperatorBuffer
}

// newReduceMulOperator creates a new [reduceMulOperator].
func newReduceMulOperator(mod []*num.Modulus, modPoly []int64, reducer *Reducer) reduceMulOperator {
	ambParams := dft.NewCyclicParameters(num.NextProdPower(2*len(modPoly)-1, []int{2}))
	ambModLen := make([]int, len(mod))
	ntt := make([]dft.Transformer, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(ambParams, mod[i]) {
			ntt[i] = dft.NewTransformer(ambParams, mod[i])
		}
		maxBits := num.Log2(ambParams.Rank()) + 2*num.Log2(mod[i].Value())
		ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
	}

	ambMod := dft.MustFindPrevNTTPrimes(ambParams, num.MaxModulusBits, vec.Max(ambModLen))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
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

		reducer: reducer.SafeCopy(),

		buf: newMulOperatorBuffer(ambParams.Rank(), max(1, vec.Max(ambModLen))),
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(op.buf.p0.Coeffs[0], e0.Coeffs[i])
				clear(op.buf.p0.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0])

				copy(op.buf.p1.Coeffs[0], e1.Coeffs[i])
				clear(op.buf.p1.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(op.buf.p1.Coeffs[0], op.buf.p1.Coeffs[0])

				vec.MMulTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0], op.buf.p1.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0])
				op.reducer.reduceTo(eOut.Coeffs[i], op.buf.p0.Coeffs[0], i)
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					clear(op.buf.p0.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])

					copy(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					clear(op.buf.p1.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], op.buf.p1.Coeffs[j])

					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(eOut.Coeffs[i], op.buf.p0.Coeffs[0], i)
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(op.buf.p0.Coeffs[0], e0.Coeffs[i])
				clear(op.buf.p0.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0])

				copy(op.buf.p1.Coeffs[0], e1.Coeffs[i])
				clear(op.buf.p1.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(op.buf.p1.Coeffs[0], op.buf.p1.Coeffs[0])

				vec.MMulTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0], op.buf.p1.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0])
				op.reducer.reduceTo(op.buf.p0.Coeffs[0][:op.rank], op.buf.p0.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					clear(op.buf.p0.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])

					copy(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					clear(op.buf.p1.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], op.buf.p1.Coeffs[j])

					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(op.buf.p0.Coeffs[0][:op.rank], op.buf.p0.Coeffs[0], i)
				vec.AddTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0][:op.rank], op.mod[i])
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

		for i := range op.mod {
			if op.ambModLen[i] == 0 {
				copy(op.buf.p0.Coeffs[0], e0.Coeffs[i])
				clear(op.buf.p0.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0])

				copy(op.buf.p1.Coeffs[0], e1.Coeffs[i])
				clear(op.buf.p1.Coeffs[0][op.rank:])
				op.ntt[i].ForwardTo(op.buf.p1.Coeffs[0], op.buf.p1.Coeffs[0])

				vec.MMulTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0], op.buf.p1.Coeffs[0], op.mod[i])
				op.ntt[i].InverseTo(op.buf.p0.Coeffs[0], op.buf.p0.Coeffs[0])
				op.reducer.reduceTo(op.buf.p0.Coeffs[0][:op.rank], op.buf.p0.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0][:op.rank], op.mod[i])
			} else {
				for j := 0; j < op.ambModLen[i]; j++ {
					copy(op.buf.p0.Coeffs[j], e0.Coeffs[i])
					clear(op.buf.p0.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])

					copy(op.buf.p1.Coeffs[j], e1.Coeffs[i])
					clear(op.buf.p1.Coeffs[j][op.rank:])
					op.ambNTT[j].ForwardTo(op.buf.p1.Coeffs[j], op.buf.p1.Coeffs[j])

					vec.MMulTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j], op.buf.p1.Coeffs[j], op.ambMod[j])
					op.ambNTT[j].InverseTo(op.buf.p0.Coeffs[j], op.buf.p0.Coeffs[j])
				}
				op.embedder[i].EmbedVecTo(op.buf.p0.Coeffs[:1], op.buf.p0.Coeffs[:op.ambModLen[i]])
				op.reducer.reduceTo(op.buf.p0.Coeffs[0][:op.rank], op.buf.p0.Coeffs[0], i)
				vec.SubTo(eOut.Coeffs[i], eOut.Coeffs[i], op.buf.p0.Coeffs[0][:op.rank], op.mod[i])
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
	ambModLenCopy := make([]int, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		ambModLenCopy[i] = op.ambModLen[idx[i]]
	}

	maxAmbModLen := vec.Max(ambModLenCopy)
	ambNTTCopy := make([]dft.Transformer, maxAmbModLen)
	for i := range ambNTTCopy {
		ambNTTCopy[i] = op.ambNTT[i].SafeCopy()
	}

	embedderCopy := make([]*Embedder, len(idx))
	for i := range idx {
		if op.embedder[idx[i]] != nil {
			embedderCopy[i] = op.embedder[idx[i]].SafeCopy()
		}
	}

	return reduceMulOperator{
		rank:    op.rank,
		ambRank: op.ambRank,
		mod:     modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    op.ambMod[:maxAmbModLen],
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: op.reducer.SubReducer(idx...),

		buf: newMulOperatorBuffer(op.ambRank, max(1, maxAmbModLen)),
	}
}

func (op *reduceMulOperator) safeCopy() reduceMulOperator {
	ambNTTCopy := make([]dft.Transformer, len(op.ambNTT))
	for i := range op.ambNTT {
		ambNTTCopy[i] = op.ambNTT[i].SafeCopy()
	}

	embedderCopy := make([]*Embedder, len(op.embedder))
	for i := range op.embedder {
		if op.embedder[i] != nil {
			embedderCopy[i] = op.embedder[i].SafeCopy()
		}
	}

	return reduceMulOperator{
		rank:    op.rank,
		ambRank: op.ambRank,
		mod:     op.mod,

		ambModLen: op.ambModLen,
		ambMod:    op.ambMod,
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: op.reducer.SafeCopy(),

		buf: newMulOperatorBuffer(op.ambRank, max(1, vec.Max(op.ambModLen))),
	}
}
