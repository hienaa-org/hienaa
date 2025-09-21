package crt

import (
	"math"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// polyMulEvaluator is the evaluator for polynomial multiplication.
type polyMulEvaluator interface {
	// Mul returns p0 * p1.
	// Panics when p0 and p1 are not both in NTT form.
	Mul(p0, p1 *Poly) *Poly
	// MulTo computes pOut = p0 * p1.
	// Panics when p0 and p1 are not both in NTT form.
	MulTo(pOut, p0, p1 *Poly)
	// MulAddTo computes pOut += p0 * p1.
	// Panics when p0 and p1 are not both in NTT form.
	MulAddTo(pOut, p0, p1 *Poly)
	// MulSubTo computes pOut -= p0 * p1.
	// Panics when p0 and p1 are not both in NTT form.
	MulSubTo(pOut, p0, p1 *Poly)
	// subEvaluator returns a evaluator for modulus of given indices.
	subEvaluator(idx ...int) polyMulEvaluator
	// SafeCopy returns a thread-safe copy.
	safeCopy() polyMulEvaluator
}

// polyMulEvaluatorBuffer is a buffer for [polyMulEvaluator].
type polyMulEvaluatorBuffer struct {
	p0 [][]uint64
	p1 [][]uint64
}

// newPolyMulEvaluatorBuffer creates a new [polyMulEvaluatorBuffer].
func newPolyMulEvaluatorBuffer(rank, modLen int) polyMulEvaluatorBuffer {
	p0 := make([][]uint64, modLen)
	p1 := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		p0[i] = make([]uint64, rank)
		p1[i] = make([]uint64, rank)
	}

	return polyMulEvaluatorBuffer{
		p0: p0,
		p1: p1,
	}
}

// polyMulEvaluatorNoReduce is a [polyMulEvaluator] for rings that polynomials are automatically reduced.
// This includes power-of-two cyclotomic rings and autfixed rings.
type polyMulEvaluatorNoReduce struct {
	rank int
	mod  []*num.Modulus

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	buf polyMulEvaluatorBuffer
}

// newPolyMulEvaluatorNoReduce creates a new [polyMulEvaluatorNoReduce].
func newPolyMulEvaluatorNoReduce(params dft.RingParameters, mod []*num.Modulus) *polyMulEvaluatorNoReduce {
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		maxBits := 2 * num.Log2(mod[i].Value())
		switch params.RingType() {
		case dft.Cyclotomic:
			maxBits += num.Log2(params.Rank())
		case dft.AutFixed:
			maxBits += num.Log2(params.CycloOrder())
		}
		ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
	}

	ambMod := dft.FindPrevNTTPrimes(params, num.MaxModulusBits, vec.Max(ambModLen))
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambNTT {
		ambNTT[i] = dft.NewTransformer(params, ambMod[i])
	}

	embedder := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			continue
		}
		embedder[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
	}

	return &polyMulEvaluatorNoReduce{
		rank: params.Rank(),
		mod:  mod,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		buf: newPolyMulEvaluatorBuffer(params.Rank(), max(1, vec.Max(ambModLen))),
	}
}

func (e *polyMulEvaluatorNoReduce) Mul(p0, p1 *Poly) *Poly {
	pOut := NewNTTPoly(e.rank, len(e.mod))
	e.MulTo(pOut, p0, p1)
	return pOut
}

func (e *polyMulEvaluatorNoReduce) MulTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			for j := 0; j < e.ambModLen[i]; j++ {
				e.ambNTT[j].ForwardTo(e.buf.p0[j], p0.Coeffs[i])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], p1.Coeffs[i])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(pOut.Coeffs[i:i+1], e.buf.p0[:e.ambModLen[i]])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorNoReduce) MulAddTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !pOut.isNTT || !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulAddTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			for j := 0; j < e.ambModLen[i]; j++ {
				e.ambNTT[j].ForwardTo(e.buf.p0[j], p0.Coeffs[i])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], p1.Coeffs[i])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorNoReduce) MulSubTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !pOut.isNTT || !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulSubTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			for j := 0; j < e.ambModLen[i]; j++ {
				e.ambNTT[j].ForwardTo(e.buf.p0[j], p0.Coeffs[i])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], p1.Coeffs[i])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorNoReduce) subEvaluator(idx ...int) polyMulEvaluator {
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

	return &polyMulEvaluatorNoReduce{
		rank: e.rank,
		mod:  modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    e.ambMod[:maxAmbModLen],
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		buf: newPolyMulEvaluatorBuffer(e.rank, max(1, maxAmbModLen)),
	}
}

func (e *polyMulEvaluatorNoReduce) safeCopy() polyMulEvaluator {
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

	return &polyMulEvaluatorNoReduce{
		rank: e.rank,
		mod:  e.mod,

		ambModLen: e.ambModLen,
		ambMod:    e.ambMod,
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		buf: newPolyMulEvaluatorBuffer(e.rank, max(1, vec.Max(e.ambModLen))),
	}
}

// polyMulEvaluatorCyclotomicNonPow2 is a [polyMulEvaluator] for non power-of-two cyclotomic rings.
type polyMulEvaluatorCyclotomicNonPow2 struct {
	params dft.RingParameters
	mod    []*num.Modulus

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	reducer *CyclotomicReducer

	buf polyMulEvaluatorBuffer
}

// newPolyMulEvaluatorCyclotomicNonPow2 creates a new [polyMulEvaluatorCyclotomicNonPow2].
func newPolyMulEvaluatorCyclotomicNonPow2(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) *polyMulEvaluatorCyclotomicNonPow2 {
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			continue
		}
		maxBits := num.Log2(params.CycloOrder()) + 2*num.Log2(mod[i].Value())
		ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
	}

	ambParams := dft.NewCyclicParameters(params.CycloOrder())
	ambMod := dft.FindPrevNTTPrimes(ambParams, num.MaxModulusBits, vec.Max(ambModLen))
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

	return &polyMulEvaluatorCyclotomicNonPow2{
		params: params,
		mod:    mod,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		reducer: reducer.SafeCopy(),

		buf: newPolyMulEvaluatorBuffer(params.CycloOrder(), max(1, vec.Max(ambModLen))),
	}
}

func (e *polyMulEvaluatorCyclotomicNonPow2) Mul(p0, p1 *Poly) *Poly {
	pOut := NewNTTPoly(e.params.Rank(), len(e.mod))
	e.MulTo(pOut, p0, p1)
	return pOut
}

func (e *polyMulEvaluatorCyclotomicNonPow2) MulTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			rank := e.params.Rank()
			for j := 0; j < e.ambModLen[i]; j++ {
				copy(e.buf.p0[j], p0.Coeffs[i])
				clear(e.buf.p0[j][rank:])
				e.ambNTT[j].ForwardTo(e.buf.p0[j], e.buf.p0[j])

				copy(e.buf.p1[j], p1.Coeffs[i])
				clear(e.buf.p1[j][rank:])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], e.buf.p1[j])

				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			e.reducer.reduceTo(pOut.Coeffs[i], e.buf.p0[0], i)
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorCyclotomicNonPow2) MulAddTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !pOut.isNTT || !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulAddTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			rank := e.params.Rank()
			for j := 0; j < e.ambModLen[i]; j++ {
				copy(e.buf.p0[j], p0.Coeffs[i])
				clear(e.buf.p0[j][rank:])
				e.ambNTT[j].ForwardTo(e.buf.p0[j], e.buf.p0[j])

				copy(e.buf.p1[j], p1.Coeffs[i])
				clear(e.buf.p1[j][rank:])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], e.buf.p1[j])

				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			e.reducer.reduceTo(e.buf.p0[0][:rank], e.buf.p0[0], i)
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:rank], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorCyclotomicNonPow2) MulSubTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !pOut.isNTT || !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulSubTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			rank := e.params.Rank()
			for j := 0; j < e.ambModLen[i]; j++ {
				copy(e.buf.p0[j], p0.Coeffs[i])
				clear(e.buf.p0[j][rank:])
				e.ambNTT[j].ForwardTo(e.buf.p0[j], e.buf.p0[j])

				copy(e.buf.p1[j], p1.Coeffs[i])
				clear(e.buf.p1[j][rank:])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], e.buf.p1[j])

				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			e.reducer.reduceTo(e.buf.p0[0][:rank], e.buf.p0[0], i)
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:rank], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorCyclotomicNonPow2) subEvaluator(idx ...int) polyMulEvaluator {
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

	return &polyMulEvaluatorCyclotomicNonPow2{
		params: e.params,
		mod:    modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    e.ambMod[:maxAmbModLen],
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: e.reducer.SafeCopy(),

		buf: newPolyMulEvaluatorBuffer(e.params.CycloOrder(), max(1, maxAmbModLen)),
	}
}

func (e *polyMulEvaluatorCyclotomicNonPow2) safeCopy() polyMulEvaluator {
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

	return &polyMulEvaluatorCyclotomicNonPow2{
		params: e.params,
		mod:    e.mod,

		ambModLen: e.ambModLen,
		ambMod:    e.ambMod,
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: e.reducer.SafeCopy(),

		buf: newPolyMulEvaluatorBuffer(e.params.CycloOrder(), max(1, vec.Max(e.ambModLen))),
	}
}

// polyMulEvaluatorReduce is a [polyMulEvaluator] for arbitrary modulo rings.
type polyMulEvaluatorReduce struct {
	rank    int
	ambRank int
	mod     []*num.Modulus

	ntt []dft.Transformer

	ambModLen []int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  []*Embedder

	reducer *Reducer

	buf polyMulEvaluatorBuffer
}

// newPolyMulEvaluatorReduce creates a new [polyMulEvaluatorCyclotomicReduce].
func newPolyMulEvaluatorReduce(mod []*num.Modulus, modPoly []int64, reducer *Reducer) *polyMulEvaluatorReduce {
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

	ambMod := dft.FindPrevNTTPrimes(ambParams, num.MaxModulusBits, vec.Max(ambModLen))
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

	return &polyMulEvaluatorReduce{
		rank:    len(modPoly) - 1,
		ambRank: ambParams.Rank(),
		mod:     mod,

		ntt: ntt,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		reducer: reducer.SafeCopy(),

		buf: newPolyMulEvaluatorBuffer(ambParams.Rank(), max(1, vec.Max(ambModLen))),
	}
}

func (e *polyMulEvaluatorReduce) Mul(p0, p1 *Poly) *Poly {
	pOut := NewNTTPoly(e.rank, len(e.mod))
	e.MulTo(pOut, p0, p1)
	return pOut
}

func (e *polyMulEvaluatorReduce) MulTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			copy(e.buf.p1[0], p1.Coeffs[i])
			clear(e.buf.p1[0][e.rank:])
			e.ntt[i].ForwardTo(e.buf.p0[0], e.buf.p0[0])

			copy(e.buf.p0[0], p0.Coeffs[i])
			clear(e.buf.p0[0][e.rank:])
			e.ntt[i].ForwardTo(e.buf.p1[0], e.buf.p1[0])

			vec.MMulTo(e.buf.p0[0], e.buf.p0[0], e.buf.p1[0], e.mod[i])
			e.ntt[i].InverseTo(e.buf.p0[0], e.buf.p0[0])
			e.reducer.reduceTo(pOut.Coeffs[i], e.buf.p0[0], i)
		} else {
			for j := 0; j < e.ambModLen[i]; j++ {
				copy(e.buf.p0[j], p0.Coeffs[i])
				clear(e.buf.p0[j][e.rank:])
				e.ambNTT[j].ForwardTo(e.buf.p0[j], e.buf.p0[j])

				copy(e.buf.p1[j], p1.Coeffs[i])
				clear(e.buf.p1[j][e.rank:])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], e.buf.p1[j])

				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			e.reducer.reduceTo(pOut.Coeffs[i], e.buf.p0[0], i)
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorReduce) MulAddTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulAddTo: inputs not consistent")
	case !pOut.isNTT || !p0.isNTT || !p1.isNTT:
		panic("MulAddTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			copy(e.buf.p1[0], p1.Coeffs[i])
			clear(e.buf.p1[0][e.rank:])
			e.ntt[i].ForwardTo(e.buf.p0[0], e.buf.p0[0])

			copy(e.buf.p0[0], p0.Coeffs[i])
			clear(e.buf.p0[0][e.rank:])
			e.ntt[i].ForwardTo(e.buf.p1[0], e.buf.p1[0])

			vec.MMulTo(e.buf.p0[0], e.buf.p0[0], e.buf.p1[0], e.mod[i])
			e.ntt[i].InverseTo(e.buf.p0[0], e.buf.p0[0])
			e.reducer.reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0], i)
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		} else {
			for j := 0; j < e.ambModLen[i]; j++ {
				copy(e.buf.p0[j], p0.Coeffs[i])
				clear(e.buf.p0[j][e.rank:])
				e.ambNTT[j].ForwardTo(e.buf.p0[j], e.buf.p0[j])

				copy(e.buf.p1[j], p1.Coeffs[i])
				clear(e.buf.p1[j][e.rank:])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], e.buf.p1[j])

				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			e.reducer.reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0], i)
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorReduce) MulSubTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulSubTo: inputs not consistent")
	case !pOut.isNTT || !p0.isNTT || !p1.isNTT:
		panic("MulSubTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			copy(e.buf.p1[0], p1.Coeffs[i])
			clear(e.buf.p1[0][e.rank:])
			e.ntt[i].ForwardTo(e.buf.p0[0], e.buf.p0[0])

			copy(e.buf.p0[0], p0.Coeffs[i])
			clear(e.buf.p0[0][e.rank:])
			e.ntt[i].ForwardTo(e.buf.p1[0], e.buf.p1[0])

			vec.MMulTo(e.buf.p0[0], e.buf.p0[0], e.buf.p1[0], e.mod[i])
			e.ntt[i].InverseTo(e.buf.p0[0], e.buf.p0[0])
			e.reducer.reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0], i)
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		} else {
			for j := 0; j < e.ambModLen[i]; j++ {
				copy(e.buf.p0[j], p0.Coeffs[i])
				clear(e.buf.p0[j][e.rank:])
				e.ambNTT[j].ForwardTo(e.buf.p0[j], e.buf.p0[j])

				copy(e.buf.p1[j], p1.Coeffs[i])
				clear(e.buf.p1[j][e.rank:])
				e.ambNTT[j].ForwardTo(e.buf.p1[j], e.buf.p1[j])

				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseTo(e.buf.p0[j], e.buf.p0[j])
			}
			e.embedder[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:e.ambModLen[i]])
			e.reducer.reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0], i)
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorReduce) subEvaluator(idx ...int) polyMulEvaluator {
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

	return &polyMulEvaluatorReduce{
		rank:    e.rank,
		ambRank: e.ambRank,
		mod:     modCopy,

		ambModLen: ambModLenCopy,
		ambMod:    e.ambMod[:maxAmbModLen],
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: e.reducer.SafeCopy(),

		buf: newPolyMulEvaluatorBuffer(e.ambRank, max(1, maxAmbModLen)),
	}
}

func (e *polyMulEvaluatorReduce) safeCopy() polyMulEvaluator {
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

	return &polyMulEvaluatorReduce{
		rank:    e.rank,
		ambRank: e.ambRank,
		mod:     e.mod,

		ambModLen: e.ambModLen,
		ambMod:    e.ambMod,
		ambNTT:    ambNTTCopy,
		embedder:  embedderCopy,

		reducer: e.reducer.SafeCopy(),

		buf: newPolyMulEvaluatorBuffer(e.ambRank, max(1, vec.Max(e.ambModLen))),
	}
}
