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
	for i := range p0 {
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
	params dft.RingParameters
	mod    []*num.Modulus

	ambNTT           []dft.Transformer
	ambMod           []*num.Modulus
	ambModLen        []int
	embedderToAmbMod []*Embedder
	embedderToMod    []*Embedder

	buf polyMulEvaluatorBuffer
}

// newPolyMulEvaluatorNoReduce creates a new [polyMulEvaluatorNoReduce].
func NewPolyMulEvaluatorNoReduce(params dft.RingParameters, mod []*num.Modulus) *polyMulEvaluatorNoReduce {
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

	embedderToAmbMod := make([]*Embedder, len(mod))
	embedderToMod := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] > 0 {
			embedderToAmbMod[i] = NewEmbedder(ambMod[:ambModLen[i]], []*num.Modulus{mod[i]})
			embedderToMod[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
		}
	}

	return &polyMulEvaluatorNoReduce{
		params: params,
		mod:    mod,

		ambNTT:           ambNTT,
		ambMod:           ambMod,
		ambModLen:        ambModLen,
		embedderToAmbMod: embedderToAmbMod,
		embedderToMod:    embedderToMod,

		buf: newPolyMulEvaluatorBuffer(params.Rank(), len(ambMod)),
	}
}

func (e *polyMulEvaluatorNoReduce) Mul(p0, p1 *Poly) *Poly {
	pOut := NewNTTPoly(e.params.Rank(), len(e.mod))
	e.MulTo(pOut, p0, p1)
	return pOut
}

func (e *polyMulEvaluatorNoReduce) MulTo(pOut, p0, p1 *Poly) {
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
			ambModLen := e.ambModLen[i]
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], p0.Coeffs)
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], p1.Coeffs)
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(pOut.Coeffs, e.buf.p0[:ambModLen])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorNoReduce) MulAddTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT || !pOut.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulAddTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			ambModLen := e.ambModLen[i]
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], p0.Coeffs)
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], p1.Coeffs)
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorNoReduce) MulSubTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT || !pOut.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulSubTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			ambModLen := e.ambModLen[i]
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], p0.Coeffs)
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], p1.Coeffs)
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorNoReduce) safeCopy() polyMulEvaluator {
	ambNTTCopy := make([]dft.Transformer, len(e.ambNTT))
	for i := range ambNTTCopy {
		ambNTTCopy[i] = e.ambNTT[i].SafeCopy()
	}

	embedderToAmbModCopy := make([]*Embedder, len(e.mod))
	embedderToModCopy := make([]*Embedder, len(e.mod))
	for i := range e.mod {
		embedderToAmbModCopy[i] = e.embedderToAmbMod[i].SafeCopy()
		embedderToModCopy[i] = e.embedderToMod[i].SafeCopy()
	}

	return &polyMulEvaluatorNoReduce{
		params: e.params,
		mod:    e.mod,

		ambNTT:           ambNTTCopy,
		ambMod:           e.ambMod,
		ambModLen:        e.ambModLen,
		embedderToAmbMod: embedderToAmbModCopy,
		embedderToMod:    embedderToModCopy,

		buf: newPolyMulEvaluatorBuffer(e.params.Rank(), len(e.ambMod)),
	}
}

// polyMulEvaluatorCyclotomicNonPow2 is a [polyMulEvaluator] for non power-of-two cyclotomic rings.
type polyMulEvaluatorCyclotomicNonPow2 struct {
	params dft.RingParameters
	mod    []*num.Modulus

	reducer []reducer

	ambNTT           []dft.Transformer
	ambMod           []*num.Modulus
	ambModLen        []int
	embedderToAmbMod []*Embedder
	embedderToMod    []*Embedder

	buf polyMulEvaluatorBuffer
}

// newPolyMulEvaluatorCyclotomicNonPow2 creates a new [polyMulEvaluatorCyclotomicNonPow2].
func newPolyMulEvaluatorCyclotomicNonPow2(params dft.RingParameters, mod []*num.Modulus, reducers []reducer) *polyMulEvaluatorCyclotomicNonPow2 {
	reducerCopy := make([]reducer, len(reducers))
	for i := range reducerCopy {
		reducerCopy[i] = reducers[i].safeCopy()
	}

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
	for i := range ambNTT {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	embedderToAmbMod := make([]*Embedder, len(mod))
	embedderToMod := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] > 0 {
			embedderToAmbMod[i] = NewEmbedder(ambMod[:ambModLen[i]], []*num.Modulus{mod[i]})
			embedderToMod[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
		}
	}

	return &polyMulEvaluatorCyclotomicNonPow2{
		params: params,
		mod:    mod,

		reducer: reducerCopy,

		ambNTT:           ambNTT,
		ambMod:           ambMod,
		ambModLen:        ambModLen,
		embedderToAmbMod: embedderToAmbMod,
		embedderToMod:    embedderToMod,

		buf: newPolyMulEvaluatorBuffer(params.CycloOrder(), len(ambMod)),
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
			ambModLen := e.ambModLen[i]

			copy(e.buf.p0[0], p0.Coeffs[i])
			clear(e.buf.p0[0][len(p0.Coeffs[i]):])
			copy(e.buf.p1[0], p1.Coeffs[i])
			clear(e.buf.p1[0][len(p1.Coeffs[i]):])

			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], e.buf.p0[:1])
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], e.buf.p1[:1])
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j+1])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[0:1], e.buf.p0[1:1+ambModLen])
			e.reducer[i].reduceTo(pOut.Coeffs[i], e.buf.p0[0])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorCyclotomicNonPow2) MulAddTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT || !pOut.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			ambModLen := e.ambModLen[i]

			copy(e.buf.p0[0], p0.Coeffs[i])
			clear(e.buf.p0[0][len(p0.Coeffs[i]):])
			copy(e.buf.p1[0], p1.Coeffs[i])
			clear(e.buf.p1[0][len(p1.Coeffs[i]):])

			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], e.buf.p0[:1])
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], e.buf.p1[:1])
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			e.reducer[i].reduceTo(e.buf.p0[0][:e.params.Rank()], e.buf.p0[0])
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.params.Rank()], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorCyclotomicNonPow2) MulSubTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT || !pOut.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		if e.ambModLen[i] == 0 {
			vec.MMulSubTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
		} else {
			ambModLen := e.ambModLen[i]

			copy(e.buf.p0[0], p0.Coeffs[i])
			clear(e.buf.p0[0][e.params.Rank():])
			copy(e.buf.p1[0], p1.Coeffs[i])
			clear(e.buf.p1[0][e.params.Rank():])

			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], e.buf.p0[:1])
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], e.buf.p1[:1])
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			e.reducer[i].reduceTo(e.buf.p0[0][:e.params.Rank()], e.buf.p0[0])
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.params.Rank()], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorCyclotomicNonPow2) safeCopy() polyMulEvaluator {
	reducerCopy := make([]reducer, len(e.reducer))
	for i := range e.reducer {
		reducerCopy[i] = e.reducer[i].safeCopy()
	}

	ambNTTCopy := make([]dft.Transformer, len(e.ambNTT))
	for i := range ambNTTCopy {
		ambNTTCopy[i] = e.ambNTT[i].SafeCopy()
	}

	embedderToAmbModCopy := make([]*Embedder, len(e.mod))
	embedderToModCopy := make([]*Embedder, len(e.mod))
	for i := range e.mod {
		embedderToAmbModCopy[i] = e.embedderToAmbMod[i].SafeCopy()
		embedderToModCopy[i] = e.embedderToMod[i].SafeCopy()
	}

	return &polyMulEvaluatorCyclotomicNonPow2{
		params: e.params,
		mod:    e.mod,

		reducer: reducerCopy,

		ambNTT:           ambNTTCopy,
		ambMod:           e.ambMod,
		ambModLen:        e.ambModLen,
		embedderToAmbMod: embedderToAmbModCopy,
		embedderToMod:    embedderToModCopy,

		buf: newPolyMulEvaluatorBuffer(e.params.CycloOrder(), len(e.ambMod)),
	}
}

// polyMulEvaluatorReduce is a [polyMulEvaluator] for arbitrary modulo rings.
type polyMulEvaluatorReduce struct {
	rank int
	mod  []*num.Modulus

	ntt     []dft.Transformer
	reducer []reducer

	ambNTT           []dft.Transformer
	ambMod           []*num.Modulus
	ambModLen        []int
	embedderToAmbMod []*Embedder
	embedderToMod    []*Embedder

	buf polyMulEvaluatorBuffer
}

// newPolyMulEvaluatorReduce creates a new [polyMulEvaluatorCyclotomicReduce].
func NewPolyMulEvaluatorReduce(mod []*num.Modulus, modPoly []int64, reducers []reducer) *polyMulEvaluatorReduce {
	reducerCopy := make([]reducer, len(reducers))
	for i := range reducerCopy {
		reducerCopy[i] = reducers[i].safeCopy()
	}

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
	for i := range ambNTT {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	embedderToAmbMod := make([]*Embedder, len(mod))
	embedderToMod := make([]*Embedder, len(mod))
	for i := range mod {
		if ambModLen[i] > 0 {
			embedderToAmbMod[i] = NewEmbedder(ambMod[:ambModLen[i]], []*num.Modulus{mod[i]})
			embedderToMod[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
		}
	}

	return &polyMulEvaluatorReduce{
		rank: len(modPoly),
		mod:  mod,

		ntt:     ntt,
		reducer: reducerCopy,

		ambNTT:           ambNTT,
		ambMod:           ambMod,
		ambModLen:        ambModLen,
		embedderToAmbMod: embedderToAmbMod,
		embedderToMod:    embedderToMod,

		buf: newPolyMulEvaluatorBuffer(ambParams.Rank(), max(len(ambMod), 1)),
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
		copy(e.buf.p0[0], p0.Coeffs[i])
		clear(e.buf.p0[0][e.rank:])
		copy(e.buf.p1[0], p1.Coeffs[i])
		clear(e.buf.p1[0][e.rank:])

		if e.ambModLen[i] == 0 {
			e.ntt[i].ForwardInPlace(e.buf.p0[0])
			e.ntt[i].ForwardInPlace(e.buf.p1[0])
			vec.MMulTo(e.buf.p0[0], e.buf.p0[0], e.buf.p1[0], e.mod[i])
			e.ntt[i].InverseInPlace(e.buf.p0[0])
			e.reducer[i].reduceTo(pOut.Coeffs[i], e.buf.p0[0])
		} else {
			ambModLen := e.ambModLen[i]
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], e.buf.p0[:1])
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], e.buf.p1[:1])
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			e.reducer[i].reduceTo(pOut.Coeffs[i], e.buf.p0[0])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorReduce) MulAddTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT || !pOut.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		copy(e.buf.p0[0], p0.Coeffs[i])
		clear(e.buf.p0[0][e.rank:])
		copy(e.buf.p1[0], p1.Coeffs[i])
		clear(e.buf.p1[0][e.rank:])

		if e.ambModLen[i] == 0 {
			e.ntt[i].ForwardInPlace(e.buf.p0[0])
			e.ntt[i].ForwardInPlace(e.buf.p1[0])
			vec.MMulTo(e.buf.p0[0], e.buf.p0[0], e.buf.p1[0], e.mod[i])
			e.ntt[i].InverseInPlace(e.buf.p0[0])
			e.reducer[i].reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0])
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		} else {
			ambModLen := e.ambModLen[i]
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], e.buf.p0[:1])
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], e.buf.p1[:1])
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			e.reducer[i].reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0])
			vec.AddTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorReduce) MulSubTo(pOut, p0, p1 *Poly) {
	switch {
	case !isTernaryToOperable(e.rank, len(e.mod), pOut, p0, p1):
		panic("MulTo: inputs not consistent")
	case !p0.isNTT || !p1.isNTT || !pOut.isNTT:
		panic("MulTo: not in NTT form")
	}

	for i := range e.mod {
		copy(e.buf.p0[0], p0.Coeffs[i])
		clear(e.buf.p0[0][e.rank:])
		copy(e.buf.p1[0], p1.Coeffs[i])
		clear(e.buf.p1[0][e.rank:])

		if e.ambModLen[i] == 0 {
			e.ntt[i].ForwardInPlace(e.buf.p0[0])
			e.ntt[i].ForwardInPlace(e.buf.p1[0])
			vec.MMulTo(e.buf.p0[0], e.buf.p0[0], e.buf.p1[0], e.mod[i])
			e.ntt[i].InverseInPlace(e.buf.p0[0])
			e.reducer[i].reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0])
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		} else {
			ambModLen := e.ambModLen[i]
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p0[:ambModLen], e.buf.p0[:1])
			e.embedderToAmbMod[i].EmbedVecTo(e.buf.p1[:ambModLen], e.buf.p1[:1])
			for j := 0; j < ambModLen; j++ {
				e.ambNTT[j].ForwardInPlace(e.buf.p0[j])
				e.ambNTT[j].ForwardInPlace(e.buf.p1[j])
				vec.MMulTo(e.buf.p0[j], e.buf.p0[j], e.buf.p1[j], e.ambMod[j])
				e.ambNTT[j].InverseInPlace(e.buf.p0[j])
			}
			e.embedderToMod[i].EmbedVecTo(e.buf.p0[:1], e.buf.p0[:ambModLen])
			e.reducer[i].reduceTo(e.buf.p0[0][:e.rank], e.buf.p0[0])
			vec.SubTo(pOut.Coeffs[i], pOut.Coeffs[i], e.buf.p0[0][:e.rank], e.mod[i])
		}
	}

	pOut.isNTT = true
}

func (e *polyMulEvaluatorReduce) safeCopy() polyMulEvaluator {
	reducerCopy := make([]reducer, len(e.reducer))
	for i := range e.reducer {
		reducerCopy[i] = e.reducer[i].safeCopy()
	}

	nttCopy := make([]dft.Transformer, len(e.ntt))
	for i := range nttCopy {
		if e.ntt[i] != nil {
			nttCopy[i] = e.ntt[i].SafeCopy()
		}
	}

	ambNTTCopy := make([]dft.Transformer, len(e.ambNTT))
	for i := range ambNTTCopy {
		ambNTTCopy[i] = e.ambNTT[i].SafeCopy()
	}

	embedderToAmbModCopy := make([]*Embedder, len(e.mod))
	embedderToModCopy := make([]*Embedder, len(e.mod))
	for i := range e.mod {
		embedderToAmbModCopy[i] = e.embedderToAmbMod[i].SafeCopy()
		embedderToModCopy[i] = e.embedderToMod[i].SafeCopy()
	}

	return &polyMulEvaluatorReduce{
		rank: e.rank,
		mod:  e.mod,

		ntt:     nttCopy,
		reducer: reducerCopy,

		ambNTT:           ambNTTCopy,
		ambMod:           e.ambMod,
		ambModLen:        e.ambModLen,
		embedderToAmbMod: embedderToAmbModCopy,
		embedderToMod:    embedderToModCopy,

		buf: newPolyMulEvaluatorBuffer(e.rank, len(e.ambMod)),
	}
}
