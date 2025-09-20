package crt

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// polyAutEvaluator is the evaluator for automorphism.
type polyAutEvaluator interface {
	// Aut returns aut_idx(p).
	// Panics when automorphism is invalid.
	// Notable cases include:
	//
	//	- In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
	//	- In any other rings, automorphism is not supported and it always panics.
	Aut(p *Poly, idx int) *Poly
	// AutTo computes pOut = aut_idx(p).
	// Panics when automorphism is invalid.
	// Notable cases include:
	//
	//	- In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
	//	- In any other rings, automorphism is not supported and it always panics.
	AutTo(pOut, p *Poly, idx int)
	// subEvaluator returns a evaluator for modulus of given indices.
	subEvaluator(idx ...int) polyAutEvaluator
	// safeCopy returns a thread-safe copy.
	safeCopy() polyAutEvaluator
}

// polyAutEvaluatorBuffer is a buffer for [polyAutEvaluator].
type polyAutEvaluatorBuffer struct {
	p []uint64
}

// newPolyAutEvaluatorBuffer creates a new [polyAutEvaluatorBuffer].
func newPolyAutEvaluatorBuffer(rank int) polyAutEvaluatorBuffer {
	return polyAutEvaluatorBuffer{
		p: make([]uint64, rank),
	}
}

// polyAutEvaluatorCyclotomicPow2 is a [polyAutEvaluator] for power-of-two cyclotomic rings.
type polyAutEvaluatorCyclotomicPow2 struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	buf polyAutEvaluatorBuffer
}

// newPolyAutEvaluatorCyclotomicPow2 creates a new [polyAutEvaluatorCyclotomicPow2].
func newPolyAutEvaluatorCyclotomicPow2(params dft.RingParameters, mod []*num.Modulus) *polyAutEvaluatorCyclotomicPow2 {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyAutEvaluatorCyclotomicPow2{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

func (e *polyAutEvaluatorCyclotomicPow2) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorCyclotomicPow2) AutTo(pOut, p *Poly, idx int) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()

	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	switch {
	case idx%2 != 1:
		panic("AutTo: idx must be odd")
	case idx == 1:
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			copy(e.buf.p, p.Coeffs[i])
			revShiftBits := 64 - int(num.Log2(rank))
			for j := 0; j < rank; j++ {
				jOut := ((2*j + 1) * idx) & (cycloOrd - 1)
				idxIn := int(bits.Reverse64((uint64(jOut)-1)/2) >> revShiftBits)
				idxOut := int(bits.Reverse64(uint64(j)) >> revShiftBits)
				pOut.Coeffs[i][idxOut] = e.buf.p[idxIn]
			}
		} else {
			clear(e.buf.p)
			for j := 0; j < rank; j++ {
				idxOut := (j * idx) & (cycloOrd - 1)
				switch {
				case idxOut >= rank:
					e.buf.p[idxOut-rank] = num.Neg(p.Coeffs[i][j], e.mod[i])
				default:
					e.buf.p[idxOut] = p.Coeffs[i][j]
				}
			}
			copy(pOut.Coeffs[i], e.buf.p)
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyAutEvaluatorCyclotomicPow2) subEvaluator(idx ...int) polyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return &polyAutEvaluatorCyclotomicPow2{
		params:        e.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

func (e *polyAutEvaluatorCyclotomicPow2) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorCyclotomicPow2{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

// polyAutEvaluatorCyclotomicNonPow2 is a [polyAutEvaluator] for non power-of-two cyclotomic rings.
type polyAutEvaluatorCyclotomicNonPow2 struct {
	params        dft.RingParameters
	cycloOrdMod   *num.Modulus
	mod           []*num.Modulus
	isNTTFriendly []bool

	reducer *CyclotomicReducer

	// primeExpMods is the prime power factors of the cyclotomic order.
	primeExpMods []*num.Modulus
	// rootExps is the generators modulo prime power factors of the cyclotomic order.
	rootExps [][]uint64
	// dims is the dimension of the hypercube structure.
	dims []int

	buf polyAutEvaluatorBuffer
}

// newPolyAutEvaluatorCyclotomicNonPow2 creates a new [polyAutEvaluatorCyclotomicNonPow2].
func newPolyAutEvaluatorCyclotomicNonPow2(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) *polyAutEvaluatorCyclotomicNonPow2 {
	isNTTFriendly := make([]bool, len(mod))
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	primes, exps := num.Factor(params.CycloOrder())
	pExpMods := make([]*num.Modulus, len(primes))
	rootExps := make([][]uint64, len(primes))
	dims := make([]int, len(primes))
	for i := range rootExps {
		pExp := 1
		for j := 0; j < int(exps[i]); j++ {
			pExp *= primes[i]
		}
		pExpMods[i] = num.NewModulus(pExp)
		dims[i] = pExp - pExp/primes[i]
		rootExps[i] = make([]uint64, dims[i])
		root := uint64(1)
		if pExp != 2 {
			root = num.GeneratorsWithFactors(pExpMods[i], []uint64{uint64(primes[i])}, []uint64{uint64(exps[i])})[0]
		}

		if primes[i] == 2 && exps[i] > 2 {
			root = 5
		}

		rootExps[i][0] = 1
		for j := 1; j < len(rootExps[i]); j++ {
			rootExps[i][j] = num.Mul(root, rootExps[i][j-1], pExpMods[i])
		}
	}

	if primes[0] == 2 {
		if exps[0] == 1 {
			rootExps = rootExps[1:]
			pExpMods = pExpMods[1:]
			dims = dims[1:]
		} else if exps[0] > 2 {
			rootExps = append([][]uint64{rootExps[0][:len(rootExps[0])/2], {1, pExpMods[0].Value() - 1}}, rootExps[1:]...)
			pExpMods = append([]*num.Modulus{pExpMods[0], pExpMods[0]}, pExpMods[1:]...)
			dims = append([]int{dims[0] / 2, 2}, dims[1:]...)
		}
	}

	return &polyAutEvaluatorCyclotomicNonPow2{
		params:        params,
		cycloOrdMod:   num.NewModulus(params.CycloOrder()),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		reducer: reducer.SafeCopy(),

		primeExpMods: pExpMods,
		rootExps:     rootExps,
		dims:         dims,

		buf: newPolyAutEvaluatorBuffer(params.CycloOrder()),
	}
}

func (e *polyAutEvaluatorCyclotomicNonPow2) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorCyclotomicNonPow2) AutTo(pOut, p *Poly, idx int) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()

	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	switch {
	case num.GCD(idx, cycloOrd) != 1:
		panic("AutTo: idx must be coprime with cyclotomic order")
	case idx == 1:
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			idxDigits := make([]int, len(e.primeExpMods))
			for cnt := 0; cnt < len(idxDigits); {
				idxRed := num.Reduce(uint64(idx), e.primeExpMods[cnt])
				if e.primeExpMods[cnt].Value()%8 == 0 {
					if idx%4 == 1 {
						idxDigits[cnt+1] = 0
					} else {
						idxDigits[cnt+1] = 1
						idxRed = num.Neg(idxRed, e.primeExpMods[cnt])
					}
					idxDigits[cnt] = slices.Index(e.rootExps[cnt], idxRed)
					cnt += 2
				} else {
					idxDigits[cnt] = slices.Index(e.rootExps[cnt], idxRed)
					cnt += 1
				}
			}

			copy(e.buf.p, p.Coeffs[i])

			idxOutDigits := make([]int, len(e.primeExpMods))
			for j := 0; j < rank; j++ {
				idxIn := j
				for k := 0; k < len(idxOutDigits); k++ {
					idxOutDigits[k] = (idxIn - idxDigits[k] + e.dims[k]) % e.dims[k]
					idxIn /= e.dims[k]
				}
				idxOut := idxOutDigits[len(idxDigits)-1]
				for k := len(idxOutDigits) - 2; k >= 0; k-- {
					idxOut *= e.dims[k]
					idxOut += idxOutDigits[k]
				}
				pOut.Coeffs[i][idxOut] = e.buf.p[j]
			}
		} else {
			clear(e.buf.p)
			for j := 0; j < rank; j++ {
				idxOut := num.Mul(uint64(j), uint64(idx), e.cycloOrdMod)
				e.buf.p[idxOut] = p.Coeffs[i][j]
			}
			e.reducer.reduceTo(pOut.Coeffs[i], e.buf.p, i)
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyAutEvaluatorCyclotomicNonPow2) subEvaluator(idx ...int) polyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return &polyAutEvaluatorCyclotomicNonPow2{
		params:        e.params,
		cycloOrdMod:   e.cycloOrdMod,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		reducer: e.reducer.SafeCopy(),

		primeExpMods: e.primeExpMods,
		rootExps:     e.rootExps,
		dims:         e.dims,

		buf: newPolyAutEvaluatorBuffer(e.params.CycloOrder()),
	}
}

func (e *polyAutEvaluatorCyclotomicNonPow2) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorCyclotomicNonPow2{
		params:        e.params,
		cycloOrdMod:   e.cycloOrdMod,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		reducer: e.reducer.SafeCopy(),

		primeExpMods: e.primeExpMods,
		rootExps:     e.rootExps,
		dims:         e.dims,

		buf: newPolyAutEvaluatorBuffer(e.params.CycloOrder()),
	}
}

// polyAutEvaluatorAutFixedPow2 is a [polyAutEvaluator] for power-of-two autfixed rings.
type polyAutEvaluatorAutFixedPow2 struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	buf polyAutEvaluatorBuffer
}

// newPolyAutEvaluatorAutFixedPow2 creates a new [polyAutEvaluatorAutFixedPow2].
func newPolyAutEvaluatorAutFixedPow2(params dft.RingParameters, mod []*num.Modulus) *polyAutEvaluatorAutFixedPow2 {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyAutEvaluatorAutFixedPow2{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

func (e *polyAutEvaluatorAutFixedPow2) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorAutFixedPow2) AutTo(pOut, p *Poly, idx int) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()

	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	switch {
	case idx%4 != 1:
		panic("AutTo: idx must be 1 mod 4")
	case idx == 1:
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			copy(e.buf.p, p.Coeffs[i])
			revShiftBits := 64 - int(num.Log2(rank)+1)
			for j := 0; j < rank; j++ {
				jOut := ((2*j + 1) * idx) & (cycloOrd - 1)
				idxIn := int(bits.Reverse64((uint64(jOut)-1)/2) >> revShiftBits)
				if idxIn >= rank {
					idxIn = 2*rank - 1 - idxIn
				}
				idxOut := int(bits.Reverse64(uint64(j)) >> revShiftBits)
				if idxOut >= rank {
					idxOut = 2*rank - 1 - idxOut
				}
				pOut.Coeffs[i][idxOut] = e.buf.p[idxIn]
			}
		} else {
			clear(e.buf.p)
			for j := 0; j < rank; j++ {
				idxOut := (j * idx) & (cycloOrd - 1)
				switch {
				case idxOut >= 3*rank:
					e.buf.p[4*rank-idxOut] = p.Coeffs[i][j]
				case idxOut >= 2*rank:
					e.buf.p[idxOut-2*rank] = num.Neg(p.Coeffs[i][j], e.mod[i])
				case idxOut >= rank:
					e.buf.p[2*rank-idxOut] = num.Neg(p.Coeffs[i][j], e.mod[i])
				default:
					e.buf.p[idxOut] = p.Coeffs[i][j]
				}
			}
			copy(pOut.Coeffs[i], e.buf.p)
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyAutEvaluatorAutFixedPow2) subEvaluator(idx ...int) polyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return &polyAutEvaluatorAutFixedPow2{
		params:        e.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

func (e *polyAutEvaluatorAutFixedPow2) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorAutFixedPow2{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

// polyAutEvaluatorAutFixedPrime is a [polyAutEvaluator] for prime-order autfixed rings.
type polyAutEvaluatorAutFixedPrime struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	// rootPow are the power of the generator modulo the cyclotomic order.
	rootPow []uint64
	// rootPowInv are the powers of the inverse of the generator modulo the cyclotomic order.
	rootPowInv []uint64

	buf polyAutEvaluatorBuffer
}

// newPolyAutEvaluatorAutFixedPrime creates a new [polyAutEvaluatorAutFixedPrime].
func newPolyAutEvaluatorAutFixedPrime(params dft.RingParameters, mod []*num.Modulus) *polyAutEvaluatorAutFixedPrime {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	cycloOrdMod := num.NewModulus(params.CycloOrder())
	root := num.Generators(cycloOrdMod)[0]
	rootInv := num.Inv(root, cycloOrdMod)

	rootPow := make([]uint64, params.Rank())
	rootPowInv := make([]uint64, params.Rank())
	rootPow[0] = 1
	rootPowInv[0] = 1
	for i := 1; i < params.Rank(); i++ {
		rootPow[i] = num.Mul(rootPow[i-1], root, cycloOrdMod)
		rootPowInv[i] = num.Mul(rootPowInv[i-1], rootInv, cycloOrdMod)
	}

	return &polyAutEvaluatorAutFixedPrime{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		rootPow:    rootPow,
		rootPowInv: rootPowInv,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

func (e *polyAutEvaluatorAutFixedPrime) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorAutFixedPrime) AutTo(pOut, p *Poly, idx int) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()

	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	rotIdx := slices.Index(e.rootPow, uint64(idx))
	if rotIdx == -1 {
		rotIdx = slices.Index(e.rootPowInv, uint64(idx))
		if rotIdx == -1 {
			panic("AutTo: idx not valid")
		}
		rotIdx = (rank - rotIdx) % rank
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			copy(e.buf.p[:rank-rotIdx], p.Coeffs[i][rotIdx:])
			copy(e.buf.p[rank-rotIdx:], p.Coeffs[i][:rotIdx])
			copy(pOut.Coeffs[i], e.buf.p)
		} else {
			copy(e.buf.p[:rotIdx], p.Coeffs[i][rank-rotIdx:])
			copy(e.buf.p[rotIdx:], p.Coeffs[i][:rank-rotIdx])
			copy(pOut.Coeffs[i], e.buf.p)
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyAutEvaluatorAutFixedPrime) subEvaluator(idx ...int) polyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return &polyAutEvaluatorAutFixedPrime{
		params:        e.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		rootPow:    e.rootPow,
		rootPowInv: e.rootPowInv,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

func (e *polyAutEvaluatorAutFixedPrime) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorAutFixedPrime{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		rootPow:    e.rootPow,
		rootPowInv: e.rootPowInv,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

// polyAutEvaluatorPanic always panics.
// Used for [dft.Cyclic] rings.
type polyAutEvaluatorPanic struct{}

func (e *polyAutEvaluatorPanic) Aut(p *Poly, idx int) *Poly {
	panic("Aut: automorphism not supported in this ring")
}

func (e *polyAutEvaluatorPanic) AutTo(pOut, p *Poly, idx int) {
	panic("AutTo: automorphism not supported in this ring")
}

func (e *polyAutEvaluatorPanic) subEvaluator(idx ...int) polyAutEvaluator {
	return &polyAutEvaluatorPanic{}
}

func (e *polyAutEvaluatorPanic) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorPanic{}
}
