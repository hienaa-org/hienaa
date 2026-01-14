package crt

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// polyAutEvaluator is the evaluator for automorphism.
type polyAutEvaluator interface {
	// CanAut returns whether the given automorphism index is valid.
	CanAut(idx int) bool
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

// pow2CyclotomicPolyAutEvaluator is a [polyAutEvaluator] for power-of-two cyclotomic ring.
type pow2CyclotomicPolyAutEvaluator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	buf polyAutEvaluatorBuffer
}

// newPow2CyclotomicPolyAutEvaluator creates a new [pow2CyclotomicPolyAutEvaluator].
func newPow2CyclotomicPolyAutEvaluator(params dft.RingParameters, mod []*num.Modulus) pow2CyclotomicPolyAutEvaluator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return pow2CyclotomicPolyAutEvaluator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

// CanAut returns whether the given automorphism index is valid.
func (e *pow2CyclotomicPolyAutEvaluator) CanAut(idx int) bool {
	cycloOrd := e.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return idx%2 == 1
}

// Aut returns aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *pow2CyclotomicPolyAutEvaluator) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *pow2CyclotomicPolyAutEvaluator) AutTo(pOut, p *Poly, idx int) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("AutTo: inputs not consistent")
	case !e.CanAut(idx):
		panic("AutTo: idx not supported")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	if idx == 1 {
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
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

	pOut.IsNTT = p.IsNTT
}

func (e *pow2CyclotomicPolyAutEvaluator) subEvaluator(idx ...int) pow2CyclotomicPolyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return pow2CyclotomicPolyAutEvaluator{
		params:        e.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

func (e *pow2CyclotomicPolyAutEvaluator) safeCopy() pow2CyclotomicPolyAutEvaluator {
	return pow2CyclotomicPolyAutEvaluator{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

// anyCyclotomicPolyAutEvaluator is a [polyAutEvaluator] for arbitrary order cyclotomic ring.
type anyCyclotomicPolyAutEvaluator struct {
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

// newAnyCyclotomicPolyAutEvaluator creates a new [polyAutEvaluatorCyclotomicNonPow2].
func newAnyCyclotomicPolyAutEvaluator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) anyCyclotomicPolyAutEvaluator {
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

	return anyCyclotomicPolyAutEvaluator{
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

// CanAut returns whether the given automorphism index is valid.
func (e *anyCyclotomicPolyAutEvaluator) CanAut(idx int) bool {
	cycloOrd := e.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return num.GCD(idx, cycloOrd) == 1
}

// Aut returns aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *anyCyclotomicPolyAutEvaluator) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *anyCyclotomicPolyAutEvaluator) AutTo(pOut, p *Poly, idx int) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("AutTo: inputs not consistent")
	case !e.CanAut(idx):
		panic("AutTo: idx not supported")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	if idx == 1 {
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
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

	pOut.IsNTT = p.IsNTT
}

func (e *anyCyclotomicPolyAutEvaluator) subEvaluator(idx ...int) anyCyclotomicPolyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return anyCyclotomicPolyAutEvaluator{
		params:        e.params,
		cycloOrdMod:   e.cycloOrdMod,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		reducer: e.reducer.SubReducer(idx...),

		primeExpMods: e.primeExpMods,
		rootExps:     e.rootExps,
		dims:         e.dims,

		buf: newPolyAutEvaluatorBuffer(e.params.CycloOrder()),
	}
}

func (e *anyCyclotomicPolyAutEvaluator) safeCopy() anyCyclotomicPolyAutEvaluator {
	return anyCyclotomicPolyAutEvaluator{
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

// pow2AutFixedPolyAutEvaluator is a [polyAutEvaluator] for power-of-two conjugate invariant ring.
type pow2AutFixedPolyAutEvaluator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	buf polyAutEvaluatorBuffer
}

// newPow2AutFixedPolyAutEvaluator creates a new [pow2AutFixedPolyAutEvaluator].
func newPow2AutFixedPolyAutEvaluator(params dft.RingParameters, mod []*num.Modulus) pow2AutFixedPolyAutEvaluator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return pow2AutFixedPolyAutEvaluator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

// CanAut returns whether the given automorphism index is valid.
func (e *pow2AutFixedPolyAutEvaluator) CanAut(idx int) bool {
	cycloOrd := e.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return idx%4 == 1
}

// Aut returns aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *pow2AutFixedPolyAutEvaluator) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *pow2AutFixedPolyAutEvaluator) AutTo(pOut, p *Poly, idx int) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("AutTo: inputs not consistent")
	case !e.CanAut(idx):
		panic("AutTo: idx not supported")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	if idx == 1 {
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
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

	pOut.IsNTT = p.IsNTT
}

func (e *pow2AutFixedPolyAutEvaluator) subEvaluator(idx ...int) pow2AutFixedPolyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return pow2AutFixedPolyAutEvaluator{
		params:        e.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

func (e *pow2AutFixedPolyAutEvaluator) safeCopy() pow2AutFixedPolyAutEvaluator {
	return pow2AutFixedPolyAutEvaluator{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

// primeAutFixedPolyAutEvaluator is a [polyAutEvaluator] for prime order autfixed ring.
type primeAutFixedPolyAutEvaluator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	// rootPow are the power of the generator modulo the cyclotomic order.
	rootPow []uint64
	// rootPowInv are the powers of the inverse of the generator modulo the cyclotomic order.
	rootPowInv []uint64

	buf polyAutEvaluatorBuffer
}

// newPrimeAutFixedPolyAutEvaluator creates a new [primeAutFixedPolyAutEvaluator].
func newPrimeAutFixedPolyAutEvaluator(params dft.RingParameters, mod []*num.Modulus) primeAutFixedPolyAutEvaluator {
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

	return primeAutFixedPolyAutEvaluator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		rootPow:    rootPow,
		rootPowInv: rootPowInv,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

// CanAut returns whether the given automorphism index is valid.
func (e *primeAutFixedPolyAutEvaluator) CanAut(idx int) bool {
	cycloOrd := e.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	rotIdx := slices.Index(e.rootPow, uint64(idx))
	if rotIdx == -1 {
		rotIdx = slices.Index(e.rootPowInv, uint64(idx))
		if rotIdx == -1 {
			return false
		}
	}
	return true
}

// Aut returns aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *primeAutFixedPolyAutEvaluator) Aut(p *Poly, idx int) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *primeAutFixedPolyAutEvaluator) AutTo(pOut, p *Poly, idx int) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("AutTo: inputs not consistent")
	case !e.CanAut(idx):
		panic("AutTo: idx not valid")
	}

	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := e.params.CycloOrder(), e.params.Rank()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd

	rotIdx := slices.Index(e.rootPow, uint64(idx))
	if rotIdx == -1 {
		rotIdx = slices.Index(e.rootPowInv, uint64(idx))
		rotIdx = (rank - rotIdx) % rank
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
			copy(e.buf.p[:rank-rotIdx], p.Coeffs[i][rotIdx:])
			copy(e.buf.p[rank-rotIdx:], p.Coeffs[i][:rotIdx])
			copy(pOut.Coeffs[i], e.buf.p)
		} else {
			copy(e.buf.p[:rotIdx], p.Coeffs[i][rank-rotIdx:])
			copy(e.buf.p[rotIdx:], p.Coeffs[i][:rank-rotIdx])
			copy(pOut.Coeffs[i], e.buf.p)
		}
	}

	pOut.IsNTT = p.IsNTT
}

func (e *primeAutFixedPolyAutEvaluator) subEvaluator(idx ...int) primeAutFixedPolyAutEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return primeAutFixedPolyAutEvaluator{
		params:        e.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		rootPow:    e.rootPow,
		rootPowInv: e.rootPowInv,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

func (e *primeAutFixedPolyAutEvaluator) safeCopy() primeAutFixedPolyAutEvaluator {
	return primeAutFixedPolyAutEvaluator{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		rootPow:    e.rootPow,
		rootPowInv: e.rootPowInv,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}

// noAutPolyAutEvaluator always panics.
// Used for cyclic and other rings.
type noAutPolyAutEvaluator struct{}

// CanAut returns whether the given automorphism index is valid.
func (e *noAutPolyAutEvaluator) CanAut(idx int) bool {
	return false
}

// Aut returns aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *noAutPolyAutEvaluator) Aut(p *Poly, idx int) *Poly {
	panic("Aut: automorphism not supported in this ring")
}

// AutTo computes pOut = aut_idx(p).
// Panics when automorphism is invalid.
// Notable cases include:
//
//   - In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
//   - In any other rings, automorphism is not supported and it always panics.
func (e *noAutPolyAutEvaluator) AutTo(pOut, p *Poly, idx int) {
	panic("AutTo: automorphism not supported in this ring")
}
