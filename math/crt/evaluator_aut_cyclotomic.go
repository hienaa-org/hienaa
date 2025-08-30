package crt

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

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
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyAutEvaluatorCyclotomicPow2{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

func (e *polyAutEvaluatorCyclotomicPow2) Aut(p *Poly, idx uint64) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorCyclotomicPow2) AutTo(pOut, p *Poly, idx uint64) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := uint64(e.params.CycloOrder()), uint64(e.params.Rank())
	cycloOrdMask := cycloOrd - 1

	idx = idx % cycloOrd

	switch {
	case idx%2 != 1:
		panic("AutTo: idx must be 1 mod 2")
	case idx == 0:
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			copy(e.buf.p, p.Coeffs[i])
			revShiftBits := 64 - uint64(num.Log2(rank)+1)
			for j := uint64(0); j < rank; j++ {
				idxIn := bits.Reverse64(((j<<1+1)*idx-1)>>1) >> revShiftBits
				idxOut := bits.Reverse64(j) >> revShiftBits
				pOut.Coeffs[i][idxOut] = e.buf.p[idxIn]
			}
		} else {
			clear(e.buf.p)
			for j := uint64(0); j < rank; j++ {
				idxFull := (j * idx) & cycloOrdMask
				switch {
				case idxFull >= rank:
					e.buf.p[idxFull-rank] = num.Neg(p.Coeffs[i][j], e.mod[i])
				default:
					e.buf.p[idxFull] = p.Coeffs[i][j]
				}
			}
			copy(pOut.Coeffs[i], e.buf.p)
		}
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
	mod           []*num.Modulus
	cycloOrdMod   *num.Modulus
	isNTTFriendly []bool

	reducer []reducer

	// primeExpMods is the prime power factors of the cyclotomic order.
	primeExpMods []*num.Modulus
	// rootExps is the generators modulo prime power factors of the cyclotomic order.
	rootExps [][]uint64
	// dims is the dimension of the hypercube structure.
	dims []uint64

	buf polyAutEvaluatorBuffer
}

// newPolyAutEvaluatorCyclotomicNonPow2 creates a new [polyAutEvaluatorCyclotomicNonPow2].
func newPolyAutEvaluatorCyclotomicNonPow2(params dft.RingParameters, mod []*num.Modulus, reducers []reducer) *polyAutEvaluatorCyclotomicNonPow2 {
	isNTTFriendly := make([]bool, len(mod))
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	reducerCopy := make([]reducer, len(reducers))
	for i := range reducerCopy {
		reducerCopy[i] = reducers[i].safeCopy()
	}

	primes, exps := num.Factor(params.CycloOrder())
	pExpMods := make([]*num.Modulus, len(primes))
	rootExps := make([][]uint64, len(primes))
	dims := make([]uint64, len(primes))
	for i := range rootExps {
		pExp := 1
		for j := 0; j < int(exps[i]); j++ {
			pExp *= primes[i]
		}
		pExpMods[i] = num.NewModulus(pExp)
		dims[i] = uint64(pExp - pExp/primes[i])
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
			dims = append([]uint64{dims[0] / 2, 2}, dims[1:]...)
		}
	}

	return &polyAutEvaluatorCyclotomicNonPow2{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
		cycloOrdMod:   num.NewModulus(params.CycloOrder()),

		reducer: reducerCopy,

		primeExpMods: pExpMods,
		rootExps:     rootExps,
		dims:         dims,

		buf: newPolyAutEvaluatorBuffer(params.CycloOrder()),
	}
}

func (e *polyAutEvaluatorCyclotomicNonPow2) Aut(p *Poly, idx uint64) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorCyclotomicNonPow2) AutTo(pOut, p *Poly, idx uint64) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := uint64(e.params.CycloOrder()), uint64(e.params.Rank())

	idx = idx % cycloOrd

	switch {
	case num.GCD(idx, cycloOrd) != 1:
		panic("AutTo: idx must be coprime to cycloOrd")
	case idx == 0:
		pOut.CopyFrom(p)
		return
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			idxDigits := make([]uint64, len(e.primeExpMods))

			for cnt := 0; cnt < len(idxDigits); {
				if e.primeExpMods[cnt].Value()%8 == 0 {
					idxRed := num.Reduce(idx, e.primeExpMods[cnt])
					if idx%4 == 1 {
						idxDigits[cnt+1] = 0
					} else {
						idxDigits[cnt+1] = 1
						idxRed = num.Neg(idxRed, e.primeExpMods[cnt])
					}
					idxDigits[cnt] = uint64(slices.Index(e.rootExps[cnt], idxRed))
					cnt += 2
				} else {
					idxRed := num.Reduce(idx, e.primeExpMods[cnt])
					idxDigits[cnt] = uint64(slices.Index(e.rootExps[cnt], idxRed))
					cnt += 1
				}
			}

			copy(e.buf.p, p.Coeffs[i])
			idxOutDigits := make([]uint64, len(idxDigits))
			for j := uint64(0); j < rank; j++ {
				idxIn := j
				for k := 0; k < len(idxOutDigits); k++ {
					idxOutDigits[k] = (idxIn - idxDigits[k] + e.dims[k]) % e.dims[k]
					idxIn /= e.dims[k]
				}
				idxOut := idxOutDigits[len(idxOutDigits)-1]
				for k := len(idxOutDigits) - 2; k >= 0; k-- {
					idxOut *= e.dims[k]
					idxOut += idxOutDigits[k]
				}
				pOut.Coeffs[i][idxOut] = e.buf.p[j]
			}
		} else {
			clear(e.buf.p)
			for j := uint64(0); j < rank; j++ {
				idxOut := num.Mul(j, idx, e.cycloOrdMod)
				e.buf.p[idxOut] = p.Coeffs[i][j]
			}
			e.reducer[i].reduceTo(pOut.Coeffs[i], e.buf.p)
		}
	}
}

func (e *polyAutEvaluatorCyclotomicNonPow2) safeCopy() polyAutEvaluator {
	reducerCopy := make([]reducer, len(e.reducer))
	for i := range reducerCopy {
		reducerCopy[i] = e.reducer[i].safeCopy()
	}

	return &polyAutEvaluatorCyclotomicNonPow2{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		reducer: reducerCopy,

		primeExpMods: e.primeExpMods,
		rootExps:     e.rootExps,
		dims:         e.dims,

		buf: newPolyAutEvaluatorBuffer(e.params.CycloOrder()),
	}
}
