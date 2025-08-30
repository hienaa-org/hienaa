package crt

import (
	"math"
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

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
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyAutEvaluatorAutFixedPow2{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(params.Rank()),
	}
}

func (e *polyAutEvaluatorAutFixedPow2) Aut(p *Poly, idx uint64) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorAutFixedPow2) AutTo(pOut, p *Poly, idx uint64) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrd, rank := uint64(e.params.CycloOrder()), uint64(e.params.Rank())
	cycloOrdMask := cycloOrd - 1

	idx = idx % cycloOrd

	switch {
	case idx%4 != 1:
		panic("AutTo: idx must be 1 mod 4")
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
				if idxIn >= rank {
					idxIn = rank<<1 - idxIn - 1
				}
				idxOut := bits.Reverse64(j) >> revShiftBits
				if idxOut >= rank {
					idxOut = rank<<1 - idxOut - 1
				}
				pOut.Coeffs[i][idxOut] = e.buf.p[idxIn]
			}
		} else {
			clear(e.buf.p)
			for j := uint64(0); j < rank; j++ {
				idxFull := (j * idx) & cycloOrdMask
				switch {
				case idxFull >= 3*rank:
					e.buf.p[rank<<2-idxFull] = p.Coeffs[i][j]
				case idxFull >= 2*rank:
					e.buf.p[idxFull-rank<<1] = num.Neg(p.Coeffs[i][j], e.mod[i])
				case idxFull >= rank:
					e.buf.p[rank<<1-idxFull] = num.Neg(p.Coeffs[i][j], e.mod[i])
				default:
					e.buf.p[idxFull] = p.Coeffs[i][j]
				}
			}
			copy(pOut.Coeffs[i], e.buf.p)
		}
	}

	pOut.isNTT = p.isNTT
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

	isNTTFriendly := make([]bool, len(mod))
	for i := 0; i < len(mod); i++ {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
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

func (e *polyAutEvaluatorAutFixedPrime) Aut(p *Poly, idx uint64) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.AutTo(pOut, p, idx)
	return pOut
}

func (e *polyAutEvaluatorAutFixedPrime) AutTo(pOut, p *Poly, idx uint64) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("AutTo: inputs not consistent")
	}

	cycloOrder, rank := uint64(e.params.CycloOrder()), uint64(e.params.Rank())

	idx = idx % cycloOrder

	if idx == 0 {
		pOut.CopyFrom(p)
		return
	}

	rotIdx := uint64(slices.Index(e.rootPow, idx))
	if rotIdx == math.MaxUint64 {
		rotIdx = uint64(slices.Index(e.rootPowInv, idx))
		if rotIdx == math.MaxUint64 {
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

func (e *polyAutEvaluatorAutFixedPrime) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorAutFixedPow2{
		params:        e.params,
		mod:           e.mod,
		isNTTFriendly: e.isNTTFriendly,

		buf: newPolyAutEvaluatorBuffer(e.params.Rank()),
	}
}
