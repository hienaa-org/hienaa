package rns

import (
	"math"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

type reducerAnyModulus struct {
}

type reducerNTTModulus struct {
	// modulus is the modulus.
	modulus *mod.Modulus

	// deg is the degree of the modulo polynomial.
	deg int

	// diffDeg is the degree of the difference between the modulo polynomial and the maximum degree of the input polynomial.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the modulo polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT singleTransformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT singleTransformer

	// modPoly is the polynomial we target to reduce to.
	modPoly []uint64
	// quotientPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is ⌊X^d/modPoly⌋ modulo the modulus, where d is the maximum degree of the input polynomial.
	quotientPoly []uint64

	// buf is the polynomial buffer for the reducer.
	buf reducerNTTBuffer
}

func NewReducerNTTModulus(modPoly []uint64, maxDeg int, modulus *mod.Modulus) *reducerNTTModulus {
	deg := len(modPoly)

	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT singleTransformer
	var quotientPoly []uint64
	var buf reducerNTTBuffer

	degNext = num.NextProdPower(uint64(deg), []uint64{2})
	diffDeg = uint64(maxDeg - deg)
	diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

	degNextParams := RingParameters{0, int(degNext), Cyclic}
	degNextNTT = newCyclicPow235Transformer(degNextParams, modulus)

	diffDegNextParams := RingParameters{0, int(diffDegNext), Cyclic}
	diffDegNextNTT = newCyclicPow235Transformer(diffDegNextParams, modulus)

	modPoly = append(modPoly, make([]uint64, int(degNext)-len(modPoly))...)

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	quotientPoly = quotientPolynomialMod(dividend, modPoly[:deg+1], modulus)
	quotientPoly = append(quotientPoly, make([]uint64, int(diffDegNext)-len(quotientPoly))...)

	degNextNTT.nttInPlace(modPoly)
	diffDegNextNTT.nttInPlace(quotientPoly)

	buf = newReducerNTTBuffer(int(deg), int(diffDegNext), int(degNext))

	return &reducerNTTModulus{
		modulus: modulus,

		deg:         int(deg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly:      modPoly,
		quotientPoly: quotientPoly,

		buf: buf,
	}
}

func (r *reducerNTTModulus) ReduceTo(pOut, pIn []uint64) {
	copy(r.buf.pIn, pIn)

	// Compute pQuo = ⌊pIn/X^deg⌋
	clear(r.buf.pQuo)
	for i := 0; i < r.diffDeg; i++ {
		r.buf.pQuo[i] = r.buf.pIn[r.deg+i]
	}

	// Compute pQuo = pQuo × ⌊X^(deg+diffDeg)/modPoly⌋
	r.diffDegNextNTT.nttInPlace(r.buf.pQuo)
	mod.MMulLazyVecTo(r.buf.pQuo, r.buf.pQuo, r.quotientPoly, r.modulus)
	r.diffDegNextNTT.invNTTInPlace(r.buf.pQuo)

	// Compute pRem = ⌊pQuo/X^diffDeg⌋ % (X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.diffDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pQuo[r.diffDeg+j] = mod.Add(r.buf.pQuo[r.diffDeg+j], r.buf.pQuo[r.diffDeg+i*r.degNext+j], r.modulus)
			r.buf.pQuo[r.diffDeg+i*r.degNext+j] = 0
		}
	}

	clear(r.buf.pRem)
	for i := 0; i < min(r.diffDeg, r.degNext); i++ {
		r.buf.pRem[i] = r.buf.pQuo[r.diffDeg+i]
	}

	// Compute pRem = pRem × quotientPoly (mod X^degNext - 1)
	r.degNextNTT.nttInPlace(r.buf.pRem)
	mod.MMulLazyVecTo(r.buf.pRem, r.buf.pRem, r.modPoly, r.modulus)
	r.degNextNTT.invNTTInPlace(r.buf.pRem)

	// Compute pIn = pIn (mod X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.diffDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pIn[j] = mod.Add(r.buf.pIn[j], r.buf.pIn[i*r.degNext+j], r.modulus)
			r.buf.pIn[i*r.degNext+j] = 0
		}
	}

	// Compute pOut = pIn - pRem
	for i := 0; i < r.deg; i++ {
		pOut[i] = mod.Sub(r.buf.pIn[i], r.buf.pRem[i], r.modulus)
	}
}
