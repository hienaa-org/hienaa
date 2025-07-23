package rns

import (
	"math"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

type cyclotomicReducerAnyModulus struct {
}

// cyclotomicReducerNTTModulus reduces a polynomial modulo a cyclotomic polynomial.
// In other words, it computes a(X) mod \Phi_m(X), where a(X) is at most degree m-1.
// It uses Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
// In this implementation, Q_sp = (Xᵐ-1)/(X^(m/p)-1) for the smallest prime factor p.
type cyclotomicReducerNTTModulus struct {
	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isPrimePow is true if the cyclotomic polynomial is a prime power.
	isPrimePow bool

	// params is the ring parameters.
	params RingParameters
	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
	// modulus is the modulus.
	modulus *mod.Modulus

	// diffDeg is the degree of the difference between the cyclotomic polynomial and the intermediate reducing polynomial Q_sp.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the cyclotomic polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT singleTransformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT singleTransformer

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly []uint64
	// quotientPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is ⌊X^d_qs/\Phi_m(X)⌋ modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quotientPoly []uint64

	// buf is the polynomial buffer for the reducer.
	buf reducerNTTBuffer
}

func NewCyclotomicReducerNTTModulus(params RingParameters, modulus *mod.Modulus) *cyclotomicReducerNTTModulus {
	cycloDeg, deg := uint64(params.CycloDegree()), uint64(params.Degree())

	factors := num.Factor(cycloDeg)
	leastFactor := cycloDeg
	for key := range factors {
		if key < leastFactor {
			leastFactor = key
		}
	}
	redDeg := cycloDeg - cycloDeg/leastFactor

	var isPrimePower bool
	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT singleTransformer
	var cycloPoly, quotientPoly []uint64
	var buf reducerNTTBuffer

	if redDeg == deg {
		isPrimePower = true
		diffDeg, diffDegNext, degNext = 0, 0, 0
		diffDegNextNTT, degNextNTT = newTrivialTransformer(modulus), newTrivialTransformer(modulus)
		cycloPoly, quotientPoly = []uint64{}, []uint64{}
		buf = newReducerNTTBuffer(0, 0, 0)
	} else {
		degNext = num.NextProdPower(uint64(deg), []uint64{2})
		diffDeg = redDeg - deg
		diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

		degNextParams := RingParameters{0, int(degNext), Cyclic}
		degNextNTT = newCyclicPow235Transformer(degNextParams, modulus)

		diffDegNextParams := RingParameters{0, int(diffDegNext), Cyclic}
		diffDegNextNTT = newCyclicPow235Transformer(diffDegNextParams, modulus)

		cycloPoly = make([]uint64, degNext)
		cycloPolyInt := computeCyclotomicPolynomial(cycloDeg)
		for i := range cycloPolyInt {
			cycloPoly[i] = reduceInt(cycloPolyInt[i], modulus)
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quotientPoly = append(quotientPolynomialMod(dividend, cycloPoly[:deg+1], modulus), make([]uint64, int(diffDegNext)-len(quotientPoly))...)

		degNextNTT.nttInPlace(cycloPoly)
		diffDegNextNTT.nttInPlace(quotientPoly)

		buf = newReducerNTTBuffer(int(cycloDeg), int(diffDegNext), int(degNext))
	}

	return &cyclotomicReducerNTTModulus{
		leastFac:   int(leastFactor),
		isPrimePow: isPrimePower,

		params:  params,
		redDeg:  int(redDeg),
		modulus: modulus,

		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly:    cycloPoly,
		quotientPoly: quotientPoly,

		buf: buf,
	}
}

func (r *cyclotomicReducerNTTModulus) ReduceTo(pOut, pIn []uint64) {
	copy(r.buf.pIn, pIn)

	leastFac := r.leastFac
	cycloDeg := r.params.CycloDegree()
	deg := r.params.Degree()
	skip := cycloDeg / leastFac

	for j := 0; j < skip; j++ {
		for i := 0; i < leastFac-1; i++ {
			r.buf.pIn[i*skip+j] = mod.Sub(r.buf.pIn[i*skip+j], r.buf.pIn[cycloDeg-skip+j], r.modulus)
		}
		r.buf.pIn[cycloDeg-skip+j] = 0
	}

	if !r.isPrimePow {
		// Compute pQuo = ⌊pIn/X^deg⌋
		clear(r.buf.pQuo)
		for i := 0; i < r.diffDeg; i++ {
			r.buf.pQuo[i] = r.buf.pIn[deg+i]
		}

		// Compute pQuo = pQuo × ⌊X^(deg+diffDeg)/\Phi_m(X)⌋
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
		mod.MMulLazyVecTo(r.buf.pRem, r.buf.pRem, r.cycloPoly, r.modulus)
		r.degNextNTT.invNTTInPlace(r.buf.pRem)

		// Compute pIn = pIn (mod X^degNext - 1)
		for i := 1; i <= int(math.Ceil(float64(r.redDeg)/float64(r.degNext))); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[j] = mod.Add(r.buf.pIn[j], r.buf.pIn[i*r.degNext+j], r.modulus)
				r.buf.pIn[i*r.degNext+j] = 0
			}
		}

		// Compute pOut = pIn - pRem
		for i := 0; i < r.params.Degree(); i++ {
			pOut[i] = mod.Sub(r.buf.pIn[i], r.buf.pRem[i], r.modulus)
		}
	} else {
		copy(pOut, r.buf.pIn[:r.params.Degree()])
	}
}
