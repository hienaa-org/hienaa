package dft

import (
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclotomicReducer reduces a polynomial modulo a cyclotomic polynomial.
// In other words, it computes a(X) mod \Phi_m(X), where a(X) is at most degree m-1.
// It uses Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
// In this implementation, Q_sp = (X^m-1)/(X^(m/p)-1) for the smallest prime factor p.
type cyclotomicReducer struct {
	params RingParameters
	mod    *num.Modulus

	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isTrivial is true if the reduction is trivial.
	isTrivial bool

	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
	// diffDeg is the degree of the difference between the cyclotomic polynomial and the intermediate reducing polynomial Q_sp.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the cyclotomic polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT *pow235CyclicTransformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT *pow235CyclicTransformer

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly []uint64
	// divPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/\Phi_m(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	divPoly []uint64

	buf reducerBuffer
}

// reducerBuffer is a buffer for [cyclotomicReducerNTTModulus].
type reducerBuffer struct {
	// pIn is a buffer for the input polynomial.
	pIn []uint64
	// pQuo is a buffer for the quotient polynomial.
	pQuo []uint64
	// pRem is a buffer for the remainder polynomial.
	pRem []uint64
}

// newReducerBuffer creates a new [reducerBuffer].
func newReducerBuffer(in, quo, rem int) reducerBuffer {
	return reducerBuffer{
		pIn:  make([]uint64, in),
		pQuo: make([]uint64, quo),
		pRem: make([]uint64, rem),
	}
}

// newCyclotomicReducer creates a new [cyclotomicReducer].
func newCyclotomicReducer(params RingParameters, mod *num.Modulus) *cyclotomicReducer {
	primes, _ := num.Factor(params.cycloOrd)

	var redDeg, leastFactor int
	if params.cycloOrd%2 == 1 {
		leastFactor = primes[0]
		redDeg = params.cycloOrd - params.cycloOrd/leastFactor
	} else {
		leastFactor = primes[1]
		redDeg = params.cycloOrd/2 - (params.cycloOrd/2)/leastFactor
	}

	isTrivial := redDeg == params.rank

	var diffDeg, diffDegNext, degNext int
	var diffDegNextNTT, degNextNTT *pow235CyclicTransformer
	var cycloPoly, divPoly []uint64

	if !isTrivial {
		degNext = num.NextProdPower(params.rank, []int{2})
		diffDeg = redDeg - params.rank
		diffDegNext = num.NextProdPower(2*diffDeg+1, []int{2})

		diffDegNextNTT = newCyclicPow235Transformer(NewCyclicParameters(diffDegNext), mod)
		degNextNTT = newCyclicPow235Transformer(NewCyclicParameters(degNext), mod)

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1

		cycloPoly = make([]uint64, degNext)
		vec.ReduceTo(cycloPoly[:params.rank+1], CyclotomicPolynomial(params.cycloOrd), mod)

		divPoly = quotient(dividend, cycloPoly[:params.rank+1], mod)
		divPoly = append(divPoly, make([]uint64, diffDegNext-len(divPoly))...)

		diffDegNextNTT.ForwardTo(divPoly, divPoly)
		degNextNTT.ForwardTo(cycloPoly, cycloPoly)
	}

	return &cyclotomicReducer{
		params: params,
		mod:    mod,

		leastFac:  leastFactor,
		isTrivial: isTrivial,

		redDeg:      redDeg,
		diffDeg:     diffDeg,
		diffDegNext: diffDegNext,
		degNext:     degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly: cycloPoly,
		divPoly:   divPoly,

		buf: newReducerBuffer(params.cycloOrd, diffDegNext, degNext),
	}
}

func (r *cyclotomicReducer) reduceTo(pOut, p []uint64) {
	cycloOrd, rank := r.params.cycloOrd, r.params.rank

	copy(r.buf.pIn, p)
	if cycloOrd%2 == 1 {
		skip := cycloOrd / r.leastFac

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				r.buf.pIn[i*skip+j] = num.Sub(r.buf.pIn[i*skip+j], r.buf.pIn[cycloOrd-skip+j], r.mod)
			}
			r.buf.pIn[cycloOrd-skip+j] = 0
		}
	} else {
		skip := (cycloOrd / 2) / r.leastFac

		for i := 0; i < cycloOrd/2; i++ {
			r.buf.pIn[i] = num.Sub(r.buf.pIn[i], r.buf.pIn[cycloOrd/2+i], r.mod)
			r.buf.pIn[cycloOrd/2+i] = 0
		}

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				if i%2 == 0 {
					r.buf.pIn[i*skip+j] = num.Sub(r.buf.pIn[i*skip+j], r.buf.pIn[cycloOrd/2-skip+j], r.mod)
				} else {
					r.buf.pIn[i*skip+j] = num.Add(r.buf.pIn[i*skip+j], r.buf.pIn[cycloOrd/2-skip+j], r.mod)
				}
			}
			r.buf.pIn[cycloOrd/2-skip+j] = 0
		}
	}

	if !r.isTrivial {
		// pQuo = floor(pIn / X^deg)
		clear(r.buf.pQuo)
		for j := 0; j < r.diffDeg; j++ {
			r.buf.pQuo[j] = r.buf.pIn[rank+j]
		}

		// pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
		r.diffDegNextNTT.ForwardTo(r.buf.pQuo, r.buf.pQuo)
		vec.MMulLazyTo(r.buf.pQuo, r.buf.pQuo, r.divPoly, r.mod)
		r.diffDegNextNTT.InverseTo(r.buf.pQuo, r.buf.pQuo)

		// pRem = floor(pQuo / X^diffDeg) % (X^degNext - 1)
		for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.diffDeg {
					break
				}
				r.buf.pQuo[r.diffDeg+j] = num.Add(r.buf.pQuo[r.diffDeg+j], r.buf.pQuo[r.diffDeg+i*r.degNext+j], r.mod)
				r.buf.pQuo[r.diffDeg+i*r.degNext+j] = 0
			}
		}

		clear(r.buf.pRem)
		for j := 0; j < min(r.diffDeg, r.degNext); j++ {
			r.buf.pRem[j] = r.buf.pQuo[r.diffDeg+j]
		}

		// pRem = pRem * cycloPoly % (X^degNext - 1)
		r.degNextNTT.ForwardTo(r.buf.pRem, r.buf.pRem)
		vec.MMulLazyTo(r.buf.pRem, r.buf.pRem, r.cycloPoly, r.mod)
		r.degNextNTT.InverseTo(r.buf.pRem, r.buf.pRem)

		// pIn = pIn % X^degNext - 1
		for i := 1; i <= num.DivCeil(r.redDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[j] = num.Add(r.buf.pIn[j], r.buf.pIn[i*r.degNext+j], r.mod)
				r.buf.pIn[i*r.degNext+j] = 0
			}
		}

		// pOut = pIn - pRem
		for i := 0; i < rank; i++ {
			pOut[i] = num.Sub(r.buf.pIn[i], r.buf.pRem[i], r.mod)
		}
	} else {
		copy(pOut, r.buf.pIn[:rank])
	}
}

func (r *cyclotomicReducer) SafeCopy() *cyclotomicReducer {
	return &cyclotomicReducer{
		params: r.params,
		mod:    r.mod,

		leastFac:  r.leastFac,
		isTrivial: r.isTrivial,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: r.diffDegNextNTT,
		degNextNTT:     r.degNextNTT,

		cycloPoly: r.cycloPoly,
		divPoly:   r.divPoly,

		buf: newReducerBuffer(r.params.cycloOrd, r.diffDegNext, r.degNext),
	}
}

// CyclotomicPolynomial computes the cyclotomic polynomial of the given cyclotomic order.
func CyclotomicPolynomial(cycloOrd int) []int64 {
	primes, _ := num.Factor(cycloOrd)

	var isEven bool
	if primes[0] == 2 {
		isEven = true
		primes = primes[1:]
	} else {
		isEven = false
	}

	pOut := make([]int64, cycloOrd+1)
	pBuf0 := make([]int64, cycloOrd+1)
	pBuf1 := make([]int64, cycloOrd+1)
	pOut[0], pOut[1] = -1, 1

	skip := int(cycloOrd)

	currDeg, prevDeg := 1, 1
	for _, prime := range primes {
		copy(pBuf0, pOut)
		clear(pBuf1)
		clear(pOut)

		for i := 0; i <= prevDeg; i++ {
			pBuf1[i*prime] = pBuf0[i]
		}

		currDeg = prevDeg*prime - prevDeg

		for i := 0; i <= (prime-1)*prevDeg; i++ {
			if pBuf1[prevDeg*prime-i] != 0 {
				pOut[currDeg-i] = pBuf1[prevDeg*prime-i] / pBuf0[prevDeg]

				for j := 0; j <= prevDeg; j++ {
					pBuf1[prevDeg*prime-i-j] -= pBuf0[prevDeg-j] * pOut[currDeg-i]
				}
			}
		}

		prevDeg = currDeg
		skip /= prime
	}

	if isEven {
		for i := 1; i <= prevDeg; i += 2 {
			pOut[i] = -pOut[i]
		}

		skip >>= 1
	}

	if skip > 1 {
		copy(pBuf0, pOut)
		clear(pOut)
		for i := 0; i <= prevDeg; i++ {
			pOut[i*skip] = pBuf0[i]
		}
	}

	return pOut[:num.Totient(cycloOrd)+1]
}

// quotient computes the quotient of two polynomials modulo a modulus.
func quotient(p0, p1 []uint64, mod *num.Modulus) []uint64 {
	if len(p0) < len(p1) {
		panic("dividend must be longer than divisor")
	} else if num.GCD(mod.Value(), p1[len(p1)-1]) != 1 {
		panic("divisor must be coprime with modulus")
	}

	quo := make([]uint64, len(p0)-len(p1)+1)
	rem := make([]uint64, len(p0))
	copy(rem, p0)

	lcInv := num.Inv(p1[len(p1)-1], mod)
	for i := 0; i <= len(p0)-len(p1); i++ {
		if rem[len(rem)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(rem[len(rem)-i-1], lcInv, mod)
			vec.MulSubScalarTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], mod)
		}
	}

	return quo
}
