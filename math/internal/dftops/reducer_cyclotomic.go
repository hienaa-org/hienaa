package dftops

import (
	"math"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// CyclotomicReducerNTTModulus reduces a polynomial modulo a cyclotomic polynomial.
// In other words, it computes a(X) mod \Phi_m(X), where a(X) is at most degree m-1.
// It uses Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
// In this implementation, Q_sp = (X^m-1)/(X^(m/p)-1) for the smallest prime factor p.
type CyclotomicReducerNTTModulus struct {
	cycloOrd int
	rank     int
	mod      *num.Modulus

	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isPrimePow is true if the cyclotomic polynomial is a prime power.
	isPrimePow bool

	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
	// diffDeg is the degree of the difference between the cyclotomic polynomial and the intermediate reducing polynomial Q_sp.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the cyclotomic polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT *CyclicPow2Transformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT *CyclicPow2Transformer

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly []uint64
	// quoPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/\Phi_m(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quoPoly []uint64

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

// newCyclotomicReducerNTTModulus creates a new [cyclotomicReducerNTTModulus].
func NewCyclotomicReducerNTTModulus(cycloOrd, rank int, mod *num.Modulus) *CyclotomicReducerNTTModulus {
	primes, _ := num.Factor(uint64(cycloOrd))

	var redDeg, leastFactor int
	if cycloOrd&1 == 1 {
		leastFactor = cycloOrd
		for _, p := range primes {
			if int(p) < leastFactor {
				leastFactor = int(p)
			}
		}
		redDeg = cycloOrd - cycloOrd/leastFactor
	} else {
		leastFactor = cycloOrd
		for _, p := range primes {
			if int(p) < leastFactor && p > 2 {
				leastFactor = int(p)
			}
		}
		redDeg = cycloOrd/2 - cycloOrd/2/leastFactor
	}

	isPrimePower := redDeg == rank

	var diffDeg, diffDegNext, degNext int
	var diffDegNextNTT, degNextNTT *CyclicPow2Transformer
	var cycloPoly, quoPoly []uint64

	if !isPrimePower {
		degNext = int(num.NextProdPower(uint64(rank), []uint64{2}))
		diffDeg = redDeg - rank
		diffDegNext = int(num.NextProdPower(2*uint64(diffDeg)+1, []uint64{2}))

		degNextNTT = NewCyclicPow2Transformer(int(degNext), mod)
		diffDegNextNTT := NewCyclicPow2Transformer(int(diffDegNext), mod)

		cycloPoly = make([]uint64, degNext)
		cycloPolySigned := CyclotomicPolynomial(int(cycloOrd))
		for i := range cycloPolySigned {
			cycloPoly[i] = reduceInt(cycloPolySigned[i], mod)
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quoPoly = append(quotient(dividend, cycloPoly[:rank+1], mod), make([]uint64, int(diffDegNext)-len(quoPoly))...)

		degNextNTT.ForwardInPlace(cycloPoly)
		diffDegNextNTT.ForwardInPlace(quoPoly)
	}

	return &CyclotomicReducerNTTModulus{
		cycloOrd: cycloOrd,
		rank:     rank,
		mod:      mod,

		leastFac:   int(leastFactor),
		isPrimePow: isPrimePower,

		redDeg:      int(redDeg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly: cycloPoly,
		quoPoly:   quoPoly,

		buf: newReducerBuffer(int(cycloOrd), int(diffDegNext), int(degNext)),
	}
}

// newReducerBuffer creates a new [reducerBuffer].
func newReducerBuffer(in, quo, rem int) reducerBuffer {
	return reducerBuffer{
		pIn:  make([]uint64, in),
		pQuo: make([]uint64, quo),
		pRem: make([]uint64, rem),
	}
}

func (r *CyclotomicReducerNTTModulus) ReduceTo(pOut, p []uint64) {
	copy(r.buf.pIn, p)

	cycloOrd := int(r.cycloOrd)
	rank := int(r.rank)

	if cycloOrd&1 == 1 {
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
				if i&1 == 0 {
					r.buf.pIn[i*skip+j] = num.Sub(r.buf.pIn[i*skip+j], r.buf.pIn[cycloOrd/2-skip+j], r.mod)
				} else {
					r.buf.pIn[i*skip+j] = num.Add(r.buf.pIn[i*skip+j], r.buf.pIn[cycloOrd/2-skip+j], r.mod)
				}
			}
			r.buf.pIn[cycloOrd-skip+j] = 0
		}
	}

	if !r.isPrimePow {
		// Compute pQuo = floor(pIn/X^deg)
		clear(r.buf.pQuo)
		for i := 0; i < r.diffDeg; i++ {
			r.buf.pQuo[i] = r.buf.pIn[rank+i]
		}

		// Compute pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
		r.diffDegNextNTT.ForwardInPlace(r.buf.pQuo)
		vec.MMulLazyTo(r.buf.pQuo, r.buf.pQuo, r.quoPoly, r.mod)
		r.diffDegNextNTT.InverseInPlace(r.buf.pQuo)

		// Compute pRem = floor(pQuo/X^diffDeg) % (X^degNext - 1)
		for i := 1; i <= int(math.Ceil(float64(r.diffDeg)/float64(r.degNext))); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.diffDeg {
					break
				}
				r.buf.pQuo[r.diffDeg+j] = num.Add(r.buf.pQuo[r.diffDeg+j], r.buf.pQuo[r.diffDeg+i*r.degNext+j], r.mod)
				r.buf.pQuo[r.diffDeg+i*r.degNext+j] = 0
			}
		}

		clear(r.buf.pRem)
		for i := 0; i < min(r.diffDeg, r.degNext); i++ {
			r.buf.pRem[i] = r.buf.pQuo[r.diffDeg+i]
		}

		// Compute pRem = pRem * quoPoly (mod X^degNext - 1)
		r.degNextNTT.ForwardInPlace(r.buf.pRem)
		vec.MMulLazyTo(r.buf.pRem, r.buf.pRem, r.cycloPoly, r.mod)
		r.degNextNTT.InverseInPlace(r.buf.pRem)

		// Compute pIn = pIn (mod X^degNext - 1)
		for i := 1; i <= int(math.Ceil(float64(r.redDeg)/float64(r.degNext))); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[j] = num.Add(r.buf.pIn[j], r.buf.pIn[i*r.degNext+j], r.mod)
				r.buf.pIn[i*r.degNext+j] = 0
			}
		}

		// Compute pOut = pIn - pRem
		for i := 0; i < rank; i++ {
			pOut[i] = num.Sub(r.buf.pIn[i], r.buf.pRem[i], r.mod)
		}
	} else {
		copy(pOut, r.buf.pIn[:rank])
	}
}

func (r *CyclotomicReducerNTTModulus) SafeCopy() *CyclotomicReducerNTTModulus {
	var buf reducerBuffer
	if !r.isPrimePow {
		buf = newReducerBuffer(int(r.cycloOrd), r.diffDegNext, r.degNext)
	} else {
		buf = newReducerBuffer(int(r.cycloOrd), 0, 0)
	}

	return &CyclotomicReducerNTTModulus{
		cycloOrd: r.cycloOrd,
		rank:     r.rank,
		mod:      r.mod,

		leastFac:   r.leastFac,
		isPrimePow: r.isPrimePow,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: r.diffDegNextNTT,
		degNextNTT:     r.degNextNTT,

		cycloPoly: r.cycloPoly,
		quoPoly:   r.quoPoly,

		buf: buf,
	}
}
