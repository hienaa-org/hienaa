package dft

import (
	"github.com/hienaa-org/hienaa/internal/pool"
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

	pool *pool.Pool[*[]uint64]
}

// newCyclotomicReducer creates a new [cyclotomicReducer].
func newCyclotomicReducer(params RingParameters, mod *num.Modulus) *cyclotomicReducer {
	primes, _ := num.Factor(params.cycloIdx)

	var redDeg, leastFactor int
	if params.cycloIdx%2 == 1 {
		leastFactor = primes[0]
		redDeg = params.cycloIdx - params.cycloIdx/leastFactor
	} else {
		leastFactor = primes[1]
		redDeg = params.cycloIdx/2 - (params.cycloIdx/2)/leastFactor
	}

	isTrivial := redDeg == params.rank

	var diffDeg, diffDegNext, degNext int
	var diffDegNextNTT, degNextNTT *pow235CyclicTransformer
	var cycloPoly, divPoly []uint64

	if !isTrivial {
		degNext = num.NextProdPower((params.rank)+1, []int{2})
		diffDeg = redDeg - params.rank
		diffDegNext = num.NextProdPower((2*diffDeg+1)+1, []int{2})

		diffDegNextNTT = newPow235CyclicTransformer(NewCyclicParameters(diffDegNext), mod)
		degNextNTT = newPow235CyclicTransformer(NewCyclicParameters(degNext), mod)

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1

		cycloPoly = make([]uint64, degNext)
		vec.ReduceTo(cycloPoly[:params.rank+1], params.modPoly, mod)

		divPoly = quotient(dividend, cycloPoly[:params.rank+1], mod)
		divPoly = append(divPoly, make([]uint64, diffDegNext-len(divPoly))...)

		diffDegNextNTT.forwardTo(divPoly, divPoly)
		degNextNTT.forwardTo(cycloPoly, cycloPoly)
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

		pool: pool.NewPool(func() *[]uint64 {
			p := make([]uint64, max(params.cycloIdx, diffDegNext, degNext))
			return &p
		}),
	}
}

func (r *cyclotomicReducer) reduceTo(pOut, p []uint64) {
	cycloIdx, rank := r.params.cycloIdx, r.params.rank

	pInPtr := r.pool.Get()
	pIn := (*pInPtr)[:cycloIdx]
	defer r.pool.Put(pInPtr)

	copy(pIn, p)
	if cycloIdx%2 == 1 {
		skip := cycloIdx / r.leastFac

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				pIn[i*skip+j] = num.Sub(pIn[i*skip+j], pIn[cycloIdx-skip+j], r.mod)
			}
			pIn[cycloIdx-skip+j] = 0
		}
	} else {
		skip := (cycloIdx / 2) / r.leastFac

		for i := 0; i < cycloIdx/2; i++ {
			pIn[i] = num.Sub(pIn[i], pIn[cycloIdx/2+i], r.mod)
			pIn[cycloIdx/2+i] = 0
		}

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				if i%2 == 0 {
					pIn[i*skip+j] = num.Sub(pIn[i*skip+j], pIn[cycloIdx/2-skip+j], r.mod)
				} else {
					pIn[i*skip+j] = num.Add(pIn[i*skip+j], pIn[cycloIdx/2-skip+j], r.mod)
				}
			}
			pIn[cycloIdx/2-skip+j] = 0
		}
	}

	if !r.isTrivial {
		pQuoPtr := r.pool.Get()
		pQuo := (*pQuoPtr)[:r.diffDegNext]
		defer r.pool.Put(pQuoPtr)

		// pQuo = floor(pIn / X^deg)
		clear(pQuo)
		for j := 0; j < r.diffDeg; j++ {
			pQuo[j] = pIn[rank+j]
		}

		// pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
		r.diffDegNextNTT.forwardTo(pQuo, pQuo)
		vec.MMulLazyTo(pQuo, pQuo, r.divPoly, r.mod)
		r.diffDegNextNTT.inverseTo(pQuo, pQuo)

		pRemPtr := r.pool.Get()
		pRem := (*pRemPtr)[:r.degNext]
		defer r.pool.Put(pRemPtr)

		// pRem = floor(pQuo / X^diffDeg) % (X^degNext - 1)
		for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.diffDeg {
					break
				}
				pQuo[r.diffDeg+j] = num.Add(pQuo[r.diffDeg+j], pQuo[r.diffDeg+i*r.degNext+j], r.mod)
				pQuo[r.diffDeg+i*r.degNext+j] = 0
			}
		}

		clear(pRem)
		for j := 0; j < min(r.diffDeg, r.degNext); j++ {
			pRem[j] = pQuo[r.diffDeg+j]
		}

		// pRem = pRem * cycloPoly % (X^degNext - 1)
		r.degNextNTT.forwardTo(pRem, pRem)
		vec.MMulLazyTo(pRem, pRem, r.cycloPoly, r.mod)
		r.degNextNTT.inverseTo(pRem, pRem)

		// pIn = pIn % X^degNext - 1
		for i := 1; i <= num.DivCeil(r.redDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				pIn[j] = num.Add(pIn[j], pIn[i*r.degNext+j], r.mod)
				pIn[i*r.degNext+j] = 0
			}
		}

		// pOut = pIn - pRem
		for i := 0; i < rank; i++ {
			pOut[i] = num.Sub(pIn[i], pRem[i], r.mod)
		}
	} else {
		copy(pOut, pIn[:rank])
	}
}
