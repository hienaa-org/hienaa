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
		vec.ReduceTo(cycloPoly[:params.rank+1], params.modPoly, mod)

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

		pool: pool.NewPool(func() *[]uint64 {
			p := make([]uint64, max(params.cycloOrd, diffDegNext, degNext))
			return &p
		}),
	}
}

func (r *cyclotomicReducer) reduceTo(pOut, p []uint64) {
	cycloOrd, rank := r.params.cycloOrd, r.params.rank

	pInPtr := r.pool.Get()
	pIn := (*pInPtr)[:cycloOrd]
	defer r.pool.Put(pInPtr)

	copy(pIn, p)
	if cycloOrd%2 == 1 {
		skip := cycloOrd / r.leastFac

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				pIn[i*skip+j] = num.Sub(pIn[i*skip+j], pIn[cycloOrd-skip+j], r.mod)
			}
			pIn[cycloOrd-skip+j] = 0
		}
	} else {
		skip := (cycloOrd / 2) / r.leastFac

		for i := 0; i < cycloOrd/2; i++ {
			pIn[i] = num.Sub(pIn[i], pIn[cycloOrd/2+i], r.mod)
			pIn[cycloOrd/2+i] = 0
		}

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				if i%2 == 0 {
					pIn[i*skip+j] = num.Sub(pIn[i*skip+j], pIn[cycloOrd/2-skip+j], r.mod)
				} else {
					pIn[i*skip+j] = num.Add(pIn[i*skip+j], pIn[cycloOrd/2-skip+j], r.mod)
				}
			}
			pIn[cycloOrd/2-skip+j] = 0
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
		r.diffDegNextNTT.ForwardTo(pQuo, pQuo)
		vec.MMulLazyTo(pQuo, pQuo, r.divPoly, r.mod)
		r.diffDegNextNTT.InverseTo(pQuo, pQuo)

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
		r.degNextNTT.ForwardTo(pRem, pRem)
		vec.MMulLazyTo(pRem, pRem, r.cycloPoly, r.mod)
		r.degNextNTT.InverseTo(pRem, pRem)

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

// CyclotomicPolynomial computes the cyclotomic polynomial of the given cyclotomic order.
func CyclotomicPolynomial(cycloOrd int) []int64 {
	switch {
	case cycloOrd <= 0:
		panic("cycloOrd must be positive")
	case cycloOrd == 1:
		return []int64{-1, 1}
	case num.IsPowerOfTwo(cycloOrd):
		cycloPoly := make([]int64, cycloOrd>>1+1)
		cycloPoly[0] = 1
		cycloPoly[cycloOrd>>1] = 1
		return cycloPoly
	case num.IsPrime(cycloOrd):
		cycloPoly := make([]int64, cycloOrd)
		for i := range cycloPoly {
			cycloPoly[i] = 1
		}
		return cycloPoly
	}

	primes, exps := num.Factor(cycloOrd)
	phi := num.TotientWithFactors(cycloOrd, primes, exps)

	var twoExp int
	if primes[0] == 2 {
		twoExp = exps[0]
		primes = primes[1:]
		exps = exps[1:]
	}

	cycloOrdSqFree := 1
	phiSqFree := 1
	for i := range primes {
		cycloOrdSqFree *= primes[i]
		phiSqFree *= primes[i] - 1
	}

	mobius := make([]int, cycloOrdSqFree+1)
	mobius[1] = 1
	for i := 1; i <= cycloOrdSqFree; i++ {
		j := 2 * i
		for j <= cycloOrdSqFree {
			mobius[j] -= mobius[i]
			j += i
		}
	}

	degSqFree := phiSqFree / 2
	cycloPolySqFree := make([]int64, phiSqFree+1)
	cycloPolySqFree[0] = 1
	for d := 1; d < cycloOrdSqFree; d++ {
		if cycloOrdSqFree%d != 0 {
			continue
		}
		if mobius[cycloOrdSqFree/d] == 1 {
			for i := degSqFree; i >= d; i-- {
				cycloPolySqFree[i] -= cycloPolySqFree[i-d]
			}
		} else {
			for i := d; i <= degSqFree; i++ {
				cycloPolySqFree[i] += cycloPolySqFree[i-d]
			}
		}
	}

	for i := degSqFree + 1; i <= phiSqFree; i++ {
		cycloPolySqFree[i] = cycloPolySqFree[phiSqFree-i]
	}

	gap := (cycloOrd / cycloOrdSqFree) >> twoExp
	if twoExp >= 1 {
		gap <<= twoExp - 1
	}
	cycloPoly := make([]int64, phi+1)
	for i := 0; i <= phiSqFree; i++ {
		cycloPoly[i*gap] = cycloPolySqFree[i]
	}

	if twoExp >= 1 {
		for i := 1 << (twoExp - 1); i <= phi; i += 1 << twoExp {
			cycloPoly[i] = -cycloPoly[i]
		}
	}

	return cycloPoly
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
