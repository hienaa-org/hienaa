package dft

import (
	"math"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// checkLength checks if all vectors have the same length,
// and panics if not.
func checkLength(xs ...int) {
	if len(xs) == 0 {
		return
	}

	for i := 1; i < len(xs); i++ {
		if xs[i] != xs[0] {
			panic("inconsistent inputs")
		}
	}
}

// CyclotomicPolynomial computes the cyclotomic polynomial of the given cyclotomic index.
func CyclotomicPolynomial(cycloIdx int) []int64 {
	switch {
	case cycloIdx <= 0:
		panic("cycloIdx must be positive")

	case cycloIdx == 1:
		return []int64{-1, 1}

	case num.IsPowerOfTwo(cycloIdx):
		cycloPoly := make([]int64, cycloIdx>>1+1)
		cycloPoly[0] = 1
		cycloPoly[cycloIdx>>1] = 1
		return cycloPoly

	case num.IsPrime(cycloIdx):
		cycloPoly := make([]int64, cycloIdx)
		for i := range cycloPoly {
			cycloPoly[i] = 1
		}
		return cycloPoly
	}

	primes, exps := num.Factor(cycloIdx)
	phi := num.Totient(cycloIdx)

	var twoExp int
	if primes[0] == 2 {
		twoExp = exps[0]
		primes = primes[1:]
	}

	cycloIdxSqFree := 1
	phiSqFree := 1
	for i := range primes {
		cycloIdxSqFree *= primes[i]
		phiSqFree *= primes[i] - 1
	}

	mobius := make([]int64, cycloIdxSqFree+1)
	mobius[1] = 1
	for i := 1; i <= cycloIdxSqFree; i++ {
		j := 2 * i
		for j <= cycloIdxSqFree {
			mobius[j] -= mobius[i]
			j += i
		}
	}

	degSqFree := phiSqFree / 2
	cycloPolySqFree := make([]int64, phiSqFree+1)
	cycloPolySqFree[0] = 1
	for d := 1; d < cycloIdxSqFree; d++ {
		if cycloIdxSqFree%d != 0 {
			continue
		}
		if mobius[cycloIdxSqFree/d] == 1 {
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

	gap := (cycloIdx / cycloIdxSqFree) >> twoExp
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
		panic("length dividend must be longer than or equal to divisor")
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
			vec.MulSubTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], mod)
		}
	}

	return quo
}

// ExpansionFactor estimates the expansion factor of the ring, defined as
// max_{f, g} ||fg||_∞ / ||f||_∞ ||g||_∞.
//
// The computation may be slow, and the result may not be exact due to numerical errors and estimations.
// If the factor is too large, it returns +Inf.
func ExpansionFactor(params RingParameters) float64 {
	return ExpansionFactorSparse(params, params.rank)
}

// ExpansionFactorSparse estimates the expansion factor of the ring, defined as
// max_{f, g} ||fg||_∞ / ||f||_∞ ||g||_∞, when the hamming weight of f is hw.
//
// The computation may be slow, and the result may not be exact due to numerical errors and estimations.
// If the factor is too large, it returns +Inf.
func ExpansionFactorSparse(params RingParameters, hw int) float64 {
	if hw < 0 || hw > params.rank {
		panic("hamming weight must be in [0, rank]")
	}

	if hw == 0 {
		return 0
	}

	switch params.ringType {
	case TypeCyclotomic:
		return cyclotomicExpansionFactor(params, hw)
	case TypeCyclic:
		return cyclicExpansionFactor(hw)
	case TypeAutFixed:
		return autFixedExpansionFactor(params, hw)
	case TypeOther:
		return otherExpansionFactor(params, hw)
	}
	return math.Inf(1)
}

func cyclotomicExpansionFactor(params RingParameters, hw int) float64 {
	if params.cycloIdx == 1 {
		return 1
	}

	if num.IsPowerOfTwo(params.cycloIdx) {
		return float64(hw)
	}

	if num.IsPrime(params.cycloIdx) {
		return min(2*float64(hw), 2*float64(params.rank)-1)
	}

	return otherExpansionFactor(params, hw)
}

func cyclicExpansionFactor(hw int) float64 {
	return float64(hw)
}

func autFixedExpansionFactor(params RingParameters, hw int) float64 {
	if num.IsPowerOfTwo(params.cycloIdx) {
		return min(2*float64(hw), 2*float64(params.rank)-1)
	}

	fold := (params.cycloIdx - 1) / params.rank
	return min(2*float64(fold)*float64(hw), 2*float64(params.cycloIdx-fold)-1)
}

func otherExpansionFactor(params RingParameters, hw int) float64 {
	modPoly := vec.Cast[float64](params.modPoly)
	rem := make([]float64, params.rank)
	sums := make([]float64, params.rank)

	for i := 0; i < params.rank; i++ {
		sums[i] = float64(min(i+1, hw))
	}
	rem[params.rank-1] = 1

	for d := params.rank; d < 2*params.rank-1; d++ {
		lc := rem[params.rank-1]
		if math.IsInf(lc, 0) || math.IsNaN(lc) {
			return math.Inf(1)
		}
		weight := float64(min((2*params.rank-1)-d, hw))
		for i := params.rank - 1; i > 0; i-- {
			rem[i] = rem[i-1] - lc*modPoly[i]
			sums[i] += math.Abs(rem[i]) * weight
		}
		rem[0] = -lc * modPoly[0]
		sums[0] += math.Abs(rem[0]) * weight
	}

	var factor float64
	for _, sum := range sums {
		factor = max(factor, sum)
	}

	if math.IsNaN(factor) {
		return math.Inf(1)
	}
	return factor
}
