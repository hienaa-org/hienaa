package dft

import (
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
			panic("inconsistent input(s)")
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
		exps = exps[1:]
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
			vec.MulSubTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], mod)
		}
	}

	return quo
}
