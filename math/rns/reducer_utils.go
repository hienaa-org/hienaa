package rns

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

func reduceInt(x int, q *mod.Modulus) uint64 {
	if x < 0 {
		return mod.Neg(uint64(-x), q)
	} else {
		return mod.Reduce(uint64(x), q)
	}
}

// quotientPolynomial computes the quotient of two polynomials.
func quotientPolynomial(dividend, divisor []int) []int {
	if len(dividend) < len(divisor) {
		panic("dividend is shorter than divisor")
	}

	quotient := make([]int, len(dividend)-len(divisor)+1)
	pBuff := make([]int, len(dividend))
	copy(pBuff, dividend)

	for i := 0; i <= len(dividend)-len(divisor); i++ {
		if pBuff[len(pBuff)-i-1] != 0 {
			quotient[len(quotient)-i-1] = pBuff[len(pBuff)-i-1] / divisor[len(divisor)-1]

			for j := 0; j < len(divisor); j++ {
				pBuff[len(pBuff)-i-j-1] -= divisor[len(divisor)-j-1] * quotient[len(quotient)-i-1]
			}
		}
	}

	return quotient
}

// quotientPolynomialMod computes the quotient of two polynomials modulo a modulus.
func quotientPolynomialMod(dividend, divisor []uint64, modulus *mod.Modulus) []uint64 {
	if len(dividend) < len(divisor) {
		panic("dividend is shorter than divisor")
	}
	if num.GCD(modulus.Value(), divisor[len(divisor)-1]) != 1 {
		panic("divisor is not coprime with modulus")
	}

	quotient := make([]uint64, len(dividend)-len(divisor)+1)
	pBuff := make([]uint64, len(dividend))
	copy(pBuff, dividend)

	for i := 0; i <= len(dividend)-len(divisor); i++ {
		if pBuff[len(pBuff)-i-1] != 0 {
			quotient[len(quotient)-i-1] = mod.Mul(pBuff[len(pBuff)-i-1], mod.Inv(divisor[len(divisor)-1], modulus), modulus)

			for j := 0; j < len(divisor); j++ {
				pBuff[len(pBuff)-i-j-1] = mod.Sub(pBuff[len(pBuff)-i-j-1], mod.Mul(divisor[len(divisor)-j-1], quotient[len(quotient)-i-1], modulus), modulus)
			}
		}
	}

	return quotient
}

// computeCyclotomicPolynomial computes the cyclotomic polynomial of the given degree.
func computeCyclotomicPolynomial(degree uint64) []int {
	factors := num.Factor(degree)
	primes := make([]int, 0, len(factors))
	for key := range factors {
		primes = append(primes, int(key))
	}
	slices.Sort(primes)

	var isEven bool
	if primes[0] == 2 {
		isEven = true
		primes = primes[1:]
	} else {
		isEven = false
	}

	pOut := make([]int, degree+1)
	pBuff0 := make([]int, degree+1)
	pBuff1 := make([]int, degree+1)
	pOut[0], pOut[1] = -1, 1

	skip := int(degree)

	currDeg := 1
	prevDeg := 1
	for _, prime := range primes {
		copy(pBuff0, pOut)
		clear(pBuff1)
		clear(pOut)

		for i := 0; i <= prevDeg; i++ {
			pBuff1[i*prime] = pBuff0[i]
		}

		currDeg = prevDeg*prime - prevDeg

		for i := 0; i <= (prime-1)*prevDeg; i++ {
			if pBuff1[prevDeg*prime-i] != 0 {
				pOut[currDeg-i] = pBuff1[prevDeg*prime-i] / pBuff0[prevDeg]

				for j := 0; j <= prevDeg; j++ {
					pBuff1[prevDeg*prime-i-j] -= pBuff0[prevDeg-j] * pOut[currDeg-i]
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
		copy(pBuff0, pOut)
		clear(pOut)
		for i := 0; i <= prevDeg; i++ {
			pOut[i*skip] = pBuff0[i]
		}
	}

	return pOut[:num.Totient(uint64(degree))+1]
}
