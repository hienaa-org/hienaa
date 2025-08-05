package dft

import (
	"github.com/hienaa-org/hienaa/math/num"
)

func reduceInt(x int, q *num.Modulus) uint64 {
	if x < 0 {
		return num.Neg(uint64(-x), q)
	} else {
		return num.Reduce(uint64(x), q)
	}
}

// quotient computes the quotient of two polynomials modulo a modulus.
func quotient(dividend, divisor []uint64, mod *num.Modulus) []uint64 {
	switch {
	case len(dividend) < len(divisor):
		panic("quotient: dividend is shorter than divisor")
	case num.GCD(mod.Value(), divisor[len(divisor)-1]) != 1:
		panic("quotient: divisor is not coprime with modulus")
	}

	quo := make([]uint64, len(dividend)-len(divisor)+1)
	pBuff := make([]uint64, len(dividend))
	copy(pBuff, dividend)

	for i := 0; i <= len(dividend)-len(divisor); i++ {
		if pBuff[len(pBuff)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(pBuff[len(pBuff)-i-1], num.Inv(divisor[len(divisor)-1], mod), mod)

			for j := 0; j < len(divisor); j++ {
				pBuff[len(pBuff)-i-j-1] = num.Sub(pBuff[len(pBuff)-i-j-1], num.Mul(divisor[len(divisor)-j-1], quo[len(quo)-i-1], mod), mod)
			}
		}
	}

	return quo
}

// CyclotomicPolynomial computes the cyclotomic polynomial of the given cyclotomic order.
func CyclotomicPolynomial(cycloOrd int) []int {
	primes, _ := num.Factor(uint64(cycloOrd))

	var isEven bool
	if primes[0] == 2 {
		isEven = true
		primes = primes[1:]
	} else {
		isEven = false
	}

	pOut := make([]int, cycloOrd+1)
	pBuff0 := make([]int, cycloOrd+1)
	pBuff1 := make([]int, cycloOrd+1)
	pOut[0], pOut[1] = -1, 1

	skip := int(cycloOrd)

	currDeg, prevDeg := 1, 1
	for _, prime := range primes {
		copy(pBuff0, pOut)
		clear(pBuff1)
		clear(pOut)

		for i := 0; i <= prevDeg; i++ {
			pBuff1[i*int(prime)] = pBuff0[i]
		}

		currDeg = prevDeg*int(prime) - prevDeg

		for i := 0; i <= (int(prime)-1)*prevDeg; i++ {
			if pBuff1[prevDeg*int(prime)-i] != 0 {
				pOut[currDeg-i] = pBuff1[prevDeg*int(prime)-i] / pBuff0[prevDeg]

				for j := 0; j <= prevDeg; j++ {
					pBuff1[prevDeg*int(prime)-i-j] -= pBuff0[prevDeg-j] * pOut[currDeg-i]
				}
			}
		}

		prevDeg = currDeg
		skip /= int(prime)
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

	return pOut[:num.Totient(uint64(cycloOrd))+1]
}
