package crt

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
