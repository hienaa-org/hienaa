package pack

import (
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// quoRem computes the quotient and remainder of two polynomials modulo a modulus.
func quoRem(p0, p1 []uint64, mod *num.Modulus) ([]uint64, []uint64) {
	switch {
	case len(p0) < len(p1):
		panic("quotient: dividend is shorter than divisor")
	case num.GCD(mod.Value(), p1[len(p1)-1]) != 1:
		panic("quotient: divisor is not coprime with modulus")
	}

	quo := make([]uint64, len(p0)-len(p1)+1)
	rem := make([]uint64, len(p0))
	copy(rem, p0)

	lcInv := num.Inv(p1[len(p1)-1], mod)
	for i := 0; i <= len(p0)-len(p1); i++ {
		if rem[len(rem)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(rem[len(rem)-i-1], lcInv, mod)
			vec.ScalarMulSubTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], mod)
		}
	}

	return quo, rem[:len(p1)-1]
}
