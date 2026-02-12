package crt

import (
	"github.com/hienaa-org/hienaa/math/num"
)

// isCoprime checks if all modulus are coprime.
func isCoprime(mod []*num.Modulus) bool {
	gcd := uint64(1)
	for i := range mod {
		gcd = num.GCD(gcd, mod[i].Value())
	}
	return gcd == 1
}

// checkShape panics if p is not consistent with given ring parameters and modulus.
func checkShape(rank, modLen int, p *Element) {
	if len(p.Coeffs) != modLen {
		panic("input(s) shape not consistent")
	}

	for i := 0; i < modLen; i++ {
		if len(p.Coeffs[i]) != rank {
			panic("input(s) shape not consistent")
		}
	}
}

// checkTernaryOperable panics if eOut, e0, e1 is not operable.
func checkTernaryOperable(rank, modLen int, eOut, e0, e1 *Element) {
	if e0.Type() == TypeScalar && e1.Type() == TypeScalar {
		if eOut.Type() != TypeScalar {
			panic("output type not consistent")
		}
	} else {
		if eOut.Type() != TypePoly {
			panic("output type not consistent")
		}
		checkShape(rank, modLen, eOut)
		if e0.Type() == TypePoly {
			checkShape(rank, modLen, e0)
		}
		if e1.Type() == TypePoly {
			checkShape(rank, modLen, e1)
		}
		if e0.Type() == TypePoly && e1.Type() == TypePoly {
			if e0.IsNTT != e1.IsNTT {
				panic("input(s) NTT flag not consistent")
			}
		}
	}
}

// checkBinaryOperable panics if eOut, e is not operable.
func checkBinaryOperable(rank, modLen int, eOut, e *Element) {
	if e.Type() == TypeScalar {
		if eOut.Type() != TypeScalar {
			panic("output type not consistent")
		}
	} else {
		if eOut.Type() != TypePoly {
			panic("output type not consistent")
		}
		checkShape(rank, modLen, eOut)
		checkShape(rank, modLen, e)
	}
}

// isEqualType checks if given elements are of the same type.
func isEqualType(e0, e1 *Element, t ElementType) bool {
	return e0.Type() == t && e1.Type() == t
}

// orderByType returns e0, e1 as the order of scalar and poly.
// Assumes that one of e0, e1 is scalar and the other is poly.
func orderByType(e0, e1 *Element) (c *Element, p *Element) {
	if e0.Type() == TypeScalar {
		return e0, e1
	}
	return e1, e0
}
