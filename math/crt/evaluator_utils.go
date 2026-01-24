package crt

import "github.com/hienaa-org/hienaa/math/num"

// isCoprime checks if all modulus are coprime.
func isCoprime(mod []*num.Modulus) bool {
	gcd := uint64(1)
	for i := range mod {
		gcd = num.GCD(gcd, mod[i].Value())
	}
	return gcd == 1
}

// mustConsistent panics if p is not consistent with given ring parameters and modulus.
func mustConsistent(rank, modLen int, p *Poly) {
	if len(p.Coeffs) != modLen {
		panic("input(s) not consistent")
	}

	for i := 0; i < modLen; i++ {
		if len(p.Coeffs[i]) != rank {
			panic("input(s) not consistent")
		}
	}
}

// mustTernaryToOperable panics if pOut, p0, p1 is not operable.
// It checks:
//
//   - pOut, p0, p1 has same shape.
//   - p0, p1 has same form.
func mustTernaryToOperable(rank, modLen int, pOut, p0, p1 *Poly) {
	mustConsistent(rank, modLen, pOut)
	mustConsistent(rank, modLen, p0)
	mustConsistent(rank, modLen, p1)
	if p0.IsNTT != p1.IsNTT {
		panic("input(s) not consistent")
	}
}

// mustBinaryToOperable panics if pOut, p is not operable.
// It checks:
//
//   - pOut, p has same shape.
func mustBinaryToOperable(rank, modLen int, pOut, p *Poly) {
	mustConsistent(rank, modLen, pOut)
	mustConsistent(rank, modLen, p)
}
