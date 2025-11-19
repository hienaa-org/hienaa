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

// isConsistent checks if p is consistent with given ring parameters and modulus.
func isConsistent(rank, modLen int, p *Poly) bool {
	if len(p.Coeffs) != modLen {
		return false
	}

	for i := 0; i < modLen; i++ {
		if len(p.Coeffs[i]) != rank {
			return false
		}
	}
	return true
}

// isTernaryToOperable checks if pOut, p0, p1 is operatable.
// It checks:
//
//   - pOut, p0, p1 has same shape.
//   - p0, p1 has same form.
func isTernaryToOperable(rank, modLen int, pOut, p0, p1 *Poly) bool {
	return isConsistent(rank, modLen, pOut) && isConsistent(rank, modLen, p0) && isConsistent(rank, modLen, p1) && p0.IsNTT == p1.IsNTT
}

// isBinaryToOperable checks if pOut, p is operatable.
// It checks:
//
//   - pOut, p has same shape.
func isBinaryToOperable(rank, modLen int, pOut, p *Poly) bool {
	return isConsistent(rank, modLen, pOut) && isConsistent(rank, modLen, p)
}
