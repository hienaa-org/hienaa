package crt

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
	return isConsistent(rank, modLen, pOut) && isConsistent(rank, modLen, p0) && isConsistent(rank, modLen, p1) && p0.isNTT == p1.isNTT
}

// isBinaryToOperable checks if pOut, p is operatable.
// It checks:
//
//   - pOut, p has same shape.
func isBinaryToOperable(rank, modLen int, pOut, p *Poly) bool {
	return isConsistent(rank, modLen, pOut) && isConsistent(rank, modLen, p)
}
