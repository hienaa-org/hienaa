package crt

// checkPolyShape panics if p is not consistent with given ring parameters and modulus.
func checkPolyShape(rank, modLen int, p *Poly) {
	if len(p.Coeffs) != modLen {
		panic("input(s) shape not consistent")
	}

	for i := 0; i < modLen; i++ {
		if len(p.Coeffs[i]) != rank {
			panic("input(s) shape not consistent")
		}
	}
}
