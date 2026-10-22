package rlwe

// // checkScalarShape panics if c is not consistent with given parameters.
// func checkScalarShape(params Parameters, c *Scalar) {
// 	if c.AuxModLen() > len(params.auxMod) || c.BaseModLen() > len(params.baseMod) {
// 		panic("inconsistent inputs")
// 	}
// }

// checkPolyShape panics if p is not consistent with given parameters.
func checkPolyShape(params Parameters, p *Poly) {
	if p.Rank() != params.Rank() || p.AuxModLen() > len(params.auxMod) || p.BaseModLen() > len(params.baseMod) {
		panic("inconsistent inputs")
	}
}
