package rlwe

// checkPolyShape panics if p is not consistent with given parameters.
func checkPolyShape(params Parameters, p *Element) {
	if p.Rank() != params.Rank() || p.AuxModLen() > len(params.auxMod) || p.BaseModLen() > len(params.baseMod) {
		panic("inconsistent inputs")
	}
}
