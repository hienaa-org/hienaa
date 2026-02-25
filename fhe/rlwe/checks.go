package rlwe

// checkShape panics if e is not consistent with given parameters.
func checkShape(baseLen, auxLen int, e *Element) {
	if e.BaseModLen() != baseLen || e.AuxModLen() != auxLen {
		panic("input(s) shape not consistent")
	}
}

// checkOperable panics if eOut, e0, e1 is not operable.
func checkOperable(eOut *Element, e ...*Element) {
	baseLen := eOut.BaseModLen()
	auxLen := eOut.AuxModLen()

	for i := range e {
		checkShape(baseLen, auxLen, e[i])
	}
}
