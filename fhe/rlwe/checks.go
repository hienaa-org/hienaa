package rlwe

// checkShape panics if e is not consistent with given parameters.
func checkShape(modLen int, hasAux bool, e *Element) {
	if e.hasAux != hasAux {
		panic("input(s) shape not consistent")
	}

	if len(e.Value.Coeffs) > modLen {
		panic("input(s) shape not consistent")
	}
}

// checkOperable panics if eOut, e0, e1 is not operable.
func checkOperable(modLen int, eOut *Element, e ...*Element) {
	checkShape(modLen, eOut.hasAux, eOut)
	for i := range e {
		checkShape(modLen, eOut.hasAux, e[i])
	}
}
