package crt

// checkScalarShape panics if c is not consistent with given modulus.
func checkScalarShape(modLen int, c *Scalar) {
	if c.ModLen() != modLen {
		panic("input(s) not consistent")
	}
}

// checkPolyShape panics if p is not consistent with given ring parameters and modulus.
func checkPolyShape(rank, modLen int, p *Poly) {
	if p.Rank() != rank || p.ModLen() != modLen {
		panic("input(s) not consistent")
	}
}

// newBinaryOpOut returns eOut given e0, e1 in binary operations.
func newBinaryOpOut[TOut, T0, T1 *Scalar | *Poly](e0 T0, e1 T1) TOut {
	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	_, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	if e0IsScalar && e1IsScalar {
		return any(NewScalar(e0Scalar.ModLen())).(TOut)
	}

	if !e0IsScalar {
		return any(NewPoly(e0Poly.Rank(), e0Poly.ModLen())).(TOut)
	}
	return any(NewPoly(e1Poly.Rank(), e1Poly.ModLen())).(TOut)
}

// isBinaryOperable panics if eOut, e0, e1 is not operable.
func isBinaryOperable[TOut, T0, T1 *Scalar | *Poly](rank, modLen int, eOut TOut, e0 T0, e1 T1) {
	eOutScalar, eOutIsScalar := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		if !eOutIsScalar {
			panic("inconsistent input(s)")
		}
		checkScalarShape(modLen, e0Scalar)
		checkScalarShape(modLen, e1Scalar)
		checkScalarShape(modLen, eOutScalar)

	case e0IsScalar && !e1IsScalar:
		if eOutIsScalar {
			panic("inconsistent input(s)")
		}
		checkScalarShape(modLen, e0Scalar)
		checkPolyShape(rank, modLen, e1Poly)
		checkPolyShape(rank, modLen, eOutPoly)

	case !e0IsScalar && e1IsScalar:
		if eOutIsScalar {
			panic("inconsistent input(s)")
		}
		checkPolyShape(rank, modLen, e0Poly)
		checkScalarShape(modLen, e1Scalar)
		checkPolyShape(rank, modLen, eOutPoly)

	case !e0IsScalar && !e1IsScalar:
		if eOutIsScalar {
			panic("inconsistent input(s)")
		}
		if e0Poly.IsNTT != e1Poly.IsNTT {
			panic("inconsistent NTT flags")
		}
		checkPolyShape(rank, modLen, e0Poly)
		checkPolyShape(rank, modLen, e1Poly)
		checkPolyShape(rank, modLen, eOutPoly)
	}
}

// newUnaryOpOut returns eOut given e in unary operations.
func newUnaryOpOut[TOut, T *Scalar | *Poly](e T) TOut {
	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	if eIsScalar {
		return any(NewScalar(eScalar.ModLen())).(TOut)
	}
	return any(NewPoly(ePoly.Rank(), ePoly.ModLen())).(TOut)
}

// isUnaryOperable panics if eOut, e is not operable.
func isUnaryOperable[TOut, T *Scalar | *Poly](rank, modLen int, eOut TOut, e T) {
	eOutScalar, eOutIsScalar := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	if eOutIsScalar != eIsScalar {
		panic("inconsistent input(s)")
	}

	if eOutIsScalar && eIsScalar {
		checkScalarShape(modLen, eOutScalar)
		checkScalarShape(modLen, eScalar)
	}

	if !eOutIsScalar && !eIsScalar {
		checkPolyShape(rank, modLen, eOutPoly)
		checkPolyShape(rank, modLen, ePoly)
	}
}

// orderByType returns e0, e1 as the order of [*Scalar] and [*Poly].
// Assumes that one of e0, e1 is [*Scalar] and the other is [*Poly].
func orderByType(e0IsScalar bool, e0Scalar *Scalar, e0Poly *Poly, e1Scalar *Scalar, e1Poly *Poly) (*Scalar, *Poly) {
	if e0IsScalar {
		return e0Scalar, e1Poly
	}
	return e1Scalar, e0Poly
}
