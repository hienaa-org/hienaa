package crt

import "github.com/hienaa-org/hienaa/math/dft"

// checkElementShape panics if e is not consistent with given ring parameters and modulus.
func checkElementShape(rank, modLen int, e *Element) {
	if e.Type == TypeScalar {
		if e.Rank() != 1 || e.ModLen() != modLen {
			panic("inconsistent inputs")
		}
	} else {
		checkPolyShape(rank, modLen, e)
	}
}

// checkPolyShape panics if p does not have the given rank and modulus length.
func checkPolyShape(rank, modLen int, p *Element) {
	if p.Rank() != rank || p.ModLen() != modLen {
		panic("inconsistent inputs")
	}
}

// newBinaryOpOut returns eOut given e0, e1 in binary operations.
func newBinaryOpOut(e0, e1 *Element) *Element {
	if e0.Type == TypeScalar && e1.Type == TypeScalar {
		return NewScalar(e0.ModLen())
	}
	if e0.Type == TypeScalar {
		return NewPoly(e1.Rank(), e1.ModLen(), dft.FormCoeff)
	}

	return NewPoly(e0.Rank(), e0.ModLen(), dft.FormCoeff)
}

// checkBinaryOperable panics if eOut, e0, e1 is not operable.
func checkBinaryOperable(rank, modLen int, eOut, e0, e1 *Element) {
	if e0.Type == TypeScalar && e1.Type == TypeScalar {
		if eOut.Type != TypeScalar {
			panic("inconsistent types")
		}
	} else {
		if eOut.Type != TypePoly {
			panic("inconsistent types")
		}
		if e0.Type == TypePoly && e1.Type == TypePoly && e0.Form != e1.Form {
			panic("inconsistent forms")
		}
	}

	checkElementShape(rank, modLen, eOut)
	checkElementShape(rank, modLen, e0)
	checkElementShape(rank, modLen, e1)
}

// newUnaryOpOut returns eOut given e in unary operations.
func newUnaryOpOut(e *Element) *Element {
	if e.Type == TypeScalar {
		return NewScalar(e.ModLen())
	}
	return NewPoly(e.Rank(), e.ModLen(), dft.FormCoeff)
}

// checkUnaryOperable panics if eOut, e is not operable.
func checkUnaryOperable(rank, modLen int, eOut, e *Element) {
	if eOut.Type != e.Type {
		panic("inconsistent types")
	}

	checkElementShape(rank, modLen, eOut)
	checkElementShape(rank, modLen, e)
}

// orderByType returns e0, e1 as the order of [TypeScalar] and [TypePoly].
// Assumes that one of e0, e1 is [TypeScalar] and the other is [TypePoly].
func orderByType(e0, e1 *Element) (eScalar *Element, ePoly *Element) {
	if e0.Type == TypeScalar {
		return e0, e1
	}
	return e1, e0
}
