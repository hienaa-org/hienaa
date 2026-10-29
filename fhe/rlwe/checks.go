package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// ciphertextShape is an interface for checking the shape of encryptions.
type ciphertextShape interface {
	Rank() int
	BaseModLen() int
	AuxModLen() int
	Form() dft.Form
}

// checkCiphertextShape panics if ct is not consistent.
func checkCiphertextShape(params Parameters, ct ciphertextShape) {
	if (ct.Rank() != params.Rank()) || !(0 <= ct.BaseModLen() && ct.BaseModLen() <= params.BaseModLen()) {
		panic("inconsistent inputs")
	}
	switch ct := ct.(type) {
	case *GadgetEncryption, *RGSWEncryption:
		if ct.AuxModLen() != params.AuxModLen() {
			panic("inconsistent inputs")
		}
	default:
		if ct.AuxModLen() != 0 && ct.AuxModLen() != params.AuxModLen() {
			panic("inconsistent inputs")
		}
	}

	ct.Form()
}

// elementShape is an interface for checking the shape of elements.
type elementShape interface {
	Rank() int
	BaseModLen() int
	AuxModLen() int
	Form() dft.Form
	Type() crt.ElementType
}

// checkElementShape panics if e is not consistent.
func checkElementShape(params Parameters, e elementShape) {
	var expectedRank int
	if e.Type() == crt.TypeScalar {
		expectedRank = 1
	} else {
		expectedRank = params.Rank()
	}

	if e.Rank() != expectedRank || !(e.AuxModLen() == 0 || e.AuxModLen() == params.AuxModLen()) ||
		!(0 <= e.BaseModLen() && e.BaseModLen() <= params.BaseModLen()) {
		panic("inconsistent inputs")
	}
}

// newBinaryOpOut returns eOut given e0, e1 in binary operations.
func newBinaryOpOut(e0, e1 *Element) *Element {
	if e0.Type() == crt.TypeScalar && e1.Type() == crt.TypeScalar {
		return NewScalarCustom(e0.BaseModLen(), e0.AuxModLen())
	}
	if e0.Type() == crt.TypeScalar {
		return NewPolyCustom(e1.Rank(), e1.BaseModLen(), e1.AuxModLen(), dft.FormCoeff)
	}

	return NewPolyCustom(e0.Rank(), e0.BaseModLen(), e0.AuxModLen(), dft.FormCoeff)
}

// checkBinaryOperable panics if eOut, e0, e1 is not operable.
func checkBinaryOperable(params Parameters, eOut, e0, e1 *Element) {
	if e0.Type() == crt.TypeScalar && e1.Type() == crt.TypeScalar {
		if eOut.Type() != crt.TypeScalar {
			panic("inconsistent types")
		}
	} else {
		if eOut.Type() != crt.TypePoly {
			panic("inconsistent types")
		}
		if e0.Type() == crt.TypePoly && e1.Type() == crt.TypePoly && e0.Form() != e1.Form() {
			panic("inconsistent forms")
		}
	}

	if eOut.BaseModLen() != e0.BaseModLen() || eOut.AuxModLen() != e0.AuxModLen() ||
		eOut.BaseModLen() != e1.BaseModLen() || eOut.AuxModLen() != e1.AuxModLen() {
		panic("inconsistent inputs")
	}

	checkElementShape(params, eOut)
	checkElementShape(params, e0)
	checkElementShape(params, e1)
}

// newUnaryOpOut returns eOut given e in unary operations.
func newUnaryOpOut(e *Element) *Element {
	if e.Type() == crt.TypeScalar {
		return NewScalarCustom(e.BaseModLen(), e.AuxModLen())
	}
	return NewPolyCustom(e.Rank(), e.BaseModLen(), e.AuxModLen(), dft.FormCoeff)
}

// checkUnaryOperable panics if eOut, e is not operable.
func checkUnaryOperable(params Parameters, eOut, e *Element) {
	if eOut.Type() != e.Type() {
		panic("inconsistent types")
	}

	if eOut.BaseModLen() != e.BaseModLen() || eOut.AuxModLen() != e.AuxModLen() {
		panic("inconsistent inputs")
	}

	checkElementShape(params, eOut)
	checkElementShape(params, e)
}
