package rlwecheck

import (
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// CiphertextShape is an interface for checking the shape of encryptions.
type CiphertextShape interface {
	Rank() int
	BaseModLen() int
	AuxModLen() int
	Form() dft.Form
}

// CheckCiphertextShape panics if ct is not consistent.
func CheckCiphertextShape(params rlwe.Parameters, ct CiphertextShape) {
	if (ct.Rank() != params.Rank()) || !(0 <= ct.BaseModLen() && ct.BaseModLen() <= params.BaseModLen()) {
		panic("inconsistent inputs")
	}
	switch ct := ct.(type) {
	case *rlwe.GadgetEncryption:
		if len(ct.Value) != params.GadgetLen() || ct.AuxModLen() != params.AuxModLen() {
			panic("inconsistent inputs")
		}
	case *rlwe.RGSWEncryption:
		if len(ct.BodyStack.Value) != params.GadgetLen() || len(ct.MaskStack.Value) != params.GadgetLen() || ct.AuxModLen() != params.AuxModLen() {
			panic("inconsistent inputs")
		}
	default:
		if ct.AuxModLen() != 0 && ct.AuxModLen() != params.AuxModLen() {
			panic("inconsistent inputs")
		}
	}

	ct.Form()
}

// ElementShape is an interface for checking the shape of elements.
type ElementShape interface {
	Rank() int
	BaseModLen() int
	AuxModLen() int
	Form() dft.Form
	Type() crt.ElementType
}

// CheckElementShape panics if e is not consistent.
func CheckElementShape(params rlwe.Parameters, e ElementShape) {
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

// CheckBinaryOperable panics if eOut, e0, e1 is not operable.
func CheckBinaryOperable(params rlwe.Parameters, eOut, e0, e1 *rlwe.Element) {
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

	CheckElementShape(params, eOut)
	CheckElementShape(params, e0)
	CheckElementShape(params, e1)
}

// CheckUnaryOperable panics if eOut, e is not operable.
func CheckUnaryOperable(params rlwe.Parameters, eOut, e *rlwe.Element) {
	if eOut.Type() != e.Type() {
		panic("inconsistent types")
	}

	if eOut.BaseModLen() != e.BaseModLen() || eOut.AuxModLen() != e.AuxModLen() {
		panic("inconsistent inputs")
	}

	CheckElementShape(params, eOut)
	CheckElementShape(params, e)
}
