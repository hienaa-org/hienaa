package crt

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type addSubType uint64

const (
	typeTrivialAddSub addSubType = iota
	typePrimeAutFixedAddSub
)

func addSubTypeOf(params dft.RingParameters) addSubType {
	if params.RingType() == dft.TypeAutFixed && num.IsPrime(params.CycloIndex()) {
		return typePrimeAutFixedAddSub
	}
	return typeTrivialAddSub
}

// trivialAddSubOperator is an AddSubOperator for [typeTrivialAddSub].
type trivialAddSubOperator struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newTrivialAddSubOperator creates a new [baseAddSubOperator].
func newTrivialAddSubOperator(params dft.RingParameters, mod []*num.Modulus) *trivialAddSubOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &trivialAddSubOperator{
		rank:          params.Rank(),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

func (op *trivialAddSubOperator) addTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Add(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		for i := range op.mod {
			vec.AddTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
		}
		eOutPoly.IsNTT = e0Poly.IsNTT

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		for i := range op.mod {
			if p.IsNTT && op.isNTTFriendly[i] {
				vec.AddTo(eOutPoly.Coeffs[i], p.Coeffs[i], num.MForm(c.Value[i], op.mod[i]), op.mod[i])
			} else {
				copy(eOutPoly.Coeffs[i], p.Coeffs[i])
				eOutPoly.Coeffs[i][0] = num.Add(eOutPoly.Coeffs[i][0], c.Value[i], op.mod[i])
			}
		}
		eOutPoly.IsNTT = p.IsNTT
	}
}

func (op *trivialAddSubOperator) subTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Sub(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		for i := range op.mod {
			vec.SubTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
		}
		eOutPoly.IsNTT = e0Poly.IsNTT

	case e0IsScalar && !e1IsScalar:
		for i := range op.mod {
			if e1Poly.IsNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOutPoly.Coeffs[i], num.MForm(e0Scalar.Value[i], op.mod[i]), e1Poly.Coeffs[i], op.mod[i])
			} else {
				vec.NegTo(eOutPoly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
				eOutPoly.Coeffs[i][0] = num.Add(eOutPoly.Coeffs[i][0], e0Scalar.Value[i], op.mod[i])
			}
		}
		eOutPoly.IsNTT = e1Poly.IsNTT

	case !e0IsScalar && e1IsScalar:
		for i := range op.mod {
			if e0Poly.IsNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], num.MForm(e1Scalar.Value[i], op.mod[i]), op.mod[i])
			} else {
				copy(eOutPoly.Coeffs[i], e0Poly.Coeffs[i])
				eOutPoly.Coeffs[i][0] = num.Sub(eOutPoly.Coeffs[i][0], e1Scalar.Value[i], op.mod[i])
			}
		}
		eOutPoly.IsNTT = e0Poly.IsNTT
	}
}

func (op *trivialAddSubOperator) withModIdx(idx ...int) *trivialAddSubOperator {
	return &trivialAddSubOperator{
		rank:          op.rank,
		mod:           vec.Gather(op.mod, idx...),
		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),
	}
}

func (op *trivialAddSubOperator) slice(lo, hi int) *trivialAddSubOperator {
	return &trivialAddSubOperator{
		rank:          op.rank,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],
	}
}

func (op *trivialAddSubOperator) append(op0 *trivialAddSubOperator) *trivialAddSubOperator {
	return &trivialAddSubOperator{
		rank:          op0.rank,
		mod:           slices.Concat(op.mod, op0.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),
	}
}

func (op *trivialAddSubOperator) appendTmpModulus(mod *num.Modulus) *trivialAddSubOperator {
	return &trivialAddSubOperator{
		rank:          op.rank,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),
	}
}

// primeAutFixedAddSubOperator is an AddSubOperator for [typePrimeAutFixedAddSub].
type primeAutFixedAddSubOperator struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newPrimeAutFixedAddSubOperator creates a new [primeAutFixedAddSubOperator].
func newPrimeAutFixedAddSubOperator(params dft.RingParameters, mod []*num.Modulus) *primeAutFixedAddSubOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &primeAutFixedAddSubOperator{
		rank:          params.Rank(),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

func (op *primeAutFixedAddSubOperator) addTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Add(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		for i := range op.mod {
			vec.AddTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
		}
		eOutPoly.IsNTT = e0Poly.IsNTT

	default:
		c, p := orderByType(e0IsScalar, e0Scalar, e0Poly, e1Scalar, e1Poly)
		for i := range op.mod {
			if p.IsNTT && op.isNTTFriendly[i] {
				vec.AddTo(eOutPoly.Coeffs[i], p.Coeffs[i], num.MForm(c.Value[i], op.mod[i]), op.mod[i])
			} else {
				vec.SubTo(eOutPoly.Coeffs[i], p.Coeffs[i], c.Value[i], op.mod[i])
			}
		}
		eOutPoly.IsNTT = p.IsNTT
	}
}

func (op *primeAutFixedAddSubOperator) subTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	e0Scalar, e0IsScalar := any(e0).(*Scalar)
	e0Poly, _ := any(e0).(*Poly)

	e1Scalar, e1IsScalar := any(e1).(*Scalar)
	e1Poly, _ := any(e1).(*Poly)

	switch {
	case e0IsScalar && e1IsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Sub(e0Scalar.Value[i], e1Scalar.Value[i], op.mod[i])
		}

	case !e0IsScalar && !e1IsScalar:
		for i := range op.mod {
			vec.SubTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Poly.Coeffs[i], op.mod[i])
		}
		eOutPoly.IsNTT = e0Poly.IsNTT

	case e0IsScalar && !e1IsScalar:
		for i := range op.mod {
			if e1Poly.IsNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOutPoly.Coeffs[i], num.MForm(e0Scalar.Value[i], op.mod[i]), e1Poly.Coeffs[i], op.mod[i])
			} else {
				vec.SubTo(eOutPoly.Coeffs[i], num.Neg(e0Scalar.Value[i], op.mod[i]), e1Poly.Coeffs[i], op.mod[i])
			}
		}
		eOutPoly.IsNTT = e1Poly.IsNTT

	case !e0IsScalar && e1IsScalar:
		for i := range op.mod {
			if e0Poly.IsNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], num.MForm(e1Scalar.Value[i], op.mod[i]), op.mod[i])
			} else {
				vec.AddTo(eOutPoly.Coeffs[i], e0Poly.Coeffs[i], e1Scalar.Value[i], op.mod[i])
			}
		}
		eOutPoly.IsNTT = e0Poly.IsNTT
	}
}

func (op *primeAutFixedAddSubOperator) withModIdx(idx ...int) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           vec.Gather(op.mod, idx...),
		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),
	}
}

func (op *primeAutFixedAddSubOperator) slice(lo, hi int) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],
	}
}

func (op *primeAutFixedAddSubOperator) append(op0 *primeAutFixedAddSubOperator) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op0.rank,
		mod:           slices.Concat(op.mod, op0.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),
	}
}

func (op *primeAutFixedAddSubOperator) appendTmpModulus(mod *num.Modulus) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),
	}
}
