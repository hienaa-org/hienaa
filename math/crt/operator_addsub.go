package crt

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type addSubType byte

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

// newTrivialAddSubOperator creates a new [trivialAddSubOperator].
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

func (op *trivialAddSubOperator) addTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		for i := range op.mod {
			vec.AddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.Form = e0.Form

	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			if p.Form == dft.FormNTT && op.isNTTFriendly[i] {
				vec.AddTo(eOut.Coeffs[i], p.Coeffs[i], num.MForm(c.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				copy(eOut.Coeffs[i], p.Coeffs[i])
				eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], c.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.Form = p.Form
	}
}

func (op *trivialAddSubOperator) subTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		for i := range op.mod {
			vec.SubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.Form = e0.Form

	case e0.Type == TypeScalar && e1.Type == TypePoly:
		for i := range op.mod {
			if e1.Form == dft.FormNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOut.Coeffs[i], num.MForm(e0.Coeffs[i][0], op.mod[i]), e1.Coeffs[i], op.mod[i])
			} else {
				vec.NegTo(eOut.Coeffs[i], e1.Coeffs[i], op.mod[i])
				eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], e0.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.Form = e1.Form

	case e0.Type == TypePoly && e1.Type == TypeScalar:
		for i := range op.mod {
			if e0.Form == dft.FormNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOut.Coeffs[i], e0.Coeffs[i], num.MForm(e1.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				copy(eOut.Coeffs[i], e0.Coeffs[i])
				eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.Form = e0.Form
	}
}

func (op *trivialAddSubOperator) slice(lo, hi int) *trivialAddSubOperator {
	return &trivialAddSubOperator{
		rank:          op.rank,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],
	}
}

func (op *trivialAddSubOperator) append(opIn *trivialAddSubOperator) *trivialAddSubOperator {
	return &trivialAddSubOperator{
		rank:          opIn.rank,
		mod:           slices.Concat(op.mod, opIn.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, opIn.isNTTFriendly),
	}
}

func (op *trivialAddSubOperator) appendTmpMod(mod *num.Modulus) *trivialAddSubOperator {
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

func (op *primeAutFixedAddSubOperator) addTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		for i := range op.mod {
			vec.AddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.Form = e0.Form

	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			if p.Form == dft.FormNTT && op.isNTTFriendly[i] {
				vec.AddTo(eOut.Coeffs[i], p.Coeffs[i], num.MForm(c.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				vec.SubTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.Form = p.Form
	}
}

func (op *primeAutFixedAddSubOperator) subTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case e0.Type == TypeScalar && e1.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}

	case e0.Type == TypePoly && e1.Type == TypePoly:
		for i := range op.mod {
			vec.SubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.Form = e0.Form

	case e0.Type == TypeScalar && e1.Type == TypePoly:
		for i := range op.mod {
			if e1.Form == dft.FormNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOut.Coeffs[i], num.MForm(e0.Coeffs[i][0], op.mod[i]), e1.Coeffs[i], op.mod[i])
			} else {
				vec.SubTo(eOut.Coeffs[i], num.Neg(e0.Coeffs[i][0], op.mod[i]), e1.Coeffs[i], op.mod[i])
			}
		}
		eOut.Form = e1.Form

	case e0.Type == TypePoly && e1.Type == TypeScalar:
		for i := range op.mod {
			if e0.Form == dft.FormNTT && op.isNTTFriendly[i] {
				vec.SubTo(eOut.Coeffs[i], e0.Coeffs[i], num.MForm(e1.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				vec.AddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.Form = e0.Form
	}
}

func (op *primeAutFixedAddSubOperator) slice(lo, hi int) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],
	}
}

func (op *primeAutFixedAddSubOperator) append(opIn *primeAutFixedAddSubOperator) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          opIn.rank,
		mod:           slices.Concat(op.mod, opIn.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, opIn.isNTTFriendly),
	}
}

func (op *primeAutFixedAddSubOperator) appendTmpMod(mod *num.Modulus) *primeAutFixedAddSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),
	}
}
