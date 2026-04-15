package crt

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// addSubOperator implements addition and subtraction operations.
type addSubOperator interface {
	// Add returns e0 + e1.
	Add(e0, e1 *Element) *Element
	// AddTo computes eOut = e0 + e1.
	AddTo(eOut, e0, e1 *Element)

	// Sub returns e0 - e1.
	Sub(e0, e1 *Element) *Element
	// SubTo computes eOut = e0 - e1.
	SubTo(eOut, e0, e1 *Element)

	subOperator(idx ...int) addSubOperator
	append(op0 addSubOperator) addSubOperator
	appendAuxModulus(mod *num.Modulus) addSubOperator
}

// baseAddSubOperator is a [addSubOperator] for every ring
// except for prime-order autfixed ring.
type baseAddSubOperator struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newBaseAddSubOperator creates a new [baseAddSubOperator].
func newBaseAddSubOperator(params dft.RingParameters, mod []*num.Modulus) *baseAddSubOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &baseAddSubOperator{
		rank:          params.Rank(),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

// Add returns e0 + e1.
func (op *baseAddSubOperator) Add(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *baseAddSubOperator) AddTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		for i := range op.mod {
			vec.AddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.IsNTT = e0.IsNTT
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			if p.IsNTT && op.isNTTFriendly[i] {
				vec.AddScalarTo(eOut.Coeffs[i], p.Coeffs[i], num.MForm(c.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				copy(eOut.Coeffs[i], p.Coeffs[i])
				eOut.Coeffs[i][0] = num.Add(eOut.Coeffs[i][0], c.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.IsNTT = p.IsNTT
	}
}

// Sub returns e0 - e1.
func (op *baseAddSubOperator) Sub(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *baseAddSubOperator) SubTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		for i := range op.mod {
			vec.SubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.IsNTT = e0.IsNTT
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			if p.IsNTT && op.isNTTFriendly[i] {
				vec.SubScalarTo(eOut.Coeffs[i], p.Coeffs[i], num.MForm(c.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				copy(eOut.Coeffs[i], p.Coeffs[i])
				eOut.Coeffs[i][0] = num.Sub(eOut.Coeffs[i][0], c.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.IsNTT = p.IsNTT
	}
}

func (op *baseAddSubOperator) subOperator(idx ...int) addSubOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
	}

	return &baseAddSubOperator{
		rank:          op.rank,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,
	}
}

func (op *baseAddSubOperator) append(op0 addSubOperator) addSubOperator {
	opOther := op0.(*baseAddSubOperator)
	return &baseAddSubOperator{
		rank:          opOther.rank,
		mod:           vec.Concat(op.mod, opOther.mod),
		isNTTFriendly: vec.Concat(op.isNTTFriendly, opOther.isNTTFriendly),
	}
}

func (op *baseAddSubOperator) appendAuxModulus(mod *num.Modulus) addSubOperator {
	return &baseAddSubOperator{
		rank:          op.rank,
		mod:           vec.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: vec.Concat(op.isNTTFriendly, []bool{false}),
	}
}

// primeAutFixedAddSubOperator is a [addSubOperator] for prime-order autfixed ring.
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

// Add returns e0 + e1.
func (op *primeAutFixedAddSubOperator) Add(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *primeAutFixedAddSubOperator) AddTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Add(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		for i := range op.mod {
			vec.AddTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.IsNTT = e0.IsNTT
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			if p.IsNTT && op.isNTTFriendly[i] {
				vec.AddScalarTo(eOut.Coeffs[i], p.Coeffs[i], num.MForm(c.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				vec.SubScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.IsNTT = p.IsNTT
	}
}

// Sub returns e0 - e1.
func (op *primeAutFixedAddSubOperator) Sub(e0, e1 *Element) *Element {
	eOut := NewPoly(max(e0.Rank(), e1.Rank()), len(op.mod))
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *primeAutFixedAddSubOperator) SubTo(eOut, e0, e1 *Element) {
	isBinaryOperable(op.rank, len(op.mod), eOut, e0, e1)

	switch {
	case isEqualType(e0, e1, TypeScalar):
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Sub(e0.Coeffs[i][0], e1.Coeffs[i][0], op.mod[i])
		}
	case isEqualType(e0, e1, TypePoly):
		for i := range op.mod {
			vec.SubTo(eOut.Coeffs[i], e0.Coeffs[i], e1.Coeffs[i], op.mod[i])
		}
		eOut.IsNTT = e0.IsNTT
	default:
		c, p := orderByType(e0, e1)
		for i := range op.mod {
			if p.IsNTT && op.isNTTFriendly[i] {
				vec.SubScalarTo(eOut.Coeffs[i], p.Coeffs[i], num.MForm(c.Coeffs[i][0], op.mod[i]), op.mod[i])
			} else {
				vec.AddScalarTo(eOut.Coeffs[i], p.Coeffs[i], c.Coeffs[i][0], op.mod[i])
			}
		}
		eOut.IsNTT = p.IsNTT
	}
}

func (op *primeAutFixedAddSubOperator) subOperator(idx ...int) addSubOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
	}

	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,
	}
}

func (op *primeAutFixedAddSubOperator) append(op0 addSubOperator) addSubOperator {
	opOther := op0.(*primeAutFixedAddSubOperator)
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           vec.Concat(op.mod, opOther.mod),
		isNTTFriendly: vec.Concat(op.isNTTFriendly, opOther.isNTTFriendly),
	}
}

func (op *primeAutFixedAddSubOperator) appendAuxModulus(mod *num.Modulus) addSubOperator {
	return &primeAutFixedAddSubOperator{
		rank:          op.rank,
		mod:           vec.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: vec.Concat(op.isNTTFriendly, []bool{false}),
	}
}
