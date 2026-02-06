package crt

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// baseOperator is the base operator for all rings.
type baseOperator struct {
	params  dft.RingParameters
	mod     []*num.Modulus
	modPoly []int64

	isNTTFriendly []bool
	ntt           []dft.Transformer
}

// newBaseOperator creates a new [baseOperator].
func newBaseOperator(params dft.RingParameters, mod []*num.Modulus, modPoly []int64) baseOperator {
	isNTTFriendly := make([]bool, len(mod))
	ntt := make([]dft.Transformer, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
		if isNTTFriendly[i] {
			ntt[i] = dft.NewTransformer(params, mod[i])
		}
	}

	return baseOperator{
		params:  params,
		mod:     mod,
		modPoly: modPoly,

		isNTTFriendly: isNTTFriendly,
		ntt:           ntt,
	}
}

// Params returns the ring parameters.
func (op *baseOperator) Params() dft.RingParameters {
	return op.params
}

// Modulus returns the modulus.
func (op *baseOperator) Modulus() []*num.Modulus {
	return op.mod
}

// ModulusPoly returns the quotient polynomial of the ring.
func (op *baseOperator) ModulusPoly() []int64 {
	return op.modPoly
}

// NewPoly creates a new polynomial element.
func (op *baseOperator) NewPoly() *Element {
	return NewPoly(op.params.Rank(), len(op.mod))
}

// NewNTTPoly creates a new polynomial element in NTT form.
func (op *baseOperator) NewNTTPoly() *Element {
	return NewNTTPoly(op.params.Rank(), len(op.mod))
}

// NewPolyCustom creates a new polynomial element.
func (op *baseOperator) NewPolyCustom(isNTT bool) *Element {
	return NewPolyCustom(op.params.Rank(), len(op.mod), isNTT)
}

// FwdNTT returns FwdNTT(e).
func (op *baseOperator) FwdNTT(e *Element) *Element {
	eOut := NewPoly(e.Rank(), len(op.mod))
	op.FwdNTTTo(eOut, e)
	return eOut
}

// FwdNTTTo computes eOut = NTT(e).
func (op *baseOperator) FwdNTTTo(eOut, e *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	if e.Type() == TypeScalar {
		eOut.CopyFrom(e)
		return
	}

	if e.IsNTT {
		panic("input(s) must be in standard form")
	}

	for i := range op.ntt {
		if op.ntt[i] != nil {
			op.ntt[i].ForwardTo(eOut.Coeffs[i], e.Coeffs[i])
		} else {
			copy(eOut.Coeffs[i], e.Coeffs[i])
		}
	}

	eOut.IsNTT = true
}

// InvNTT returns InvNTT(e).
func (op *baseOperator) InvNTT(e *Element) *Element {
	eOut := NewPoly(e.Rank(), len(op.mod))
	op.InvNTTTo(eOut, e)
	return eOut
}

// InvNTTTo computes eOut = InvNTT(e).
func (op *baseOperator) InvNTTTo(eOut, e *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	if e.Type() == TypeScalar {
		eOut.CopyFrom(e)
		return
	}

	if !e.IsNTT {
		panic("input(s) must be in NTT form")
	}

	for i := range op.ntt {
		if op.ntt[i] != nil {
			op.ntt[i].InverseTo(eOut.Coeffs[i], e.Coeffs[i])
		} else {
			copy(eOut.Coeffs[i], e.Coeffs[i])
		}
	}

	eOut.IsNTT = false
}

// Neg returns -e.
func (op *baseOperator) Neg(e *Element) *Element {
	eOut := NewPolyCustom(e.Rank(), len(op.mod), e.IsNTT)
	op.NegTo(eOut, e)
	return eOut
}

// NegTo computes eOut = -e.
func (op *baseOperator) NegTo(eOut, e *Element) {
	checkBinaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	switch e.Type() {
	case TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Neg(e.Coeffs[i][0], op.mod[i])
		}
	case TypePoly:
		for i := range op.mod {
			vec.NegTo(eOut.Coeffs[i], e.Coeffs[i], op.mod[i])
		}
		eOut.IsNTT = e.IsNTT
	}
}

// AsBig returns p as *[big.Int] vector.
func (op *baseOperator) AsBig(e *Element) []*big.Int {
	if e.Type() == TypePoly {
		checkShape(op.params.Rank(), len(op.mod), e)
		if e.IsNTT {
			panic("input(s) must be in standard form")
		}
	}

	modBig := make([]*big.Int, len(op.mod))
	modProd := big.NewInt(1)
	for i := range op.mod {
		modBig[i] = new(big.Int).SetUint64(op.mod[i].Value())
		modProd.Mul(modProd, modBig[i])
	}
	modProdHalf := new(big.Int).Rsh(modProd, 1)

	gadget := make([]*big.Int, len(op.mod))
	for i := range modBig {
		qStar := new(big.Int).Div(modProd, modBig[i])
		qStarInv := new(big.Int).ModInverse(qStar, modBig[i])
		gadget[i] = new(big.Int).Mul(qStar, qStarInv)
		gadget[i].Mod(gadget[i], modProd)
	}

	pBig := make([]*big.Int, e.Rank())
	for j := 0; j < e.Rank(); j++ {
		pBig[j] = big.NewInt(0)
		for i := range op.mod {
			c := new(big.Int).SetUint64(e.Coeffs[i][j])
			c.Mul(c, gadget[i])
			pBig[j].Add(pBig[j], c)
		}
		pBig[j].Mod(pBig[j], modProd)
		if pBig[j].Cmp(modProdHalf) > 0 {
			pBig[j].Sub(pBig[j], modProd)
		}
	}

	return pBig
}

func (op *baseOperator) subOperator(idx ...int) baseOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	nttCopy := make([]dft.Transformer, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
		if op.ntt[idx[i]] != nil {
			nttCopy[i] = op.ntt[idx[i]].SafeCopy()
		}
	}

	return baseOperator{
		params:  op.params,
		mod:     modCopy,
		modPoly: op.modPoly,

		isNTTFriendly: isNTTFriendlyCopy,
		ntt:           nttCopy,
	}
}

func (op *baseOperator) safeCopy() baseOperator {
	nttCopy := make([]dft.Transformer, len(op.ntt))
	for i := range op.ntt {
		if op.ntt[i] != nil {
			nttCopy[i] = op.ntt[i].SafeCopy()
		}
	}

	return baseOperator{
		params:  op.params,
		mod:     op.mod,
		modPoly: op.modPoly,

		isNTTFriendly: op.isNTTFriendly,
		ntt:           nttCopy,
	}
}
