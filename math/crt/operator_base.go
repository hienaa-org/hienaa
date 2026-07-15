package crt

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// baseOperator is the base operator for all rings.
type baseOperator struct {
	params dft.RingParameters
	mod    []*num.Modulus

	isNTTFriendly []bool
	ntt           []dft.Transformer
}

// newBaseOperator creates a new [baseOperator].
func newBaseOperator(params dft.RingParameters, mod []*num.Modulus) *baseOperator {
	isNTTFriendly := make([]bool, len(mod))
	ntt := make([]dft.Transformer, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
		if isNTTFriendly[i] {
			ntt[i] = dft.NewTransformer(params, mod[i])
		}
	}

	return &baseOperator{
		params: params,
		mod:    mod,

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
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

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
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

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
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

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

func (op *baseOperator) withModIdx(idx ...int) *baseOperator {
	return &baseOperator{
		params: op.params,
		mod:    vec.Gather(op.mod, idx...),

		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),
		ntt:           vec.Gather(op.ntt, idx...),
	}
}

func (op *baseOperator) slice(lo, hi int) *baseOperator {
	return &baseOperator{
		params: op.params,
		mod:    op.mod[lo:hi:hi],

		isNTTFriendly: op.isNTTFriendly[lo:hi:hi],
		ntt:           op.ntt[lo:hi:hi],
	}
}

func (op *baseOperator) append(op0 *baseOperator) *baseOperator {
	return &baseOperator{
		params: op.params,
		mod:    vec.Concat(op.mod, op0.mod),

		isNTTFriendly: vec.Concat(op.isNTTFriendly, op0.isNTTFriendly),
		ntt:           vec.Concat(op.ntt, op0.ntt),
	}
}

func (op *baseOperator) appendTmpModulus(mod *num.Modulus) *baseOperator {
	return &baseOperator{
		params: op.params,
		mod:    vec.Concat(op.mod, []*num.Modulus{mod}),

		isNTTFriendly: vec.Concat(op.isNTTFriendly, []bool{false}),
		ntt:           vec.Concat(op.ntt, []dft.Transformer{nil}),
	}
}
