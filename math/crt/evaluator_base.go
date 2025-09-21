package crt

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// polyEvaluatorBase is the base evaluator for all rings.
type polyEvaluatorBase struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool
	ntt           []dft.Transformer
}

// newPolyEvaluatorBase creates a new [polyEvaluatorBase].
func newPolyEvaluatorBase(params dft.RingParameters, mod []*num.Modulus) *polyEvaluatorBase {
	isNTTFriendly := make([]bool, len(mod))
	ntt := make([]dft.Transformer, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
		if isNTTFriendly[i] {
			ntt[i] = dft.NewTransformer(params, mod[i])
		}
	}

	return &polyEvaluatorBase{
		params: params,
		mod:    mod,

		isNTTFriendly: isNTTFriendly,
		ntt:           ntt,
	}
}

// NewPoly returns a new [Poly].
func (e *polyEvaluatorBase) NewPoly() *Poly {
	return NewPolyCustom(e.params.Rank(), len(e.mod), false)
}

// NewNTTPoly returns a new [Poly] in NTT form.
func (e *polyEvaluatorBase) NewNTTPoly() *Poly {
	return NewPolyCustom(e.params.Rank(), len(e.mod), true)
}

// NewPolyCustom creates a new [Poly] with the given parameters.
func (e *polyEvaluatorBase) NewPolyCustom(isNTT bool) *Poly {
	return NewPolyCustom(e.params.Rank(), len(e.mod), isNTT)
}

// Params returns the ring parameters.
func (e *polyEvaluatorBase) Params() dft.RingParameters {
	return e.params
}

// Modulus returns the modulus.
func (e *polyEvaluatorBase) Modulus() []*num.Modulus {
	return e.mod
}

// Add returns p0 + p1.
func (e *polyEvaluatorBase) Add(p0, p1 *Poly) *Poly {
	pOut := e.NewPoly()
	e.AddTo(pOut, p0, p1)
	return pOut
}

// AddTo computes pOut = p0 + p1.
func (e *polyEvaluatorBase) AddTo(pOut, p0, p1 *Poly) {
	if !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1) {
		panic("AddTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.AddTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
	}

	pOut.isNTT = p0.isNTT
}

// Sub returns p0 - p1.
func (e *polyEvaluatorBase) Sub(p0, p1 *Poly) *Poly {
	pOut := e.NewPoly()
	e.SubTo(pOut, p0, p1)
	return pOut
}

// SubTo computes pOut = p0 - p1.
func (e *polyEvaluatorBase) SubTo(pOut, p0, p1 *Poly) {
	if !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1) {
		panic("SubTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.SubTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
	}

	pOut.isNTT = p0.isNTT
}

// Neg returns -p.
func (e *polyEvaluatorBase) Neg(p *Poly) *Poly {
	pOut := NewPolyCustom(e.params.Rank(), len(e.mod), p.isNTT)
	e.NegTo(pOut, p)
	return pOut
}

// NegTo computes pOut = -p.
func (e *polyEvaluatorBase) NegTo(pOut, p *Poly) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("NegTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.NegTo(pOut.Coeffs[i], p.Coeffs[i], e.mod[i])
	}

	pOut.isNTT = p.isNTT
}

// ScalarMul returns p * c.
func (e *polyEvaluatorBase) ScalarMul(p *Poly, c Scalar) *Poly {
	pOut := e.NewPoly()
	e.ScalarMulTo(pOut, p, c)
	return pOut
}

// ScalarMulTo computes pOut = p * c.
func (e *polyEvaluatorBase) ScalarMulTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarMulTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.ScalarMulTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
	}

	pOut.isNTT = p.isNTT
}

// ScalarMulAddTo computes pOut += p * c.
func (e *polyEvaluatorBase) ScalarMulAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarMulAddTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.ScalarMulAddTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
	}

	pOut.isNTT = p.isNTT
}

// ScalarMulSubTo computes pOut -= p * c.
func (e *polyEvaluatorBase) ScalarMulSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarMulAddTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.ScalarMulSubTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
	}

	pOut.isNTT = p.isNTT
}

// NTT returns NTT(p).
func (e *polyEvaluatorBase) NTT(p *Poly) *Poly {
	pOut := e.NewPoly()
	e.NTTTo(pOut, p)
	return pOut
}

// NTTTo computes pOut = NTT(p).
func (e *polyEvaluatorBase) NTTTo(pOut, p *Poly) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("NTTTo: inputs not consistent")
	case p.isNTT:
		panic("NTTTo: already in NTT form")
	}

	for i := range e.ntt {
		if e.ntt[i] != nil {
			e.ntt[i].ForwardTo(pOut.Coeffs[i], p.Coeffs[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
		}
	}

	pOut.isNTT = true
}

// InvNTT returns InvNTT(p).
func (e *polyEvaluatorBase) InvNTT(p *Poly) *Poly {
	pOut := e.NewPoly()
	e.InvNTTTo(pOut, p)
	return pOut
}

// InvNTTTo computes pOut = InvNTT(p).
func (e *polyEvaluatorBase) InvNTTTo(pOut, p *Poly) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("InvNTTTo: inputs not consistent")
	case !p.isNTT:
		panic("InvNTTTo: already in Standard form")
	}

	for i := range e.ntt {
		if e.ntt[i] != nil {
			e.ntt[i].InverseTo(pOut.Coeffs[i], p.Coeffs[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
		}
	}

	pOut.isNTT = false
}

func (e *polyEvaluatorBase) subEvaluator(idx ...int) *polyEvaluatorBase {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	nttCopy := make([]dft.Transformer, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
		if e.ntt[idx[i]] != nil {
			nttCopy[i] = e.ntt[idx[i]].SafeCopy()
		}
	}

	return &polyEvaluatorBase{
		params: e.params,
		mod:    modCopy,

		isNTTFriendly: isNTTFriendlyCopy,
		ntt:           nttCopy,
	}
}

func (e *polyEvaluatorBase) safeCopy() *polyEvaluatorBase {
	nttCopy := make([]dft.Transformer, len(e.ntt))
	for i := range e.ntt {
		if e.ntt[i] != nil {
			nttCopy[i] = e.ntt[i].SafeCopy()
		}
	}

	return &polyEvaluatorBase{
		params: e.params,
		mod:    e.mod,

		isNTTFriendly: e.isNTTFriendly,
		ntt:           nttCopy,
	}
}
