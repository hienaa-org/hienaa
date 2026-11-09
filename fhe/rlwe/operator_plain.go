package rlwe

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// PlainOperator evaluates ring operations over scalar and polynomial [Element] values.
//
// Operations panic when inputs are not consistent with the operator's parameters,
// or the operations are not valid. Arithmetic operations require matching base
// and auxiliary modulus lengths.
type PlainOperator struct {
	params Parameters
	crtOp  *crt.Operator

	poolEmb *pool.Pool[*[]uint64]
}

// NewPlainOperator creates a new [PlainOperator].
func NewPlainOperator(params Parameters) *PlainOperator {
	return &PlainOperator{
		params: params,
		crtOp:  params.crtOp,

		poolEmb: crt.NewEmbedderPool(params.Rank(), params.FullModLen()),
	}
}

// Params returns the [Parameters] of the operator.
func (op *PlainOperator) Params() Parameters {
	return op.params
}

// CRTOp returns the [crt.Operator] of the operator.
func (op *PlainOperator) CRTOp() *crt.Operator {
	return op.crtOp
}

// NewScalar creates a new scalar [Element].
func (op *PlainOperator) NewScalar(withAux bool) *Element {
	return NewScalar(op.params, withAux)
}

// NewPoly creates a new polynomial [Element].
func (op *PlainOperator) NewPoly(withAux bool, form dft.Form) *Element {
	return NewPoly(op.params, withAux, form)
}

// withModLen returns a [crt.Operator] for given modulus length.
func (op *PlainOperator) withModLen(baseLen, auxLen int) *crt.Operator {
	return op.crtOp.WithModSlice(op.params.AuxModLen()-auxLen, op.params.AuxModLen()+baseLen)
}

// FwdNTT returns FwdNTT(p).
//
// Panics if p is [crt.TypeScalar].
func (op *PlainOperator) FwdNTT(p *Element) *Element {
	pOut := NewPolyCustom(op.params.Rank(), p.BaseModLen(), p.AuxModLen(), dft.FormNTT)
	op.FwdNTTTo(pOut, p)
	return pOut
}

// FwdNTTTo computes pOut = NTT(p).
//
// Panics if p is [crt.TypeScalar].
func (op *PlainOperator) FwdNTTTo(pOut, p *Element) {
	checkUnaryOperable(op.params, pOut, p)
	op.withModLen(pOut.BaseModLen(), pOut.AuxModLen()).FwdNTTTo(pOut.Value, p.Value)
}

// InvNTT returns InvNTT(p).
//
// Panics if p is [crt.TypeScalar].
func (op *PlainOperator) InvNTT(p *Element) *Element {
	pOut := NewPolyCustom(op.params.Rank(), p.BaseModLen(), p.AuxModLen(), dft.FormCoeff)
	op.InvNTTTo(pOut, p)
	return pOut
}

// InvNTTTo computes pOut = InvNTT(p).
//
// Panics if p is [crt.TypeScalar].
func (op *PlainOperator) InvNTTTo(pOut, p *Element) {
	checkUnaryOperable(op.params, pOut, p)
	op.withModLen(pOut.BaseModLen(), pOut.AuxModLen()).InvNTTTo(pOut.Value, p.Value)
}

// Add returns e0 + e1.
func (op *PlainOperator) Add(e0, e1 *Element) *Element {
	eOut := newBinaryOpOut(e0, e1)
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *PlainOperator) AddTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params, eOut, e0, e1)
	op.withModLen(eOut.BaseModLen(), eOut.AuxModLen()).AddTo(eOut.Value, e0.Value, e1.Value)
}

// Sub returns e0 - e1.
func (op *PlainOperator) Sub(e0, e1 *Element) *Element {
	eOut := newBinaryOpOut(e0, e1)
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *PlainOperator) SubTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params, eOut, e0, e1)
	op.withModLen(eOut.BaseModLen(), eOut.AuxModLen()).SubTo(eOut.Value, e0.Value, e1.Value)
}

// Neg returns -e.
func (op *PlainOperator) Neg(e *Element) *Element {
	eOut := newUnaryOpOut(e)
	op.NegTo(eOut, e)
	return eOut
}

// NegTo computes eOut = -e.
func (op *PlainOperator) NegTo(eOut, e *Element) {
	checkUnaryOperable(op.params, eOut, e)
	op.withModLen(eOut.BaseModLen(), eOut.AuxModLen()).NegTo(eOut.Value, e.Value)
}

// Mul returns e0 * e1.
func (op *PlainOperator) Mul(e0, e1 *Element) *Element {
	eOut := newBinaryOpOut(e0, e1)
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
func (op *PlainOperator) MulTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params, eOut, e0, e1)
	op.withModLen(eOut.BaseModLen(), eOut.AuxModLen()).MulTo(eOut.Value, e0.Value, e1.Value)
}

// MulAddTo computes eOut += e0 * e1.
func (op *PlainOperator) MulAddTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params, eOut, e0, e1)
	op.withModLen(eOut.BaseModLen(), eOut.AuxModLen()).MulAddTo(eOut.Value, e0.Value, e1.Value)
}

// MulSubTo computes eOut -= e0 * e1.
func (op *PlainOperator) MulSubTo(eOut, e0, e1 *Element) {
	checkBinaryOperable(op.params, eOut, e0, e1)
	op.withModLen(eOut.BaseModLen(), eOut.AuxModLen()).MulSubTo(eOut.Value, e0.Value, e1.Value)
}

// CanAut returns if automorphism X -> X^idx is valid.
func (op *PlainOperator) CanAut(idx int) bool {
	return op.crtOp.CanAut(idx)
}

// Aut returns p(X^idx).
//
// Panics if p is [crt.TypeScalar].
func (op *PlainOperator) Aut(p *Element, idx int) *Element {
	pOut := NewPolyCustom(op.params.Rank(), p.BaseModLen(), p.AuxModLen(), dft.FormCoeff)
	op.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = p(X^idx).
//
// Panics if p is [crt.TypeScalar].
func (op *PlainOperator) AutTo(pOut, p *Element, idx int) {
	checkUnaryOperable(op.params, pOut, p)
	op.withModLen(pOut.BaseModLen(), pOut.AuxModLen()).AutTo(pOut.Value, p.Value, idx)
}

// Embed returns an embedding of e to modulus length outBaseLen and outAuxLen.
func (op *PlainOperator) Embed(e *Element, outBaseLen, outAuxLen int, outForm dft.Form) *Element {
	eOut := NewElementCustom(e.Rank(), outBaseLen, outAuxLen, e.Type(), outForm)
	op.EmbedTo(eOut, e, outForm)
	return eOut
}

// EmbedTo embeds e to eOut.
func (op *PlainOperator) EmbedTo(eOut, e *Element, outForm dft.Form) {
	checkElementShape(op.params, eOut)
	checkElementShape(op.params, e)

	opIn := op.withModLen(e.BaseModLen(), e.AuxModLen())
	opOut := op.withModLen(eOut.BaseModLen(), eOut.AuxModLen())

	emb := crt.NewEmbedderUnsafe(opOut, opIn, op.poolEmb)
	emb.EmbedTo(eOut.Value, e.Value, outForm)
}

// Scale returns a scaling of e to modulus length outBaseLen and outAuxLen.
//
// Panics when e has auxiliary modulus.
func (op *PlainOperator) Scale(e *Element, outBaseLen, outAuxLen int, outForm dft.Form) *Element {
	eOut := NewElementCustom(e.Rank(), outBaseLen, outAuxLen, e.Type(), outForm)
	op.ScaleTo(eOut, e, outForm)
	return eOut
}

// ScaleTo scales e to eOut.
func (op *PlainOperator) ScaleTo(eOut, e *Element, outForm dft.Form) {
	checkElementShape(op.params, eOut)
	checkElementShape(op.params, e)

	opIn := op.withModLen(e.BaseModLen(), e.AuxModLen())
	opOut := op.withModLen(eOut.BaseModLen(), eOut.AuxModLen())

	sc := crt.NewScalerUnsafe(opOut, opIn, op.poolEmb)
	sc.ScaleTo(eOut.Value, e.Value, outForm)
}
