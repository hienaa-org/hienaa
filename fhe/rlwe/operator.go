package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

//////////////////////////////////////////
///////// PlainScalar operations /////////
//////////////////////////////////////////

type PlainOperator struct {
	Params Parameters
	Eval   crt.PolyEvaluator
}

func NewPlainOperator(p Parameters) *PlainOperator {
	eval := crt.NewPolyEvaluator(p.ringParams, append(p.auxModulus, p.modulus...))

	return &PlainOperator{
		Params: p,
		Eval:   eval,
	}
}

// SafeCopy returns a thread-safe copy.
func (o *PlainOperator) SafeCopy() *PlainOperator {
	eval := o.Eval.SafeCopy()
	return &PlainOperator{
		Params: o.Params,
		Eval:   eval,
	}
}

// SubEvaluatorAt returns a sub evaluator at the given indices.
func (o *PlainOperator) SubEvaluatorAt(hasAux bool, modLen int) crt.PolyEvaluator {
	if hasAux {
		if o.Params.auxModulus == nil {
			panic("auxiliary modulus is not set")
		} else if modLen > len(o.Params.modulus)+len(o.Params.auxModulus) {
			panic("modLen is too large")
		}
		return o.Eval.SubEvaluator(vec.Range(0, modLen)...)
	} else {
		if modLen > len(o.Params.modulus) {
			panic("modLen is too large")
		}
		return o.Eval.SubEvaluator(vec.Range(len(o.Params.auxModulus), modLen+len(o.Params.auxModulus))...)
	}
}

// AddScalar returns c0 + c1.
func (o *PlainOperator) AddScalar(c0 *PlainScalar, c1 *PlainScalar) *PlainScalar {
	modLen := c0.ModLen()
	var auxLen int
	if c0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewPlainScalarCustom(modLen, auxLen)
	o.AddScalarTo(cOut, c0, c1)

	return cOut
}

// AddScalarTo performs c0 + c1 and stores the result in cOut.
func (o *PlainOperator) AddScalarTo(cOut *PlainScalar, c0 *PlainScalar, c1 *PlainScalar) {
	if !(cOut.IsConsistent(c0) && cOut.IsConsistent(c1)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	crt.AddScalarTo(cOut.Value, c0.Value, c1.Value, eval.Modulus())
	cOut.HasAux = c0.HasAux
}

// SubScalar returns c0 - c1.
func (o *PlainOperator) SubScalar(c0 *PlainScalar, c1 *PlainScalar) *PlainScalar {
	modLen := c0.ModLen()
	var auxLen int
	if c0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewPlainScalarCustom(modLen, auxLen)
	o.SubScalarTo(cOut, c0, c1)

	return cOut
}

// SubScalarTo performs c0 - c1 and stores the result in cOut.
func (o *PlainOperator) SubScalarTo(cOut *PlainScalar, c0 *PlainScalar, c1 *PlainScalar) {
	if !(cOut.IsConsistent(c0) && cOut.IsConsistent(c1)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	crt.SubScalarTo(cOut.Value, c0.Value, c1.Value, eval.Modulus())
	cOut.HasAux = c0.HasAux
}

// NegScalar returns -c.
func (o *PlainOperator) NegScalar(c *PlainScalar) *PlainScalar {
	modLen := c.ModLen()
	var auxLen int
	if c.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewPlainScalarCustom(modLen, auxLen)
	o.NegScalarTo(cOut, c)

	return cOut
}

// NegScalarTo performs -c and stores the result in cOut.
func (o *PlainOperator) NegScalarTo(cOut *PlainScalar, c *PlainScalar) {
	if !cOut.IsConsistent(c) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	crt.NegScalarTo(cOut.Value, c.Value, eval.Modulus())
	cOut.HasAux = c.HasAux
}

// MulScalar returns c0 * c1.
func (o *PlainOperator) MulScalar(c0 *PlainScalar, c1 *PlainScalar) *PlainScalar {
	modLen := c0.ModLen()
	var auxLen int
	if c0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewPlainScalarCustom(modLen, auxLen)
	o.MulScalarTo(cOut, c0, c1)

	return cOut
}

// MulScalarTo performs c0 * c1 and stores the result in cOut.
func (o *PlainOperator) MulScalarTo(cOut *PlainScalar, c0 *PlainScalar, c1 *PlainScalar) {
	if !(cOut.IsConsistent(c0) && cOut.IsConsistent(c1)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	crt.MulScalarTo(cOut.Value, c0.Value, c1.Value, eval.Modulus())
	cOut.HasAux = c0.HasAux
}

// MulAddScalarTo performs cOut += c0 * c1.
func (o *PlainOperator) MulAddScalarTo(cOut *PlainScalar, c0 *PlainScalar, c1 *PlainScalar) {
	if !(cOut.IsConsistent(c0) && cOut.IsConsistent(c1)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	crt.MulAddScalarTo(cOut.Value, c0.Value, c1.Value, eval.Modulus())
	cOut.HasAux = c0.HasAux
}

// MulSubScalarTo performs cOut -= c0 * c1.
func (o *PlainOperator) MulSubScalarTo(cOut *PlainScalar, c0 *PlainScalar, c1 *PlainScalar) {
	if !(cOut.IsConsistent(c0) && cOut.IsConsistent(c1)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	crt.MulSubScalarTo(cOut.Value, c0.Value, c1.Value, eval.Modulus())
	cOut.HasAux = c0.HasAux
}

////////////////////////////////////////
///////// PlainPoly operations /////////
////////////////////////////////////////

// FwdNTT performs FwdNTT(pt).
func (o *PlainOperator) FwdNTT(pt *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := pt.ModLen()
	var auxLen int
	if pt.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, true)
	o.FwdNTTTo(ptOut, pt)

	return ptOut
}

// FwdNTTTo performs FwdNTT(pt) and stores the result in pOut.
func (o *PlainOperator) FwdNTTTo(ptOut *PlainPoly, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn)) {
		panic("FwdNTTTo: inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.FwdNTTTo(ptOut.Value, ptIn.Value)
}

// InvNTT performs InvNTT(pt).
func (o *PlainOperator) InvNTT(pt *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := pt.ModLen()
	var auxLen int
	if pt.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.InvNTTTo(ptOut, pt)

	return ptOut
}

// InvNTTTo performs InvNTT(pt) and stores the result in pOut.
func (o *PlainOperator) InvNTTTo(ptOut *PlainPoly, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.InvNTTTo(ptOut.Value, ptIn.Value)
}

// ScalarAdd returns c + pt.
func (o *PlainOperator) ScalarAdd(cIn *PlainScalar, ptIn *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := cIn.ModLen()
	var auxLen int
	if cIn.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.ScalarAddTo(ptOut, cIn, ptIn)

	return ptOut
}

// ScalarAddTo performs c + pt and stores the result in ptOut.
func (o *PlainOperator) ScalarAddTo(ptOut *PlainPoly, cIn *PlainScalar, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn) && (ptOut.ModLen() == cIn.ModLen())) {
		panic("inconsistent plain polynomials or scalars")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.ScalarAddTo(ptOut.Value, ptIn.Value, cIn.Value)
}

// ScalarSub returns c - pt.
func (o *PlainOperator) ScalarSub(cIn *PlainScalar, ptIn *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := cIn.ModLen()
	var auxLen int
	if cIn.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.ScalarSubTo(ptOut, cIn, ptIn)

	return ptOut
}

// ScalarSubTo performs c - pt and stores the result in ptOut.
func (o *PlainOperator) ScalarSubTo(ptOut *PlainPoly, cIn *PlainScalar, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn) && (ptOut.ModLen() == cIn.ModLen())) {
		panic("inconsistent plain polynomials or scalars")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.ScalarSubTo(ptOut.Value, ptIn.Value, cIn.Value)
}

// ScalarMul returns c * pt.
func (o *PlainOperator) ScalarMul(cIn *PlainScalar, ptIn *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := cIn.ModLen()
	var auxLen int
	if cIn.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.ScalarMulTo(ptOut, cIn, ptIn)

	return ptOut
}

// ScalarMulTo performs c * pt and stores the result in ptOut.
func (o *PlainOperator) ScalarMulTo(ptOut *PlainPoly, cIn *PlainScalar, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn) && (ptOut.ModLen() == cIn.ModLen())) {
		panic("inconsistent plain polynomials or scalars")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.ScalarMulTo(ptOut.Value, ptIn.Value, cIn.Value)
}

// ScalarMulAddTo performs ptOut += c * pt.
func (o *PlainOperator) ScalarMulAddTo(ptOut *PlainPoly, cIn *PlainScalar, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn) && (ptOut.ModLen() == cIn.ModLen())) {
		panic("inconsistent plain polynomials or scalars")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.ScalarMulAddTo(ptOut.Value, ptIn.Value, cIn.Value)
}

// ScalarMulSubTo performs ptOut -= c * pt.
func (o *PlainOperator) ScalarMulSubTo(ptOut *PlainPoly, cIn *PlainScalar, ptIn *PlainPoly) {
	if !(ptOut.IsConsistent(ptIn) && (ptOut.ModLen() == cIn.ModLen())) {
		panic("inconsistent plain polynomials or scalars")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.ScalarMulSubTo(ptOut.Value, ptIn.Value, cIn.Value)
}

// Add performs pt0 + pt1.
func (o *PlainOperator) Add(pt0, pt1 *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := pt0.ModLen()
	var auxLen int
	if pt0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.AddTo(ptOut, pt0, pt1)

	return ptOut
}

// AddTo performs pt0 + pt1 and stores the result in ptOut.
func (o *PlainOperator) AddTo(ptOut *PlainPoly, pt0, pt1 *PlainPoly) {
	if !(ptOut.IsConsistent(pt0) && ptOut.IsConsistent(pt1)) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(pt0.HasAux, pt0.ModLen())
	eval.AddTo(ptOut.Value, pt0.Value, pt1.Value)
}

// Sub performs pt0 - pt1.
func (o *PlainOperator) Sub(pt0, pt1 *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := pt0.ModLen()
	var auxLen int
	if pt0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.SubTo(ptOut, pt0, pt1)

	return ptOut
}

// SubTo performs pt0 - pt1 and stores the result in ptOut.
func (o *PlainOperator) SubTo(ptOut *PlainPoly, pt0, pt1 *PlainPoly) {
	if !(ptOut.IsConsistent(pt0) && ptOut.IsConsistent(pt1)) {
		panic("SubTo: inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(pt0.HasAux, pt0.ModLen())
	eval.SubTo(ptOut.Value, pt0.Value, pt1.Value)
}

// Neg performs -pt.
func (o *PlainOperator) Neg(pt *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := pt.ModLen()
	var auxLen int
	if pt.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.NegTo(ptOut, pt)

	return ptOut
}

// NegTo performs -pt and stores the result in ptOut.
func (o *PlainOperator) NegTo(ptOut *PlainPoly, ptIn *PlainPoly) {
	if !ptOut.IsConsistent(ptIn) {
		panic("inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(ptIn.HasAux, ptIn.ModLen())
	eval.NegTo(ptOut.Value, ptIn.Value)
}

// Mul returns pt0 * pt1.
func (o *PlainOperator) Mul(pt0, pt1 *PlainPoly) *PlainPoly {
	rank := o.Params.ringParams.Rank()
	modLen := pt0.ModLen()
	var auxLen int
	if pt0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	ptOut := NewPlainPolyCustom(rank, modLen, auxLen, false)
	o.MulTo(ptOut, pt0, pt1)

	return ptOut
}

// MulTo performs pt0 * pt1 and stores the result in ptOut.
func (o *PlainOperator) MulTo(ptOut *PlainPoly, pt0, pt1 *PlainPoly) {
	if !(ptOut.IsConsistent(pt0) && ptOut.IsConsistent(pt1)) {
		panic("MulTo: inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(pt0.HasAux, pt0.ModLen())
	eval.MulTo(ptOut.Value, pt0.Value, pt1.Value)
}

// MulAddTo performs ptOut += pt0 * pt1.
func (o *PlainOperator) MulAddTo(ptOut *PlainPoly, pt0, pt1 *PlainPoly) {
	if !(ptOut.IsConsistent(pt0) && ptOut.IsConsistent(pt1)) {
		panic("MulAddTo: inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(pt0.HasAux, pt0.ModLen())
	eval.MulAddTo(ptOut.Value, pt0.Value, pt1.Value)
}

// MulSubTo performs ptOut -= pt0 * pt1.
func (o *PlainOperator) MulSubTo(ptOut *PlainPoly, pt0, pt1 *PlainPoly) {
	if !(ptOut.IsConsistent(pt0) && ptOut.IsConsistent(pt1)) {
		panic("MulSubTo: inconsistent plaintexts")
	}

	eval := o.SubEvaluatorAt(pt0.HasAux, pt0.ModLen())
	eval.MulSubTo(ptOut.Value, pt0.Value, pt1.Value)
}

////////////////////////////////////////
///////// Ciphertext operations ////////
////////////////////////////////////////

// operatorBuffer is a buffer for the [Operator].
type operatorBuffer struct {
	// pNTT is the buffer for the NTT of the plaintext.
	pNTT *crt.Poly
	// pDiv is the buffer for division by auxiliary modulus.
	pDiv *crt.Poly
	// pScale is the buffer for ciphertext scaling.
	pScale *crt.Poly
	// pKsw is the buffer for key switching.
	pKsw *crt.Poly

	// ctGad is the buffer for gadget product.
	ctGad *Ciphertext
	// ctExt is the buffer for external product.
	ctExt *Ciphertext

	// decmp is the buffer for the gadget decomposition.
	decmp *Tensor
}

// newOperatorBuffer creates a new [operatorBuffer] for the given parameters.
func newOperatorBuffer(p Parameters) operatorBuffer {
	pNTT := crt.NewPoly(p.ringParams.Rank(), len(p.modulus)+len(p.auxModulus))
	var pDiv *crt.Poly
	if p.auxModulus != nil {
		pDiv = crt.NewPoly(p.ringParams.Rank(), len(p.modulus)+len(p.auxModulus))
	}
	pScale := crt.NewPoly(p.ringParams.Rank(), len(p.modulus))
	pKsw := crt.NewPoly(p.ringParams.Rank(), len(p.modulus))

	ctProd := NewCiphertext(p, true, false)
	ctExt := NewCiphertext(p, true, false)

	decmp := NewTensor(p, true, p.gadgetParams.gadgetLen(p), false)

	return operatorBuffer{
		pNTT:   pNTT,
		pDiv:   pDiv,
		pScale: pScale,
		pKsw:   pKsw,

		ctGad: ctProd,
		ctExt: ctExt,

		decmp: decmp,
	}
}

// Operator is a struct that performs RLWE operations.
type Operator struct {
	// Params is the RLWE parameters.
	Params Parameters

	// PlainOp is the plaintext operator.
	PlainOp *PlainOperator

	// Eval is the polynomial Eval.
	Eval crt.PolyEvaluator
	// Decmp is the decompose.
	Decmp Decomposer

	// buf is the buffer for the operator.
	buf operatorBuffer
}

// NewOperator creates a new [Operator] for the given parameters.
func NewOperator(p Parameters) *Operator {
	plainOp := NewPlainOperator(p)
	eval := plainOp.Eval
	decomposer := NewDecomposer(p)

	return &Operator{
		Params: p,

		PlainOp: plainOp,

		Eval:  eval,
		Decmp: decomposer,

		buf: newOperatorBuffer(p),
	}
}

// SafeCopy returns a thread-safe copy.
func (o *Operator) SafeCopy() *Operator {
	plainOp := o.PlainOp.SafeCopy()
	eval := plainOp.Eval
	decomposer := o.Decmp.SafeCopy()

	return &Operator{
		Params: o.Params,

		PlainOp: plainOp,

		Eval:  eval,
		Decmp: decomposer,

		buf: newOperatorBuffer(o.Params),
	}
}

// SubEvaluatorAt returns a sub evaluator at the given indices.
func (o *Operator) SubEvaluatorAt(hasAux bool, modLen int) crt.PolyEvaluator {
	return o.PlainOp.SubEvaluatorAt(hasAux, modLen)
}

// FwdNTT performs FwdNTT(c).
func (o *Operator) FwdNTT(c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c.ModLen()
	var auxLen int
	if c.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c.HasAux)
	o.FwdNTTTo(cOut, c)
	return cOut
}

// FwdNTTTo performs FwdNTT(c) and stores the result in cOut.
func (o *Operator) FwdNTTTo(cOut *Ciphertext, c *Ciphertext) {
	if !cOut.IsConsistent(c) {
		panic("inconsistent ciphertexts")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.FwdNTTTo(cOut.Body, c.Body)
	eval.FwdNTTTo(cOut.Mask, c.Mask)
}

// InvNTT performs InvNTT(c).
func (o *Operator) InvNTT(c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c.ModLen()
	var auxLen int
	if c.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c.HasAux)
	o.InvNTTTo(cOut, c)
	return cOut
}

// InvNTTTo performs InvNTT(c) and stores the result in cOut.
func (o *Operator) InvNTTTo(cOut *Ciphertext, c *Ciphertext) {
	if !cOut.IsConsistent(c) {
		panic("inconsistent ciphertexts")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.InvNTTTo(cOut.Body, c.Body)
	eval.InvNTTTo(cOut.Mask, c.Mask)
}

// Add performs c0 + c1.
func (o *Operator) Add(c0, c1 *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c0.ModLen()
	var auxLen int
	if c0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c0.HasAux)
	o.AddTo(cOut, c0, c1)
	return cOut
}

// AddTo performs c0 + c1 and stores the result in cOut.
func (o *Operator) AddTo(cOut *Ciphertext, c0, c1 *Ciphertext) {
	if !cOut.IsConsistent(c0) || !cOut.IsConsistent(c1) {
		panic("inconsistent ciphertexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	eval.AddTo(cOut.Body, c0.Body, c1.Body)
	eval.AddTo(cOut.Mask, c0.Mask, c1.Mask)
}

// ScalarAdd performs s + c.
func (o *Operator) ScalarAdd(s *PlainScalar, c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c.ModLen()
	var auxLen int
	if c.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c.HasAux)
	o.ScalarAddTo(cOut, s, c)
	return cOut
}

// ScalarAddTo performs c + s and stores the result in cOut.
func (o *Operator) ScalarAddTo(cOut *Ciphertext, s *PlainScalar, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == s.ModLen())) {
		panic("inconsistent ciphertexts or scalars")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.ScalarAddTo(cOut.Body, c.Body, s.Value)
	eval.ScalarAddTo(cOut.Mask, c.Mask, s.Value)
}

// PolyAdd performs p + c.
func (o *Operator) PolyAdd(p *PlainPoly, c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := p.ModLen()
	var auxLen int
	if p.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, p.HasAux)
	o.PolyAddTo(cOut, p, c)
	return cOut
}

// PolyAddTo performs p + c and stores the result in cOut.
func (o *Operator) PolyAddTo(cOut *Ciphertext, p *PlainPoly, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == p.ModLen())) {
		panic("inconsistent ciphertexts or polynomials")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.AddTo(cOut.Body, c.Body, p.Value)
	eval.AddTo(cOut.Mask, c.Mask, p.Value)
}

// Sub performs c0 - c1.
func (o *Operator) Sub(c0, c1 *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c0.ModLen()
	var auxLen int
	if c0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c0.HasAux)
	o.SubTo(cOut, c0, c1)
	return cOut
}

// SubTo performs c0 - c1 and stores the result in cOut.
func (o *Operator) SubTo(cOut *Ciphertext, c0, c1 *Ciphertext) {
	if !cOut.IsConsistent(c0) || !cOut.IsConsistent(c1) {
		panic("inconsistent ciphertexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	eval.SubTo(cOut.Body, c0.Body, c1.Body)
	eval.SubTo(cOut.Mask, c0.Mask, c1.Mask)
}

// ScalarSub performs s - c.
func (o *Operator) ScalarSub(s *PlainScalar, c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := s.ModLen()
	var auxLen int
	if s.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, s.HasAux)
	o.ScalarSubTo(cOut, s, c)
	return cOut
}

// ScalarSubTo performs s - c and stores the result in cOut.
func (o *Operator) ScalarSubTo(cOut *Ciphertext, s *PlainScalar, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == s.ModLen())) {
		panic("ScalarSubTo: inconsistent ciphertexts or scalars")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.ScalarSubTo(cOut.Body, c.Body, s.Value)
	eval.ScalarSubTo(cOut.Mask, c.Mask, s.Value)
}

// SubPoly performs c0 - c1.
func (o *Operator) SubPoly(c0, c1 *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c0.ModLen()
	var auxLen int
	if c0.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c0.HasAux)
	o.SubPolyTo(cOut, c0, c1)
	return cOut
}

// SubPolyTo performs c0 - c1 and stores the result in cOut.
func (o *Operator) SubPolyTo(cOut *Ciphertext, c0, c1 *Ciphertext) {
	if !cOut.IsConsistent(c0) || !cOut.IsConsistent(c1) {
		panic("inconsistent ciphertexts")
	}

	eval := o.SubEvaluatorAt(c0.HasAux, c0.ModLen())
	eval.SubTo(cOut.Body, c0.Body, c1.Body)
	eval.SubTo(cOut.Mask, c0.Mask, c1.Mask)
}

// Neg performs -c.
func (o *Operator) Neg(c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := c.ModLen()
	var auxLen int
	if c.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, c.HasAux)
	o.NegTo(cOut, c)
	return cOut
}

// NegTo performs -c and stores the result in cOut.
func (o *Operator) NegTo(cOut *Ciphertext, c *Ciphertext) {
	if !cOut.IsConsistent(c) {
		panic("inconsistent ciphertexts")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.NegTo(cOut.Body, c.Body)
	eval.NegTo(cOut.Mask, c.Mask)
}

// ScalarMul performs s * c.
func (o *Operator) ScalarMul(s *PlainScalar, c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := s.ModLen()
	var auxLen int
	if s.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, s.HasAux)
	o.ScalarMulTo(cOut, s, c)
	return cOut
}

// ScalarMulTo performs s * c and stores the result in cOut.
func (o *Operator) ScalarMulTo(cOut *Ciphertext, s *PlainScalar, c *Ciphertext) {
	if !cOut.IsConsistent(c) || !(cOut.ModLen() == s.ModLen()) {
		panic("inconsistent ciphertexts or scalars")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.ScalarMulTo(cOut.Body, c.Body, s.Value)
	eval.ScalarMulTo(cOut.Mask, c.Mask, s.Value)
}

// MulPoly performs p * c.
func (o *Operator) PolyMul(p *PlainPoly, c *Ciphertext) *Ciphertext {
	rank := o.Params.ringParams.Rank()
	modLen := p.ModLen()
	var auxLen int
	if p.HasAux {
		auxLen = len(o.Params.auxModulus)
	}

	cOut := NewCiphertextCustom(rank, modLen, auxLen, p.HasAux)
	o.PolyMulTo(cOut, p, c)
	return cOut
}

// PolyMulTo performs p * c and stores the result in cOut.
func (o *Operator) PolyMulTo(cOut *Ciphertext, p *PlainPoly, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == p.ModLen())) {
		panic("inconsistent ciphertexts or polynomials")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.MulTo(cOut.Body, c.Body, p.Value)
	eval.MulTo(cOut.Mask, c.Mask, p.Value)
}

// ScalarMulAddTo performs cOut += s * c.
func (o *Operator) ScalarMulAddTo(cOut *Ciphertext, s *PlainScalar, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == s.ModLen())) {
		panic("inconsistent ciphertexts or scalars")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.ScalarMulAddTo(cOut.Body, c.Body, s.Value)
	eval.ScalarMulAddTo(cOut.Mask, c.Mask, s.Value)
}

// PolyMulAddTo performs cOut += p * c.
func (o *Operator) PolyMulAddTo(cOut *Ciphertext, p *PlainPoly, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == p.ModLen())) {
		panic("inconsistent ciphertexts or polynomials")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.MulAddTo(cOut.Body, c.Body, p.Value)
	eval.MulAddTo(cOut.Mask, c.Mask, p.Value)
}

// ScalarMulSubTo performs cOut -= s * c.
func (o *Operator) ScalarMulSubTo(cOut *Ciphertext, s *PlainScalar, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == s.ModLen())) {
		panic("inconsistent ciphertexts or scalars")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.ScalarMulSubTo(cOut.Body, c.Body, s.Value)
	eval.ScalarMulSubTo(cOut.Mask, c.Mask, s.Value)
}

// PolyMulSubTo performs cOut -= p * c.
func (o *Operator) PolyMulSubTo(cOut *Ciphertext, p *PlainPoly, c *Ciphertext) {
	if !(cOut.IsConsistent(c) && (cOut.ModLen() == p.ModLen())) {
		panic("inconsistent ciphertexts or polynomials")
	}

	eval := o.SubEvaluatorAt(c.HasAux, c.ModLen())
	eval.MulSubTo(cOut.Body, c.Body, p.Value)
	eval.MulSubTo(cOut.Mask, c.Mask, p.Value)
}

// DivByAux performs round(cIn / auxModulus).
func (o *Operator) DivByAux(cIn *Ciphertext, isNTT bool) *Ciphertext {
	if !cIn.HasAux {
		panic("input ciphertext must have auxiliary modulus")
	}

	rank := o.Params.ringParams.Rank()
	modLen := cIn.ModLen() - len(o.Params.auxModulus)
	cOut := NewCiphertextCustom(rank, modLen, 0, false)
	o.DivByAuxTo(cOut, cIn, isNTT)
	return cOut
}

// DivByAuxTo performs cOut = round(cIn / auxModulus) and stores the result in cOut.
func (o *Operator) DivByAuxTo(cOut, cIn *Ciphertext, isNTT bool) {
	if o.Params.auxModulus == nil {
		panic("Auxiliary modulus is not set.")
	} else if !cIn.HasAux {
		panic("Input ciphertext must have auxiliary modulus.")
	} else if cOut.ModLen() != cIn.ModLen()-len(o.Params.auxModulus) {
		panic("Inconsistent ciphertext output.")
	} else if !(cOut.Body.Rank() == cIn.Body.Rank() && cOut.Mask.Rank() == cIn.Mask.Rank()) {
		panic("Inconsistent rank.")
	}

	auxLen := len(o.Params.auxModulus)
	modLen := cIn.ModLen() - auxLen
	evalMod := o.SubEvaluatorAt(false, modLen)
	evalAux := o.SubEvaluatorAt(true, auxLen)

	buf := o.buf.pDiv.WithModIdx(vec.Range(0, modLen+auxLen)...)
	bufMod := buf.WithModIdx(vec.Range(auxLen, modLen+auxLen)...)
	bufAux := buf.WithModIdx(vec.Range(0, auxLen)...)

	mod := evalMod.Modulus()
	auxMod := evalAux.Modulus()
	ss := crt.NewScaler(mod, auxMod)

	// Compute auxInvMod and modInvAux
	auxInvMod := crt.NewScalar(1, mod)
	modInvAux := crt.NewScalar(1, auxMod)
	for i, modi := range mod {
		for j, auxj := range auxMod {
			auxInvMod[i] = num.Mul(auxInvMod[i], num.Inv(auxj.Value(), modi), modi)
			modInvAux[j] = num.Mul(modInvAux[j], num.Inv(modi.Value(), auxj), auxj)
		}
	}

	// scale the body.
	buf.CopyFrom(cIn.Body)
	bufMod.IsNTT = buf.IsNTT
	bufAux.IsNTT = buf.IsNTT

	evalAux.ScalarMulTo(bufAux, bufAux, modInvAux)
	evalMod.ScalarMulTo(bufMod, bufMod, auxInvMod)

	if bufAux.IsNTT {
		evalAux.InvNTTTo(bufAux, bufAux)
	}

	ss.ScaleTo(cOut.Body, bufAux)
	cOut.Body.IsNTT = false

	if isNTT {
		evalMod.FwdNTTTo(cOut.Body, cOut.Body)
		if !bufMod.IsNTT {
			evalMod.FwdNTTTo(bufMod, bufMod)
		}
	} else {
		if bufMod.IsNTT {
			evalMod.InvNTTTo(bufMod, bufMod)
		}
	}
	evalMod.AddTo(cOut.Body, cOut.Body, bufMod)

	// scale the mask.
	buf.CopyFrom(cIn.Mask)
	bufMod.IsNTT = buf.IsNTT
	bufAux.IsNTT = buf.IsNTT

	evalAux.ScalarMulTo(bufAux, bufAux, modInvAux)
	evalMod.ScalarMulTo(bufMod, bufMod, auxInvMod)

	if bufAux.IsNTT {
		evalAux.InvNTTTo(bufAux, bufAux)
	}

	ss.ScaleTo(cOut.Mask, bufAux)
	cOut.Mask.IsNTT = false

	if isNTT {
		evalMod.FwdNTTTo(cOut.Mask, cOut.Mask)
		if !bufMod.IsNTT {
			evalMod.FwdNTTTo(bufMod, bufMod)
		}
	} else {
		if bufMod.IsNTT {
			evalMod.InvNTTTo(bufMod, bufMod)
		}
	}
	evalMod.AddTo(cOut.Mask, cOut.Mask, bufMod)

	// Set the auxiliary flag to false.
	cOut.HasAux = false
}

// Scale performs round(cIn / modulus).
func (o *Operator) Scale(cIn *Ciphertext, newLen int, isNTT bool) *Ciphertext {
	if newLen < 1 {
		panic("Invalid output modulus length.")
	}

	rank := o.Params.ringParams.Rank()
	cOut := NewCiphertextCustom(rank, newLen, 0, false)
	o.ScaleTo(cOut, cIn, newLen, isNTT)
	return cOut
}

// ScaleTo scales cIn to cOut and stores the result in cOut.
func (o *Operator) ScaleTo(cOut *Ciphertext, cIn *Ciphertext, newLen int, isNTT bool) {
	if cIn.HasAux {
		panic("Input ciphertext must not have auxiliary modulus.")
	} else if cIn.ModLen() > len(o.Params.modulus) {
		panic("Invalid input modulus length.")
	} else if cOut.ModLen() != newLen {
		panic("Inconsistent output modulus length.")
	} else if !(cOut.Body.Rank() == cIn.Body.Rank() && cOut.Mask.Rank() == cIn.Mask.Rank()) {
		panic("Inconsistent rank.")
	}

	oldLen := cIn.ModLen()

	if oldLen > newLen {
		evalOld := o.SubEvaluatorAt(false, oldLen)
		evalNew := evalOld.SubEvaluator(vec.Range(0, newLen)...)
		evalSc := evalOld.SubEvaluator(vec.Range(newLen, oldLen)...)

		bufOld := o.buf.pScale.WithModIdx(vec.Range(0, oldLen)...)
		bufNew := bufOld.WithModIdx(vec.Range(0, newLen)...)
		bufSc := bufOld.WithModIdx(vec.Range(newLen, oldLen)...)

		modOld := evalOld.Modulus()
		modNew := modOld[:newLen]
		modSc := modOld[newLen:]
		ss := crt.NewScaler(modNew, modSc)

		// Compute scInvNew and newInvSc
		scInvNew := crt.NewScalar(1, modNew)
		newInvSc := crt.NewScalar(1, modSc)
		for i, newi := range modNew {
			for j, scj := range modSc {
				scInvNew[i] = num.Mul(scInvNew[i], num.Inv(scj.Value(), newi), newi)
				newInvSc[j] = num.Mul(newInvSc[j], num.Inv(newi.Value(), scj), scj)
			}
		}

		// scale the body.
		bufOld.CopyFrom(cIn.Body)
		bufNew.IsNTT = bufOld.IsNTT
		bufSc.IsNTT = bufOld.IsNTT

		evalNew.ScalarMulTo(bufNew, bufNew, scInvNew)
		evalSc.ScalarMulTo(bufSc, bufSc, newInvSc)

		if bufSc.IsNTT {
			evalSc.InvNTTTo(bufSc, bufSc)
		}

		ss.ScaleTo(cOut.Body, bufSc)
		cOut.Body.IsNTT = false

		if isNTT {
			evalNew.FwdNTTTo(cOut.Body, cOut.Body)
			if !bufNew.IsNTT {
				evalNew.FwdNTTTo(bufNew, bufNew)
			}
		} else {
			if bufNew.IsNTT {
				evalNew.InvNTTTo(bufNew, bufNew)
			}
		}
		evalNew.AddTo(cOut.Body, cOut.Body, bufNew)

		// scale the mask.
		bufOld.CopyFrom(cIn.Mask)
		bufNew.IsNTT = bufOld.IsNTT
		bufSc.IsNTT = bufOld.IsNTT

		evalNew.ScalarMulTo(bufNew, bufNew, scInvNew)
		evalSc.ScalarMulTo(bufSc, bufSc, newInvSc)

		if bufSc.IsNTT {
			evalSc.InvNTTTo(bufSc, bufSc)
		}

		ss.ScaleTo(cOut.Mask, bufSc)
		cOut.Mask.IsNTT = false

		if isNTT {
			evalNew.FwdNTTTo(cOut.Mask, cOut.Mask)
			if !bufNew.IsNTT {
				evalNew.FwdNTTTo(bufNew, bufNew)
			}
		} else {
			if bufNew.IsNTT {
				evalNew.InvNTTTo(bufNew, bufNew)
			}
		}
		evalNew.AddTo(cOut.Mask, cOut.Mask, bufNew)
	} else if oldLen < newLen {
		evalNew := o.SubEvaluatorAt(false, newLen)
		modOld := evalNew.Modulus()[:oldLen]
		modSc := evalNew.Modulus()[oldLen:]

		sc := crt.NewScalar(1, modOld)
		for i, modi := range modOld {
			for _, modj := range modSc {
				sc[i] = num.Mul(sc[i], modj.Value(), modi)
			}
		}

		cOut.Body.IsNTT = cIn.Body.IsNTT
		cOut.Mask.IsNTT = cIn.Mask.IsNTT

		for i := 0; i < oldLen; i++ {
			vec.ScalarMulTo(cOut.Body.Coeffs[i], cIn.Body.Coeffs[i], sc[i], evalNew.Modulus()[i])
			vec.ScalarMulTo(cOut.Mask.Coeffs[i], cIn.Mask.Coeffs[i], sc[i], evalNew.Modulus()[i])
		}
		for i := oldLen; i < newLen; i++ {
			clear(cOut.Body.Coeffs[i])
			clear(cOut.Mask.Coeffs[i])
		}

		// perform NTT if needed.
		if isNTT && !cOut.Body.IsNTT {
			evalNew.FwdNTTTo(cOut.Body, cOut.Body)
			evalNew.FwdNTTTo(cOut.Mask, cOut.Mask)
		} else if !isNTT && cOut.Body.IsNTT {
			evalNew.InvNTTTo(cOut.Body, cOut.Body)
			evalNew.InvNTTTo(cOut.Mask, cOut.Mask)
		}
	} else {
		eval := o.SubEvaluatorAt(false, newLen)
		cOut.CopyFrom(cIn)

		// perform NTT if needed.
		if isNTT && !cOut.Body.IsNTT {
			eval.FwdNTTTo(cOut.Body, cOut.Body)
			eval.FwdNTTTo(cOut.Mask, cOut.Mask)
		} else if !isNTT && cOut.Body.IsNTT {
			eval.InvNTTTo(cOut.Body, cOut.Body)
			eval.InvNTTTo(cOut.Mask, cOut.Mask)
		}
	}
}

// HoistedGadgetProdLazy performs a lazy hoisted gadget product and returns the result.
func (o *Operator) HoistedGadgetProdLazy(decmp *Tensor, gadenc *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := decmp.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.HoistedGadgetProdLazyTo(cOut, decmp, gadenc, isNTT)
	return cOut
}

// HoistedGadgetProdLazyTo performs a lazy hoisted gadget product and stores the result in cOut.
func (o *Operator) HoistedGadgetProdLazyTo(cOut *Ciphertext, decmp *Tensor, gadenc *GadgetEncryption, isNTT bool) {
	if !(decmp.HasAux == gadenc.Value[0].HasAux && decmp.HasAux == (o.Params.auxModulus != nil)) {
		panic("Inconsistent auxiliary flag.")
	} else if cOut.HasAux != decmp.HasAux {
		panic("Inconsistent auxiliary flag.")
	}

	gadLen := decmp.Degree()
	modLen := decmp.ModLen()
	eval := o.SubEvaluatorAt(decmp.HasAux, modLen)

	if cOut.ModLen() != modLen {
		panic("Inconsistent ciphertext output.")
	} else if gadLen != o.Decmp.DecomposeLen(modLen) {
		panic("Inconsistent decomposition length.")
	}

	cOut.Clear()
	cOut.Body.IsNTT = true
	cOut.Mask.IsNTT = true
	cOut.HasAux = false

	for i := 0; i < gadLen; i++ {
		if !decmp.Value[i].IsNTT {
			panic("Decomposition must be in NTT form.")
		}

		gadenci := gadenc.Value[i].WithModIdx(vec.Range(0, modLen)...)
		eval.MulAddTo(cOut.Body, decmp.Value[i], gadenci.Body)
		eval.MulAddTo(cOut.Mask, decmp.Value[i], gadenci.Mask)
	}

	if !isNTT {
		eval.InvNTTTo(cOut.Body, cOut.Body)
		eval.InvNTTTo(cOut.Mask, cOut.Mask)
	}
}

// HoistedGadgetProd performs a hoisted gadget product and returns the result.
func (o *Operator) HoistedGadgetProd(decmp *Tensor, gadenc *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := decmp.ModLen()
	if o.Params.auxModulus != nil {
		modLen -= len(o.Params.auxModulus)
	}
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.HoistedGadgetProdTo(cOut, decmp, gadenc, isNTT)
	return cOut
}

// HoistedGadgetProdTo performs a hoisted gadget product and stores the result in cOut.
func (o *Operator) HoistedGadgetProdTo(cOut *Ciphertext, decmp *Tensor, gadenc *GadgetEncryption, isNTT bool) {
	modLen := decmp.ModLen()
	buf := o.buf.ctGad.WithModIdx(vec.Range(0, modLen)...)
	o.HoistedGadgetProdLazyTo(buf, decmp, gadenc, true)

	if decmp.HasAux {
		cOut.CopyFrom(buf)
	} else {
		o.DivByAuxTo(cOut, buf, isNTT)
	}

	eval := o.SubEvaluatorAt(false, cOut.ModLen())
	if isNTT && !cOut.Body.IsNTT {
		eval.FwdNTTTo(cOut.Body, cOut.Body)
		eval.FwdNTTTo(cOut.Mask, cOut.Mask)
	} else if !isNTT && cOut.Body.IsNTT {
		eval.InvNTTTo(cOut.Body, cOut.Body)
		eval.InvNTTTo(cOut.Mask, cOut.Mask)
	}
}

// GadgetProdLazy performs a lazy gadget product and returns the result.
func (o *Operator) GadgetProdLazy(pIn *PlainPoly, gadenc *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := pIn.ModLen()
	if o.Params.auxModulus != nil {
		modLen += len(o.Params.auxModulus)
	}
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.GadgetProdLazyTo(cOut, pIn, gadenc, isNTT)
	return cOut
}

// GadgetProdLazyTo performs a lazy gadget product and stores the result in cOut.
func (o *Operator) GadgetProdLazyTo(cOut *Ciphertext, pIn *PlainPoly, gadenc *GadgetEncryption, isNTT bool) {
	if pIn.HasAux {
		panic("Plain polynomial must not have auxiliary modulus.")
	}

	modLen := pIn.ModLen()
	decmpLen := o.Decmp.DecomposeLen(modLen)
	decmp := o.buf.decmp.WithDegreeAndModIdx(decmpLen, vec.Range(0, modLen)...)
	o.Decmp.DecomposeTo(decmp, pIn.Value)

	o.HoistedGadgetProdLazyTo(cOut, decmp, gadenc, isNTT)
}

// GadgetProd performs a gadget product and returns the result.
func (o *Operator) GadgetProd(pIn *PlainPoly, gadenc *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := pIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.GadgetProdTo(cOut, pIn, gadenc, isNTT)
	return cOut
}

// GadgetProdTo performs a gadget product and stores the result in cOut.
func (o *Operator) GadgetProdTo(cOut *Ciphertext, pIn *PlainPoly, gadenc *GadgetEncryption, isNTT bool) {
	if pIn.HasAux {
		panic("Plain polynomial must not have auxiliary modulus.")
	}

	modLen := pIn.ModLen()
	decmpLen := o.Decmp.DecomposeLen(modLen)
	decmp := o.buf.decmp.WithDegreeAndModIdx(decmpLen, vec.Range(0, modLen)...)
	o.Decmp.DecomposeTo(decmp, pIn.Value)

	o.HoistedGadgetProdTo(cOut, decmp, gadenc, isNTT)
}

// Relin performs a relinearisation and returns the result.
func (o *Operator) Relin(cIn *Tensor, rlk *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := cIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.RelinTo(cOut, cIn, rlk, isNTT)
	return cOut
}

// RelinTo performs a relinearisation and stores the result in cOut.
func (o *Operator) RelinTo(cOut *Ciphertext, cIn *Tensor, rlk *GadgetEncryption, isNTT bool) {
	if cIn.HasAux {
		panic("Ciphertext must not have auxiliary modulus.")
	} else if cIn.Degree() != 3 {
		panic("Ciphertext must have degree 3.")
	}

	modLen := cIn.ModLen()
	eval := o.SubEvaluatorAt(false, modLen)
	buf := o.buf.pNTT.WithModIdx(vec.Range(0, modLen)...)

	o.GadgetProdTo(cOut, &PlainPoly{Value: cIn.Value[2]}, rlk, isNTT)

	if cIn.Value[0].IsNTT && !isNTT {
		eval.InvNTTTo(buf, cIn.Value[0])
	} else if !cIn.Value[0].IsNTT && isNTT {
		eval.FwdNTTTo(buf, cIn.Value[0])
	} else {
		buf.CopyFrom(cIn.Value[0])
	}

	eval.AddTo(cOut.Body, cOut.Body, buf)

	if cIn.Value[1].IsNTT && !isNTT {
		eval.InvNTTTo(buf, cIn.Value[1])
	} else if !cIn.Value[1].IsNTT && isNTT {
		eval.FwdNTTTo(buf, cIn.Value[1])
	} else {
		buf.CopyFrom(cIn.Value[1])
	}

	eval.AddTo(cOut.Mask, cOut.Mask, buf)
}

// KeySwitch performs a key switch and returns the result.
func (o *Operator) KeySwitch(cIn *Ciphertext, ksk *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := cIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.KeySwitchTo(cOut, cIn, ksk, isNTT)
	return cOut
}

// KeySwitchTo performs a key switch and stores the result in cOut.
func (o *Operator) KeySwitchTo(cOut *Ciphertext, cIn *Ciphertext, ksk *GadgetEncryption, isNTT bool) {
	if cIn.HasAux {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	modLen := cIn.ModLen()
	decmpLen := o.Decmp.DecomposeLen(modLen)
	decmp := o.buf.decmp.WithDegreeAndModIdx(decmpLen, vec.Range(0, modLen)...)
	o.Decmp.DecomposeTo(decmp, cIn.Mask)

	o.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)
}

// HoistedKeySwitch performs a hoisted key switch and returns the result.
func (o *Operator) HoistedKeySwitch(decmp *Tensor, cIn *Ciphertext, ksk *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := cIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)
	return cOut
}

// HoistedKeySwitchTo performs a hoisted key switch and stores the result in cOut.
func (o *Operator) HoistedKeySwitchTo(cOut *Ciphertext, decmp *Tensor, cIn *Ciphertext, ksk *GadgetEncryption, isNTT bool) {
	if cIn.HasAux {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	modLen := cIn.ModLen()
	eval := o.SubEvaluatorAt(false, modLen)
	buf := o.buf.pKsw.WithModIdx(vec.Range(0, modLen)...)

	if cIn.Body.IsNTT && !isNTT {
		eval.InvNTTTo(buf, cIn.Body)
	} else if !cIn.Body.IsNTT && isNTT {
		eval.FwdNTTTo(buf, cIn.Body)
	} else {
		buf.CopyFrom(cIn.Body)
	}

	o.HoistedGadgetProdTo(cOut, decmp, ksk, isNTT)
	eval.AddTo(cOut.Body, cOut.Body, buf)
}

// Aut performs an automorphism and returns the result.
func (o *Operator) Aut(idx int, cIn *Ciphertext, atk *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := cIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.AutTo(cOut, idx, cIn, atk, isNTT)
	return cOut
}

// AutTo performs an automorphism and stores the result in cOut.
func (o *Operator) AutTo(cOut *Ciphertext, idx int, cIn *Ciphertext, atk *GadgetEncryption, isNTT bool) {
	modLen := cIn.ModLen()
	eval := o.SubEvaluatorAt(false, modLen)

	o.KeySwitchTo(cOut, cIn, atk, isNTT)
	eval.AutTo(cOut.Body, cOut.Body, idx)
	eval.AutTo(cOut.Mask, cOut.Mask, idx)
}

// HoistedAut performs a hoisted automorphism and returns the result.
func (o *Operator) HoistedAut(idx int, decmp *Tensor, cIn *Ciphertext, atk *GadgetEncryption, isNTT bool) *Ciphertext {
	modLen := cIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.HoistedAutTo(cOut, idx, decmp, cIn, atk, isNTT)
	return cOut
}

// HoistedAutTo performs a hoisted automorphism and stores the result in cOut.
func (o *Operator) HoistedAutTo(cOut *Ciphertext, idx int, decmp *Tensor, cIn *Ciphertext, atk *GadgetEncryption, isNTT bool) {
	modLen := cIn.ModLen()
	eval := o.SubEvaluatorAt(false, modLen)

	o.HoistedKeySwitchTo(cOut, decmp, cIn, atk, isNTT)
	eval.AutTo(cOut.Body, cOut.Body, idx)
	eval.AutTo(cOut.Mask, cOut.Mask, idx)
}

// ExtProd performs an external product and returns the result.
func (o *Operator) ExtProd(cIn *Ciphertext, gsw *RGSW, isNTT bool) *Ciphertext {
	modLen := cIn.ModLen()
	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
	o.ExtProdTo(cOut, cIn, gsw, isNTT)
	return cOut
}

// ExtProdTo performs an external product and stores the result in cOut.
func (o *Operator) ExtProdTo(cOut *Ciphertext, cIn *Ciphertext, gsw *RGSW, isNTT bool) {
	if cIn.HasAux {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	modLen := cIn.ModLen()
	buf := o.buf.ctExt.WithModIdx(vec.Range(0, modLen)...)

	o.GadgetProdTo(buf, &PlainPoly{Value: cIn.Body}, gsw.Body, isNTT)
	o.GadgetProdTo(cOut, &PlainPoly{Value: cIn.Mask}, gsw.Mask, isNTT)

	o.AddTo(cOut, cOut, buf)
}

// HoistedExtProdTo performs a hoisted external product and stores the result in cOut.
func (o *Operator) HoistedExtProdTo(cOut *Ciphertext, decmpBody *Tensor, decmpMask *Tensor, gsw *RGSW, isNTT bool) {
	modLen := decmpBody.ModLen()
	if o.Params.auxModulus != nil {
		modLen -= len(o.Params.auxModulus)
	}
	buf := o.buf.ctExt.WithModIdx(vec.Range(0, modLen)...)

	o.HoistedGadgetProdLazyTo(buf, decmpBody, gsw.Body, isNTT)
	o.HoistedGadgetProdLazyTo(cOut, decmpMask, gsw.Mask, isNTT)

	o.AddTo(cOut, cOut, buf)
}
