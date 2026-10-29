package rlwe

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// Operator evaluates homomorphic operations over [Ciphertext] values.
type Operator struct {
	params Parameters

	plainOp *PlainOperator
	dcmp    *Decomposer

	poolPoly       *pool.Pool[*Element]
	poolRLWE       *pool.Pool[*Ciphertext]
	poolPolyDcmp   *pool.Pool[*Vector]
	poolScalarDcmp *pool.Pool[*Vector]
}

// NewOperator creates a new [Operator].
func NewOperator(params Parameters) *Operator {
	return &Operator{
		params: params,

		plainOp: NewPlainOperator(params),
		dcmp:    NewDecomposer(params),

		poolPoly: pool.NewPool(func() *Element {
			return NewPoly(params, params.AuxModLen() > 0, dft.FormCoeff)
		}),
		poolRLWE: pool.NewPool(func() *Ciphertext {
			return NewCiphertext(params, params.AuxModLen() > 0, dft.FormCoeff)
		}),
		poolPolyDcmp: pool.NewPool(func() *Vector {
			return NewVector(params, params.GadgetLen(), params.AuxModLen() > 0, crt.TypePoly, dft.FormCoeff)
		}),
		poolScalarDcmp: pool.NewPool(func() *Vector {
			return NewVector(params, params.GadgetLen(), params.AuxModLen() > 0, crt.TypeScalar, dft.FormCoeff)
		}),
	}
}

// Params returns the [Parameters] of the operator.
func (op *Operator) Params() Parameters {
	return op.params
}

// PlainOp returns the [PlainOperator] of the operator.
func (op *Operator) PlainOp() *PlainOperator {
	return op.plainOp
}

// CRTOp returns the [crt.Operator] of the operator.
func (op *Operator) CRTOp() *crt.Operator {
	return op.params.crtOp
}

// Decomposer returns the [Decomposer] of the operator.
func (op *Operator) Decomposer() *Decomposer {
	return op.dcmp
}

// NewCiphertext creates a new [Ciphertext].
func (op *Operator) NewCiphertext(withAux bool, form dft.Form) *Ciphertext {
	return NewCiphertext(op.params, withAux, form)
}

// FwdNTT returns FwdNTT(ct).
func (op *Operator) FwdNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormNTT)
	op.FwdNTTTo(ctOut, ct)
	return ctOut
}

// FwdNTTTo computes ctOut = NTT(ct).
func (op *Operator) FwdNTTTo(ctOut, ct *Ciphertext) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.FwdNTTTo(ctOut.Body, ct.Body)
	op.plainOp.FwdNTTTo(ctOut.Mask, ct.Mask)
}

// InvNTT returns InvNTT(ct).
func (op *Operator) InvNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.InvNTTTo(ctOut, ct)
	return ctOut
}

// InvNTTTo computes ctOut = InvNTT(ct).
func (op *Operator) InvNTTTo(ctOut, ct *Ciphertext) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.InvNTTTo(ctOut.Body, ct.Body)
	op.plainOp.InvNTTTo(ctOut.Mask, ct.Mask)
}

// Add returns ct + e.
func (op *Operator) Add[T *Element | *Ciphertext](ct *Ciphertext, e T) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.AddTo(ctOut, ct, e)
	return ctOut
}

// AddTo computes ctOut = ct + e.
func (op *Operator) AddTo[T *Element | *Ciphertext](ctOut, ct *Ciphertext, e T) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	switch e := any(e).(type) {
	case *Element:
		op.plainOp.AddTo(ctOut.Body, ct.Body, e)
		ctOut.Mask.CopyFrom(ct.Mask)
	case *Ciphertext:
		checkCiphertextShape(op.params, e)
		op.plainOp.AddTo(ctOut.Body, ct.Body, e.Body)
		op.plainOp.AddTo(ctOut.Mask, ct.Mask, e.Mask)
	}
}

// Sub returns ct - e.
func (op *Operator) Sub[T *Element | *Ciphertext](ct *Ciphertext, e T) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.SubTo(ctOut, ct, e)
	return ctOut
}

// SubTo computes ctOut = ct - e.
func (op *Operator) SubTo[T *Element | *Ciphertext](ctOut, ct *Ciphertext, e T) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	switch e := any(e).(type) {
	case *Element:
		op.plainOp.SubTo(ctOut.Body, ct.Body, e)
		ctOut.Mask.CopyFrom(ct.Mask)
	case *Ciphertext:
		checkCiphertextShape(op.params, e)
		op.plainOp.SubTo(ctOut.Body, ct.Body, e.Body)
		op.plainOp.SubTo(ctOut.Mask, ct.Mask, e.Mask)
	}
}

// Neg returns -ct.
func (op *Operator) Neg(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.NegTo(ctOut, ct)
	return ctOut
}

// NegTo computes ctOut = -ct.
func (op *Operator) NegTo(ctOut, ct *Ciphertext) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.NegTo(ctOut.Body, ct.Body)
	op.plainOp.NegTo(ctOut.Mask, ct.Mask)
}

// Mul returns ct * e.
//
// Assumes that inputs do not overlap.
func (op *Operator) Mul(ct *Ciphertext, e *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.MulTo(ctOut, ct, e)
	return ctOut
}

// MulTo computes ctOut = ct * e.
//
// Assumes that inputs do not overlap.
func (op *Operator) MulTo(ctOut, ct *Ciphertext, e *Element) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.MulTo(ctOut.Body, ct.Body, e)
	op.plainOp.MulTo(ctOut.Mask, ct.Mask, e)
}

// MulAddTo computes ctOut += ct * e.
//
// Assumes that inputs do not overlap.
func (op *Operator) MulAddTo(ctOut, ct *Ciphertext, e *Element) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.MulAddTo(ctOut.Body, ct.Body, e)
	op.plainOp.MulAddTo(ctOut.Mask, ct.Mask, e)
}

// MulSubTo computes ctOut -= ct * e.
//
// Assumes that inputs do not overlap.
func (op *Operator) MulSubTo(ctOut, ct *Ciphertext, e *Element) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.MulSubTo(ctOut.Body, ct.Body, e)
	op.plainOp.MulSubTo(ctOut.Mask, ct.Mask, e)
}

// Tensor returns ct0 * ct1.
//
// Assumes that inputs do not overlap.
func (op *Operator) Tensor(ct0, ct1 *Ciphertext) *Vector {
	ctVecOut := NewVectorCustom(ct0.Rank(), ct0.BaseModLen(), ct0.AuxModLen(), 3, crt.TypePoly, dft.FormCoeff)
	op.TensorTo(ctVecOut, ct0, ct1)
	return ctVecOut
}

// TensorTo computes ctVecOut = ct0 * ct1.
// Length of ctVecOut must be 3.
//
// Assumes that inputs do not overlap.
func (op *Operator) TensorTo(ctVecOut *Vector, ct0, ct1 *Ciphertext) {
	checkCiphertextShape(op.params, ctVecOut)
	checkCiphertextShape(op.params, ct0)
	checkCiphertextShape(op.params, ct1)
	if ctVecOut.Len() != 3 {
		panic("inconsistent inputs")
	}

	op.plainOp.MulTo(ctVecOut.Value[0], ct0.Body, ct1.Body)
	op.plainOp.MulTo(ctVecOut.Value[1], ct0.Body, ct1.Mask)
	op.plainOp.MulAddTo(ctVecOut.Value[1], ct0.Mask, ct1.Body)
	op.plainOp.MulTo(ctVecOut.Value[2], ct0.Mask, ct1.Mask)
}

// CanAut returns if automorphism X -> X^idx is valid.
func (op *Operator) CanAut(idx int) bool {
	return op.plainOp.CanAut(idx)
}

// Embed returns an embedding of ct to modulus lengths outBaseLen and outAuxLen.
func (op *Operator) Embed(ct *Ciphertext, outBaseLen, outAuxLen int, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), outBaseLen, outAuxLen, outForm)
	op.EmbedTo(ctOut, ct, outForm)
	return ctOut
}

// EmbedTo embeds ct to ctOut.
func (op *Operator) EmbedTo(ctOut, ct *Ciphertext, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.EmbedTo(ctOut.Body, ct.Body, outForm)
	op.plainOp.EmbedTo(ctOut.Mask, ct.Mask, outForm)
}

// Scale returns a scaling of ct to base modulus length outBaseLen without auxiliary modulus.
func (op *Operator) Scale(ct *Ciphertext, outBaseLen int, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), outBaseLen, 0, outForm)
	op.ScaleTo(ctOut, ct, outForm)
	return ctOut
}

// ScaleTo scales ct to ctOut.
func (op *Operator) ScaleTo(ctOut, ct *Ciphertext, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.plainOp.ScaleTo(ctOut.Body, ct.Body, outForm)
	op.plainOp.ScaleTo(ctOut.Mask, ct.Mask, outForm)
}

// GadgetProd returns p * ctGadEnc.
func (op *Operator) GadgetProd(p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), p.BaseModLen(), 0, dft.FormCoeff)
	op.GadgetProdTo(ctOut, p, ctGadEnc, outForm)
	return ctOut
}

// GadgetProdTo computes ctOut = p * ctGadEnc.
func (op *Operator) GadgetProdTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	if ctOut.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctOutAuxFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutAuxFull)

	ctOutAux := ctOutAuxFull.WithModLen(ctOut.BaseModLen(), op.params.AuxModLen())
	op.GadgetProdLazyTo(ctOutAux, p, ctGadEnc, outForm)

	if op.params.AuxModLen() > 0 {
		op.ScaleTo(ctOut, ctOutAux, outForm)
		return
	}

	ctOut.CopyFrom(ctOutAux)
}

// GadgetProdHoisted returns pDcmp * ctGadEnc.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadgetProdHoisted(pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), pDcmp.BaseModLen(), 0, dft.FormCoeff)
	op.GadgetProdHoistedTo(ctOut, pDcmp, ctGadEnc, outForm)
	return ctOut
}

// GadgetProdHoistedTo computes ctOut = pDcmp * ctGadEnc.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadgetProdHoistedTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	if ctOut.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctOutAuxFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutAuxFull)

	ctOutAux := ctOutAuxFull.WithModLen(ctOut.BaseModLen(), pDcmp.AuxModLen())
	op.GadgetProdHoistedLazyTo(ctOutAux, pDcmp, ctGadEnc, outForm)

	if op.params.AuxModLen() > 0 {
		op.ScaleTo(ctOut, ctOutAux, outForm)
		return
	}

	ctOut.CopyFrom(ctOutAux)
}

// GadgetProdLazy returns p * ctGadEnc,
// but the output modulus includes the auxiliary modulus.
func (op *Operator) GadgetProdLazy(p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), p.BaseModLen(), op.params.AuxModLen(), dft.FormCoeff)
	op.GadgetProdLazyTo(ctOut, p, ctGadEnc, outForm)
	return ctOut
}

// GadgetProdLazyTo computes ctOut = p * ctGadEnc,
// but the modulus of ctOut includes the auxiliary modulus.
func (op *Operator) GadgetProdLazyTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkElementShape(op.params, p)
	checkCiphertextShape(op.params, ctGadEnc)
	if p.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	var pDcmpFull *Vector
	switch {
	case p.Type() == crt.TypeScalar:
		pDcmpFull = op.poolScalarDcmp.Get()
		defer op.poolScalarDcmp.Put(pDcmpFull)
	case p.Type() == crt.TypePoly:
		pDcmpFull = op.poolPolyDcmp.Get()
		defer op.poolPolyDcmp.Put(pDcmpFull)
	}
	pDcmp := pDcmpFull.Slice(0, op.dcmp.DecomposeLen(p.BaseModLen())).WithModLen(p.BaseModLen(), op.params.AuxModLen())

	op.dcmp.DecomposeTo(pDcmp, p, dft.FormNTT)
	op.GadgetProdHoistedLazyTo(ctOut, pDcmp, ctGadEnc, outForm)
}

// GadgetProdHoistedLazy returns pDcmp * ctGadEnc,
// but the output modulus includes the auxiliary modulus.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadgetProdHoistedLazy(pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), pDcmp.BaseModLen(), pDcmp.AuxModLen(), dft.FormCoeff)
	op.GadgetProdHoistedLazyTo(ctOut, pDcmp, ctGadEnc, outForm)
	return ctOut
}

// GadgetProdHoistedLazyTo computes ctOut = pDcmp * ctGadEnc,
// but the modulus of ctOut includes the auxiliary modulus.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadgetProdHoistedLazyTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkElementShape(op.params, pDcmp)
	checkCiphertextShape(op.params, ctGadEnc)
	if ctOut.BaseModLen() != pDcmp.BaseModLen() || ctOut.AuxModLen() != op.params.AuxModLen() ||
		pDcmp.Len() != op.dcmp.DecomposeLen(pDcmp.BaseModLen()) {
		panic("inconsistent inputs")
	}

	ctOutBufFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutBufFull)

	ctOutBuf := ctOutBufFull.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen())
	op.plainOp.MulTo(ctOutBuf.Body, pDcmp.Value[0], ctGadEnc.Value[0].Body.WithModLen(pDcmp.BaseModLen(), pDcmp.AuxModLen()))
	op.plainOp.MulTo(ctOutBuf.Mask, pDcmp.Value[0], ctGadEnc.Value[0].Mask.WithModLen(pDcmp.BaseModLen(), pDcmp.AuxModLen()))
	for i := 1; i < pDcmp.Len(); i++ {
		op.plainOp.MulAddTo(ctOutBuf.Body, pDcmp.Value[i], ctGadEnc.Value[i].Body.WithModLen(pDcmp.BaseModLen(), pDcmp.AuxModLen()))
		op.plainOp.MulAddTo(ctOutBuf.Mask, pDcmp.Value[i], ctGadEnc.Value[i].Mask.WithModLen(pDcmp.BaseModLen(), pDcmp.AuxModLen()))
	}

	switch {
	case ctOutBuf.Form() == dft.FormNTT && outForm == dft.FormCoeff:
		op.InvNTTTo(ctOut, ctOutBuf)
	case ctOutBuf.Form() == dft.FormCoeff && outForm == dft.FormNTT:
		op.FwdNTTTo(ctOut, ctOutBuf)
	case ctOutBuf.Form() == outForm:
		ctOut.CopyFrom(ctOutBuf)
	}
}

// KeySwitch returns the key-switched ct using ksk.
func (op *Operator) KeySwitch(ct *Ciphertext, ksk *KeySwitchKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), 0, dft.FormCoeff)
	op.KeySwitchTo(ctOut, ct, ksk, outForm)
	return ctOut
}

// KeySwitchTo key-switches ct to ctOut using ksk.
func (op *Operator) KeySwitchTo(ctOut, ct *Ciphertext, ksk *KeySwitchKey, outForm dft.Form) {
	ctBodyFull := op.poolPoly.Get()
	defer op.poolPoly.Put(ctBodyFull)

	ctBody := ctBodyFull.WithModLen(ct.BaseModLen(), ct.AuxModLen())
	switch {
	case ct.Form() == dft.FormNTT && outForm == dft.FormCoeff:
		op.plainOp.InvNTTTo(ctBody, ct.Body)
	case ct.Form() == dft.FormCoeff && outForm == dft.FormNTT:
		op.plainOp.FwdNTTTo(ctBody, ct.Body)
	case ct.Form() == outForm:
		ctBody.CopyFrom(ct.Body)
	}

	op.GadgetProdTo(ctOut, ct.Mask, (*GadgetEncryption)(ksk), outForm)
	op.plainOp.AddTo(ctOut.Body, ctOut.Body, ctBody)
}

// KeySwitchHoisted returns the key-switched ct using ksk.
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) KeySwitchHoisted(ct *Ciphertext, ctMaskDcmp *Vector, ksk *KeySwitchKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), 0, dft.FormCoeff)
	op.KeySwitchHoistedTo(ctOut, ct, ctMaskDcmp, ksk, outForm)
	return ctOut
}

// KeySwitchHoistedTo key-switches ct to ctOut using ksk.
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) KeySwitchHoistedTo(ctOut, ct *Ciphertext, ctMaskDcmp *Vector, ksk *KeySwitchKey, outForm dft.Form) {
	ctBodyFull := op.poolPoly.Get()
	defer op.poolPoly.Put(ctBodyFull)

	ctBody := ctBodyFull.WithModLen(ct.BaseModLen(), ct.AuxModLen())
	switch {
	case ct.Form() == dft.FormNTT && outForm == dft.FormCoeff:
		op.plainOp.InvNTTTo(ctBody, ct.Body)
	case ct.Form() == dft.FormCoeff && outForm == dft.FormNTT:
		op.plainOp.FwdNTTTo(ctBody, ct.Body)
	case ct.Form() == outForm:
		ctBody.CopyFrom(ct.Body)
	}

	op.GadgetProdHoistedTo(ctOut, ctMaskDcmp, (*GadgetEncryption)(ksk), outForm)
	op.plainOp.AddTo(ctOut.Body, ctOut.Body, ctBody)
}

// KeySwitchLazy returns the key-switched ct using ksk,
// but the output modulus includes the auxiliary modulus.
func (op *Operator) KeySwitchLazy(ct *Ciphertext, ksk *KeySwitchKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), op.params.AuxModLen(), dft.FormCoeff)
	op.KeySwitchLazyTo(ctOut, ct, ksk, outForm)
	return ctOut
}

// KeySwitchLazyTo key-switches ct to ctOut using ksk,
// but the modulus of ctOut includes auxiliary modulus.
func (op *Operator) KeySwitchLazyTo(ctOut, ct *Ciphertext, ksk *KeySwitchKey, outForm dft.Form) {
	ctBodyFull := op.poolPoly.Get()
	defer op.poolPoly.Put(ctBodyFull)

	ctBody := ctBodyFull.WithModLen(ct.BaseModLen(), op.params.AuxModLen())

	opIn := op.plainOp.withModLen(ct.BaseModLen(), 0)
	opOut := op.plainOp.withModLen(ct.BaseModLen(), op.params.AuxModLen())
	sc := crt.NewScaler(opOut, opIn).WithPool(op.plainOp.poolEmb)
	sc.ScaleTo(ctBody.Value, ct.Body.Value, outForm)

	op.GadgetProdLazyTo(ctOut, ct.Mask, (*GadgetEncryption)(ksk), outForm)
	op.plainOp.AddTo(ctOut.Body, ctOut.Body, ctBody)
}

// KeySwitchHoistedLazy returns the key-switched ct using ksk,
// but the output modulus includes the auxiliary modulus.
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) KeySwitchHoistedLazy(ct *Ciphertext, ctMaskDcmp *Vector, ksk *KeySwitchKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ctMaskDcmp.AuxModLen(), dft.FormCoeff)
	op.KeySwitchHoistedLazyTo(ctOut, ct, ctMaskDcmp, ksk, outForm)
	return ctOut
}

// KeySwitchHoistedLazyTo key-switches ct to ctOut using ksk,
// but the modulus of ctOut includes auxiliary modulus.
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) KeySwitchHoistedLazyTo(ctOut, ct *Ciphertext, ctMaskDcmp *Vector, ksk *KeySwitchKey, outForm dft.Form) {
	ctBodyFull := op.poolPoly.Get()
	defer op.poolPoly.Put(ctBodyFull)

	ctBody := ctBodyFull.WithModLen(ct.BaseModLen(), ctMaskDcmp.AuxModLen())

	opIn := op.plainOp.withModLen(ct.BaseModLen(), 0)
	opOut := op.plainOp.withModLen(ct.BaseModLen(), ctMaskDcmp.AuxModLen())
	sc := crt.NewScaler(opOut, opIn).WithPool(op.plainOp.poolEmb)
	sc.ScaleTo(ctBody.Value, ct.Body.Value, outForm)

	op.GadgetProdHoistedLazyTo(ctOut, ctMaskDcmp, (*GadgetEncryption)(ksk), outForm)
	op.plainOp.AddTo(ctOut.Body, ctOut.Body, ctBody)
}

// Aut returns ct(X^idx).
func (op *Operator) Aut(ct *Ciphertext, atk *AutKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.AutTo(ctOut, ct, atk, outForm)
	return ctOut
}

// AutTo computes ctOut = ct(X^idx).
func (op *Operator) AutTo(ctOut, ct *Ciphertext, atk *AutKey, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.KeySwitchTo(ctOut, ct, atk.Value, outForm)

	op.plainOp.AutTo(ctOut.Body, ctOut.Body, atk.Idx)
	op.plainOp.AutTo(ctOut.Mask, ctOut.Mask, atk.Idx)
}

// AutHoisted returns ct(X^idx).
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) AutHoisted(ct *Ciphertext, ctMaskDcmp *Vector, atk *AutKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), 0, dft.FormCoeff)
	op.AutHoistedTo(ctOut, ct, ctMaskDcmp, atk, outForm)
	return ctOut
}

// AutHoistedTo computes ctOut = ct(X^idx).
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) AutHoistedTo(ctOut, ct *Ciphertext, ctMaskDcmp *Vector, atk *AutKey, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.KeySwitchHoistedTo(ctOut, ct, ctMaskDcmp, atk.Value, outForm)

	op.plainOp.AutTo(ctOut.Body, ctOut.Body, atk.Idx)
	op.plainOp.AutTo(ctOut.Mask, ctOut.Mask, atk.Idx)
}

// AutLazy returns ct(X^idx),
// but the output modulus includes the auxiliary modulus.
func (op *Operator) AutLazy(ct *Ciphertext, atk *AutKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), op.params.AuxModLen(), dft.FormCoeff)
	op.AutLazyTo(ctOut, ct, atk, outForm)
	return ctOut
}

// AutLazyTo computes ctOut = ct(X^idx),
// but the modulus of ctOut includes the auxiliary modulus.
func (op *Operator) AutLazyTo(ctOut, ct *Ciphertext, atk *AutKey, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.KeySwitchLazyTo(ctOut, ct, atk.Value, outForm)

	op.plainOp.AutTo(ctOut.Body, ctOut.Body, atk.Idx)
	op.plainOp.AutTo(ctOut.Mask, ctOut.Mask, atk.Idx)
}

// AutHoistedLazy returns ct(X^idx),
// but the output modulus includes the auxiliary modulus.
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) AutHoistedLazy(ct *Ciphertext, ctMaskDcmp *Vector, atk *AutKey, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ctMaskDcmp.AuxModLen(), dft.FormCoeff)
	op.AutHoistedLazyTo(ctOut, ct, ctMaskDcmp, atk, outForm)
	return ctOut
}

// AutHoistedLazyTo computes ctOut = ct(X^idx),
// but the modulus of ctOut includes the auxiliary modulus.
// ctMaskDcmp is a precomputed gadget decomposition of ct.Mask.
func (op *Operator) AutHoistedLazyTo(ctOut, ct *Ciphertext, ctMaskDcmp *Vector, atk *AutKey, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	op.KeySwitchHoistedLazyTo(ctOut, ct, ctMaskDcmp, atk.Value, outForm)

	op.plainOp.AutTo(ctOut.Body, ctOut.Body, atk.Idx)
	op.plainOp.AutTo(ctOut.Mask, ctOut.Mask, atk.Idx)
}
