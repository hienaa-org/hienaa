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
func (op *Operator) Mul(ct *Ciphertext, e *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	op.MulTo(ctOut, ct, e)
	return ctOut
}

// MulTo computes ctOut = ct * e.
func (op *Operator) MulTo(ctOut, ct *Ciphertext, e *Element) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	switch {
	case e.Type() == crt.TypeScalar:
		op.plainOp.MulTo(ctOut.Body, ct.Body, e)
		op.plainOp.MulTo(ctOut.Mask, ct.Mask, e)

	case e.Type() == crt.TypePoly:
		eBufFull := op.poolPoly.Get()
		defer op.poolPoly.Put(eBufFull)

		eBuf := eBufFull.WithModLen(e.BaseModLen(), e.AuxModLen())
		eBuf.CopyFrom(e)

		op.plainOp.MulTo(ctOut.Body, ct.Body, eBuf)
		op.plainOp.MulTo(ctOut.Mask, ct.Mask, eBuf)
	}
}

// MulAddTo computes ctOut += ct * e.
func (op *Operator) MulAddTo(ctOut, ct *Ciphertext, e *Element) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	switch {
	case e.Type() == crt.TypeScalar:
		op.plainOp.MulAddTo(ctOut.Body, ct.Body, e)
		op.plainOp.MulAddTo(ctOut.Mask, ct.Mask, e)

	case e.Type() == crt.TypePoly:
		eBufFull := op.poolPoly.Get()
		defer op.poolPoly.Put(eBufFull)

		eBuf := eBufFull.WithModLen(e.BaseModLen(), e.AuxModLen())
		eBuf.CopyFrom(e)

		op.plainOp.MulAddTo(ctOut.Body, ct.Body, eBuf)
		op.plainOp.MulAddTo(ctOut.Mask, ct.Mask, eBuf)
	}
}

// MulSubTo computes ctOut -= ct * e.
func (op *Operator) MulSubTo(ctOut, ct *Ciphertext, e *Element) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)

	switch {
	case e.Type() == crt.TypeScalar:
		op.plainOp.MulSubTo(ctOut.Body, ct.Body, e)
		op.plainOp.MulSubTo(ctOut.Mask, ct.Mask, e)

	case e.Type() == crt.TypePoly:
		eBufFull := op.poolPoly.Get()
		defer op.poolPoly.Put(eBufFull)

		eBuf := eBufFull.WithModLen(e.BaseModLen(), e.AuxModLen())
		eBuf.CopyFrom(e)

		op.plainOp.MulSubTo(ctOut.Body, ct.Body, eBuf)
		op.plainOp.MulSubTo(ctOut.Mask, ct.Mask, eBuf)
	}
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

// GadProd returns p * ctGadEnc.
func (op *Operator) GadProd(p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), p.BaseModLen(), 0, dft.FormCoeff)
	op.GadProdTo(ctOut, p, ctGadEnc, outForm)
	return ctOut
}

// GadProdTo computes ctOut = p * ctGadEnc.
func (op *Operator) GadProdTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	if ctOut.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctOutAuxFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutAuxFull)

	ctOutAux := ctOutAuxFull.WithModLen(ctOut.BaseModLen(), op.params.AuxModLen())
	op.GadProdLazyTo(ctOutAux, p, ctGadEnc, outForm)

	if op.params.AuxModLen() > 0 {
		op.ScaleTo(ctOut, ctOutAux, outForm)
		return
	}

	ctOut.CopyFrom(ctOutAux)
}

// GadProdHoisted returns pDcmp * ctGadEnc.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadProdHoisted(pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), pDcmp.BaseModLen(), 0, dft.FormCoeff)
	op.GadProdHoistedTo(ctOut, pDcmp, ctGadEnc, outForm)
	return ctOut
}

// GadProdHoistedTo computes ctOut = pDcmp * ctGadEnc.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadProdHoistedTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	if ctOut.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctOutAuxFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutAuxFull)

	ctOutAux := ctOutAuxFull.WithModLen(ctOut.BaseModLen(), pDcmp.AuxModLen())
	op.GadProdHoistedLazyTo(ctOutAux, pDcmp, ctGadEnc, outForm)

	if op.params.AuxModLen() > 0 {
		op.ScaleTo(ctOut, ctOutAux, outForm)
		return
	}

	ctOut.CopyFrom(ctOutAux)
}

// GadProdLazy returns p * ctGadEnc,
// but the output modulus includes the auxiliary modulus.
func (op *Operator) GadProdLazy(p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), p.BaseModLen(), op.params.AuxModLen(), dft.FormCoeff)
	op.GadProdLazyTo(ctOut, p, ctGadEnc, outForm)
	return ctOut
}

// GadProdLazyTo computes ctOut = p * ctGadEnc,
// but the modulus of ctOut includes the auxiliary modulus.
func (op *Operator) GadProdLazyTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, outForm dft.Form) {
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
	op.GadProdHoistedLazyTo(ctOut, pDcmp, ctGadEnc, outForm)
}

// GadProdHoistedLazy returns pDcmp * ctGadEnc,
// but the output modulus includes the auxiliary modulus.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadProdHoistedLazy(pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), pDcmp.BaseModLen(), pDcmp.AuxModLen(), dft.FormCoeff)
	op.GadProdHoistedLazyTo(ctOut, pDcmp, ctGadEnc, outForm)
	return ctOut
}

// GadProdHoistedLazyTo computes ctOut = pDcmp * ctGadEnc,
// but the modulus of ctOut includes the auxiliary modulus.
// pDcmp is a precomputed gadget decomposition.
func (op *Operator) GadProdHoistedLazyTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, outForm dft.Form) {
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

	op.GadProdTo(ctOut, ct.Mask, (*GadgetEncryption)(ksk), outForm)
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

	op.GadProdHoistedTo(ctOut, ctMaskDcmp, (*GadgetEncryption)(ksk), outForm)
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

	op.GadProdLazyTo(ctOut, ct.Mask, (*GadgetEncryption)(ksk), outForm)
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

	op.GadProdHoistedLazyTo(ctOut, ctMaskDcmp, (*GadgetEncryption)(ksk), outForm)
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

// ExtProd returns ct * ctRGSW.
// ct must not have auxiliary modulus, and ctRGSW must be in NTT form.
func (op *Operator) ExtProd(ct *Ciphertext, ctRGSW *RGSWEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), 0, dft.FormCoeff)
	op.ExtProdTo(ctOut, ct, ctRGSW, outForm)
	return ctOut
}

// ExtProdTo computes ctOut = ct * ctRGSW.
func (op *Operator) ExtProdTo(ctOut, ct *Ciphertext, ctRGSW *RGSWEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	if ctOut.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctOutAuxFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutAuxFull)

	ctOutAux := ctOutAuxFull.WithModLen(ctOut.BaseModLen(), op.params.AuxModLen())
	op.ExtProdLazyTo(ctOutAux, ct, ctRGSW, outForm)

	if op.params.AuxModLen() > 0 {
		op.ScaleTo(ctOut, ctOutAux, outForm)
		return
	}

	ctOut.CopyFrom(ctOutAux)
}

// ExtProdHoisted returns ct * ctRGSW.
// ctBodyDcmp and ctMaskDcmp are precomputed gadget decompositions of ct.Body and ct.Mask in NTT form.
func (op *Operator) ExtProdHoisted(ctBodyDcmp, ctMaskDcmp *Vector, ctRGSW *RGSWEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ctBodyDcmp.BaseModLen(), 0, dft.FormCoeff)
	op.ExtProdHoistedTo(ctOut, ctBodyDcmp, ctMaskDcmp, ctRGSW, outForm)
	return ctOut
}

// ExtProdHoistedTo computes ctOut = ct * ctRGSW.
// ctBodyDcmp and ctMaskDcmp are precomputed gadget decompositions of ct.Body and ct.Mask in NTT form.
func (op *Operator) ExtProdHoistedTo(ctOut *Ciphertext, ctBodyDcmp, ctMaskDcmp *Vector, ctRGSW *RGSWEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	if ctOut.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctOutAuxFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctOutAuxFull)

	ctOutAux := ctOutAuxFull.WithModLen(ctOut.BaseModLen(), ctBodyDcmp.AuxModLen())
	op.ExtProdHoistedLazyTo(ctOutAux, ctBodyDcmp, ctMaskDcmp, ctRGSW, outForm)

	if op.params.AuxModLen() > 0 {
		op.ScaleTo(ctOut, ctOutAux, outForm)
		return
	}

	ctOut.CopyFrom(ctOutAux)
}

// ExtProdLazy returns ct * ctRGSW,
// but the output modulus includes the auxiliary modulus.
// ct must not have auxiliary modulus, and ctRGSW must be in NTT form.
func (op *Operator) ExtProdLazy(ct *Ciphertext, ctRGSW *RGSWEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.BaseModLen(), op.params.AuxModLen(), dft.FormCoeff)
	op.ExtProdLazyTo(ctOut, ct, ctRGSW, outForm)
	return ctOut
}

// ExtProdLazyTo computes ctOut = ct * ctRGSW,
// but the modulus of ctOut includes the auxiliary modulus.
func (op *Operator) ExtProdLazyTo(ctOut, ct *Ciphertext, ctRGSW *RGSWEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ct)
	checkCiphertextShape(op.params, ctRGSW)
	if ct.AuxModLen() > 0 {
		panic("inconsistent inputs")
	}

	ctBodyDcmpFull := op.poolPolyDcmp.Get()
	defer op.poolPolyDcmp.Put(ctBodyDcmpFull)
	ctMaskDcmpFull := op.poolPolyDcmp.Get()
	defer op.poolPolyDcmp.Put(ctMaskDcmpFull)

	ctBodyDcmp := ctBodyDcmpFull.Slice(0, op.dcmp.DecomposeLen(ct.BaseModLen())).WithModLen(ct.BaseModLen(), op.params.AuxModLen())
	ctMaskDcmp := ctMaskDcmpFull.Slice(0, op.dcmp.DecomposeLen(ct.BaseModLen())).WithModLen(ct.BaseModLen(), op.params.AuxModLen())

	op.dcmp.DecomposeTo(ctBodyDcmp, ct.Body, dft.FormNTT)
	op.dcmp.DecomposeTo(ctMaskDcmp, ct.Mask, dft.FormNTT)
	op.ExtProdHoistedLazyTo(ctOut, ctBodyDcmp, ctMaskDcmp, ctRGSW, outForm)
}

// ExtProdHoistedLazy returns ct * ctRGSW,
// but the output modulus includes the auxiliary modulus.
// ctBodyDcmp and ctMaskDcmp are precomputed gadget decompositions of ct.Body and ct.Mask in NTT form.
func (op *Operator) ExtProdHoistedLazy(ctBodyDcmp, ctMaskDcmp *Vector, ctRGSW *RGSWEncryption, outForm dft.Form) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ctBodyDcmp.BaseModLen(), ctBodyDcmp.AuxModLen(), dft.FormCoeff)
	op.ExtProdHoistedLazyTo(ctOut, ctBodyDcmp, ctMaskDcmp, ctRGSW, outForm)
	return ctOut
}

// ExtProdHoistedLazyTo computes ctOut = ct * ctRGSW,
// but the modulus of ctOut includes the auxiliary modulus.
// ctBodyDcmp and ctMaskDcmp are precomputed gadget decompositions of ct.Body and ct.Mask in NTT form.
func (op *Operator) ExtProdHoistedLazyTo(ctOut *Ciphertext, ctBodyDcmp, ctMaskDcmp *Vector, ctRGSW *RGSWEncryption, outForm dft.Form) {
	checkCiphertextShape(op.params, ctOut)
	checkCiphertextShape(op.params, ctRGSW)
	if ctBodyDcmp.Type() != crt.TypePoly || ctMaskDcmp.Type() != crt.TypePoly {
		panic("inconsistent types")
	}

	ctBodyFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctBodyFull)
	ctMaskFull := op.poolRLWE.Get()
	defer op.poolRLWE.Put(ctMaskFull)

	ctBody := ctBodyFull.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen())
	ctMask := ctMaskFull.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen())
	op.GadProdHoistedLazyTo(ctBody, ctBodyDcmp, ctRGSW.BodyStack, dft.FormNTT)
	op.GadProdHoistedLazyTo(ctMask, ctMaskDcmp, ctRGSW.MaskStack, dft.FormNTT)

	op.AddTo(ctOut, ctBody, ctMask)
	if outForm == dft.FormCoeff {
		op.InvNTTTo(ctOut, ctOut)
	}
}
