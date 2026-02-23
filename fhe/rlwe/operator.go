package rlwe

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// PlainOperator evaluates operations over [*Element].
type PlainOperator struct {
	Params Parameters
	crtOp  crt.Operator
}

// NewPlainOperator creates a new [PlainOperator].
func NewPlainOperator(params Parameters) *PlainOperator {
	return &PlainOperator{
		Params: params,
		crtOp:  params.crtOp,
	}
}

// Add returns e0 + e1.
func (op *PlainOperator) Add(e0, e1 *Element) *Element {
	eOut := NewPolyCustom(e0.Rank(), e0.ModLen(), e0.HasAuxModulus(), e0.IsNTT())
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *PlainOperator) AddTo(eOut, e0, e1 *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e0, e1)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).AddTo(eOut.Value, e0.Value, e1.Value)
}

// Sub returns e0 - e1.
func (op *PlainOperator) Sub(e0, e1 *Element) *Element {
	eOut := NewPolyCustom(e0.Rank(), e0.ModLen(), e0.HasAuxModulus(), e0.IsNTT())
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *PlainOperator) SubTo(eOut, e0, e1 *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e0, e1)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).SubTo(eOut.Value, e0.Value, e1.Value)
}

// Neg returns -e.
func (op *PlainOperator) Neg(e *Element) *Element {
	eOut := NewPolyCustom(e.Rank(), e.ModLen(), e.HasAuxModulus(), e.IsNTT())
	op.NegTo(eOut, e)
	return eOut
}

// NegTo computes eOut = -e.
func (op *PlainOperator) NegTo(eOut, e *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).NegTo(eOut.Value, e.Value)
}

// Mul returns e0 * e1.
func (op *PlainOperator) Mul(e0, e1 *Element) *Element {
	eOut := NewPolyCustom(e0.Rank(), e0.ModLen(), e0.HasAuxModulus(), e0.IsNTT())
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
func (op *PlainOperator) MulTo(eOut, e0, e1 *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e0, e1)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).MulTo(eOut.Value, e0.Value, e1.Value)
}

// MulAdd returns eOut += e0 * e1.
func (op *PlainOperator) MulAdd(e0, e1 *Element) *Element {
	eOut := NewPolyCustom(e0.Rank(), e0.ModLen(), e0.HasAuxModulus(), e0.IsNTT())
	op.MulAddTo(eOut, e0, e1)
	return eOut
}

// MulAddTo computes eOut += e0 * e1.
func (op *PlainOperator) MulAddTo(eOut, e0, e1 *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e0, e1)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).MulAddTo(eOut.Value, e0.Value, e1.Value)
}

// MulSub returns eOut -= e0 * e1.
func (op *PlainOperator) MulSub(e0, e1 *Element) *Element {
	eOut := NewPolyCustom(e0.Rank(), e0.ModLen(), e0.HasAuxModulus(), e0.IsNTT())
	op.MulSubTo(eOut, e0, e1)
	return eOut
}

// MulSubTo computes eOut -= e0 * e1.
func (op *PlainOperator) MulSubTo(eOut, e0, e1 *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e0, e1)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).MulSubTo(eOut.Value, e0.Value, e1.Value)
}

// FwdNTTT returns FwdNTT(e).
func (op *PlainOperator) FwdNTT(e *Element) *Element {
	eOut := NewPolyCustom(e.Rank(), e.ModLen(), e.HasAuxModulus(), true)
	op.FwdNTTTo(eOut, e)
	return eOut
}

// FwdNTTTo computes FwdNTT(e) and stores the result in eOut.
func (op *PlainOperator) FwdNTTTo(eOut, e *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).FwdNTTTo(eOut.Value, e.Value)
}

// InvNTT returns InvNTT(e).
func (op *PlainOperator) InvNTT(e *Element) *Element {
	eOut := NewPolyCustom(e.Rank(), e.ModLen(), e.HasAuxModulus(), true)
	op.InvNTTTo(eOut, e)
	return eOut
}

// InvNTTTo computes InvNTT(e) and stores the result in eOut.
func (op *PlainOperator) InvNTTTo(eOut, e *Element) {
	checkOperable(len(op.Params.fullMod), eOut, e)

	op.subCRTOperator(eOut.ModLen(), eOut.HasAuxModulus()).InvNTTTo(eOut.Value, e.Value)
}

// subCRTOperator returns a [crt.Operator] for modulus up to given modulus length.
func (op *PlainOperator) subCRTOperator(modLen int, hasAux bool) crt.Operator {
	var idx []int
	if hasAux {
		idx = vec.Range(0, modLen)
	} else {
		idx = vec.Range(len(op.Params.auxMod), len(op.Params.auxMod)+modLen)
	}
	return op.crtOp.SubOperator(idx...)
}

// Operator evaluates operations over [*Ciphertext].
type Operator struct {
	Params Parameters

	plainOp *PlainOperator
	// TODO: We should eliminate the need for crtOp,
	// and replace it with plainOp instead
	crtOp crt.Operator
	dcmp  Decomposer

	pPool    *sync.Pool
	ctPool   *sync.Pool
	dcmpPool *sync.Pool
}

// NewOperator creates a new [Operator].
func NewOperator(params Parameters) *Operator {
	return &Operator{
		Params: params,

		plainOp: NewPlainOperator(params),
		crtOp:   params.crtOp,
		dcmp:    NewDecomposer(params),

		pPool: &sync.Pool{
			New: func() any {
				return crt.NewPoly(params.RingParams().Rank(), len(params.fullMod))
			},
		},
		ctPool: &sync.Pool{
			New: func() any {
				return NewCiphertext(params, params.HasAuxModulus(), true)
			},
		},
		dcmpPool: &sync.Pool{
			New: func() any {
				return NewVector(params, params.GadgetLen(), params.HasAuxModulus(), false)
			},
		},
	}
}

// AddPlain returns ct + pt.
func (op *Operator) AddPlain(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), pt.HasAuxModulus(), ct.IsNTT())
	op.AddPlainTo(ctOut, ct, pt)
	return ctOut
}

// AddPlainTo computes ctOut = ct + pt.
func (op *Operator) AddPlainTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.AddTo(ctOut.Body, ct.Body, pt)
	ctOut.Mask.CopyFrom(ct.Mask)
}

// SubPlain returns ct - pt.
func (op *Operator) SubPlain(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), pt.HasAuxModulus(), ct.IsNTT())
	op.SubPlainTo(ctOut, ct, pt)
	return ctOut
}

// SubPlainTo computes ctOut = ct - pt.
func (op *Operator) SubPlainTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.SubTo(ctOut.Body, ct.Body, pt)
	ctOut.Mask.CopyFrom(ct.Mask)
}

// MulPlain returns ct * pt.
func (op *Operator) MulPlain(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), pt.HasAuxModulus(), ct.IsNTT())
	op.MulPlainTo(ctOut, ct, pt)
	return ctOut
}

// MulPlainTo computes ctOut = ct * pt.
func (op *Operator) MulPlainTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.MulTo(ctOut.Body, ct.Body, pt)
	op.plainOp.MulTo(ctOut.Mask, ct.Mask, pt)
}

// MulAddPlain returns ct += pt.
func (op *Operator) MulAddPlain(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), pt.HasAuxModulus(), ct.IsNTT())
	op.MulAddPlainTo(ctOut, ct, pt)
	return ctOut
}

// MulAddPlainTo computes ctOut += ct * pt.
func (op *Operator) MulAddPlainTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.MulAddTo(ctOut.Body, ct.Body, pt)
	op.plainOp.MulAddTo(ctOut.Mask, ct.Mask, pt)
}

// MulSub returns ct0 - ct1.
func (op *Operator) MulSub(ct0, ct1 *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), ct1.HasAuxModulus(), ct1.IsNTT())
	op.MulSubTo(ctOut, ct0, ct1)
	return ctOut
}

// MulSubTo computes ctOut = ct0 - ct1.
func (op *Operator) MulSubTo(ctOut, ct0, ct1 *Ciphertext) {
	op.plainOp.MulSubTo(ctOut.Body, ct0.Body, ct1.Body)
	op.plainOp.MulSubTo(ctOut.Mask, ct0.Mask, ct1.Mask)
}

// Add returns ct0 + ct1.
func (op *Operator) Add(ct0, ct1 *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), ct1.HasAuxModulus(), ct1.IsNTT())
	op.AddTo(ctOut, ct0, ct1)
	return ctOut
}

// AddTo computes ctOut = ct0 + ct1.
func (op *Operator) AddTo(ctOut, ct0, ct1 *Ciphertext) {
	op.plainOp.AddTo(ctOut.Body, ct0.Body, ct1.Body)
	op.plainOp.AddTo(ctOut.Mask, ct0.Mask, ct1.Mask)
}

// Sub returns ct0 - ct1.
func (op *Operator) Sub(ct0, ct1 *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), ct1.HasAuxModulus(), ct1.IsNTT())
	op.SubTo(ctOut, ct0, ct1)
	return ctOut
}

// SubTo computes ctOut = ct0 - ct1.
func (op *Operator) SubTo(ctOut, ct0, ct1 *Ciphertext) {
	op.plainOp.SubTo(ctOut.Body, ct0.Body, ct1.Body)
	op.plainOp.SubTo(ctOut.Mask, ct0.Mask, ct1.Mask)
}

// TODO: FwdNTT/InvNTT of ciphertexts are very common, especially with isNTT flags.
// We should create a seperate structure for this. (something like `rlwe.Transformer`.)

// FwdNTT returns FwdNTT(ct).
func (op *Operator) FwdNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), ct.HasAuxModulus(), ct.IsNTT())
	op.FwdNTTTo(ctOut, ct)
	return ctOut
}

// FwdNTTTo computes ctOut = FwdNTT(ct).
func (op *Operator) FwdNTTTo(ctOut, ct *Ciphertext) {
	op.plainOp.FwdNTTTo(ctOut.Body, ct.Body)
	op.plainOp.FwdNTTTo(ctOut.Mask, ct.Mask)
}

// InvNTT returns InvNTT(ct).
func (op *Operator) InvNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), ct.HasAuxModulus(), ct.IsNTT())
	op.InvNTTTo(ctOut, ct)
	return ctOut
}

// InvNTTTo computes ctOut = InvNTT(ct).
func (op *Operator) InvNTTTo(ctOut, ct *Ciphertext) {
	op.plainOp.InvNTTTo(ctOut.Body, ct.Body)
	op.plainOp.InvNTTTo(ctOut.Mask, ct.Mask)
}

// DivByAuxModulus returns round(ct / AuxModulus).
func (op *Operator) DivByAuxModulus(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), isNTT, ct.IsNTT())
	op.DivByAuxModulusTo(ctOut, ct, isNTT)
	return ctOut
}

// DivByAuxModulusTo computes ctOut = round(ct / AuxModulus).
func (op *Operator) DivByAuxModulusTo(ctOut, ct *Ciphertext, isNTT bool) {
	if !op.Params.HasAuxModulus() || !ct.HasAuxModulus() {
		panic("input(s) must have auxiliary modulus")
	} else if ctOut.ModLen()+len(op.Params.auxMod) != ct.ModLen() {
		panic("inconsistent output")
	}

	baseLen, auxLen := ct.ModLen()-len(op.Params.auxMod), len(op.Params.auxMod)

	crtOpBase := op.crtOp.SubOperator(vec.Range(auxLen, auxLen+baseLen)...)
	crtOpAux := op.crtOp.SubOperator(vec.Range(0, auxLen)...)

	scaler := crt.NewScaler(op.Params.baseMod[:baseLen], op.Params.auxMod)

	auxInvBase := crt.NewScalar(baseLen) // auxMod^{-1} mod baseMod
	baseInvAux := crt.NewScalar(auxLen)  // baseMod^{-1} mod auxMod
	for i := 0; i < baseLen; i++ {
		for j := 0; j < auxLen; j++ {
			auxInvBase.Coeffs[i][0] = num.Mul(auxInvBase.Coeffs[i][0], num.Inv(op.Params.auxMod[j].Value(), op.Params.baseMod[i]), op.Params.baseMod[i])
			baseInvAux.Coeffs[j][0] = num.Mul(baseInvAux.Coeffs[j][0], num.Inv(op.Params.baseMod[i].Value(), op.Params.auxMod[j]), op.Params.auxMod[j])
		}
	}

	pDiv := op.pPool.Get().(*crt.Element)
	defer op.pPool.Put(pDiv)
	pDivBase := pDiv.WithModIdx(vec.Range(auxLen, auxLen+baseLen)...)
	pDivAux := pDiv.WithModIdx(vec.Range(0, auxLen)...)
	pDiv = pDiv.WithModIdx(vec.Range(0, baseLen+auxLen)...)

	pDiv.CopyFrom(ct.Body.Value)
	pDivBase.IsNTT = pDiv.IsNTT
	pDivAux.IsNTT = pDiv.IsNTT

	crtOpBase.MulTo(pDivBase, pDivBase, auxInvBase)
	crtOpAux.MulTo(pDivAux, pDivAux, baseInvAux)

	if pDivAux.IsNTT {
		crtOpAux.InvNTTTo(pDivAux, pDivAux)
	}

	scaler.ScaleTo(ctOut.Body.Value, pDivAux)
	ctOut.Body.Value.IsNTT = false

	if isNTT {
		crtOpBase.FwdNTTTo(ctOut.Body.Value, ctOut.Body.Value)
		if !pDivBase.IsNTT {
			crtOpBase.FwdNTTTo(pDivBase, pDivBase)
		}
	} else {
		if pDivBase.IsNTT {
			crtOpBase.InvNTTTo(pDivBase, pDivBase)
		}
	}
	crtOpBase.AddTo(ctOut.Body.Value, ctOut.Body.Value, pDivBase)

	pDiv.CopyFrom(ct.Mask.Value)
	pDivBase.IsNTT = pDiv.IsNTT
	pDivAux.IsNTT = pDiv.IsNTT

	crtOpBase.MulTo(pDivBase, pDivBase, auxInvBase)
	crtOpAux.MulTo(pDivAux, pDivAux, baseInvAux)

	if pDivAux.IsNTT {
		crtOpAux.InvNTTTo(pDivAux, pDivAux)
	}

	scaler.ScaleTo(ctOut.Mask.Value, pDivAux)
	ctOut.Mask.Value.IsNTT = false

	if isNTT {
		crtOpBase.FwdNTTTo(ctOut.Mask.Value, ctOut.Mask.Value)
		if !pDivBase.IsNTT {
			crtOpBase.FwdNTTTo(pDivBase, pDivBase)
		}
	} else {
		if pDivBase.IsNTT {
			crtOpBase.InvNTTTo(pDivBase, pDivBase)
		}
	}
	crtOpBase.AddTo(ctOut.Mask.Value, ctOut.Mask.Value, pDivBase)

	ctOut.Body.hasAux = false
	ctOut.Mask.hasAux = false
}

// Scale scales ct to l-th modulus.
//
//   - When l < ct.ModLen, it returns round(ct / Modulus[l:ct.ModLen]).
//   - When l > ct.ModLen, it returns ct * Modulus[l:ct.ModLen].
//   - When l == ct.ModLen, it returns a copy of ct.
func (op *Operator) Scale(ct *Ciphertext, l int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), l, false, ct.IsNTT())
	op.ScaleTo(ctOut, ct, l, isNTT)
	return ctOut
}

// ScaleTo scales ct to l-th modulus and writes the result to ctOut.
//
//   - When l < ct.ModLen, ctOut = round(ct / Modulus[l:ct.ModLen]).
//   - When l > ct.ModLen, ctOut = ct * Modulus[l:ct.ModLen].
//   - When l == ct.ModLen, ctOut = ct.
func (op *Operator) ScaleTo(ctOut, ct *Ciphertext, l int, isNTT bool) {
	if ct.HasAuxModulus() {
		panic("inconsistent input(s)")
	} else if ctOut.ModLen() != l {
		panic("inconsistent output")
	}

	inMod := op.Params.baseMod[:ct.ModLen()]
	outMod := op.Params.baseMod[:l]

	switch {
	case len(inMod) > len(outMod):
		scaleMod := inMod[l:]

		crtOpOutMod := op.crtOp.SubOperator(vec.Range(0, len(outMod))...)
		crtOpScale := op.crtOp.SubOperator(vec.Range(len(outMod), len(inMod))...)

		pDiv := op.pPool.Get().(*crt.Element)
		defer op.pPool.Put(pDiv)

		pDivInMod := pDiv.WithModIdx(vec.Range(0, len(inMod))...)
		pDivOutMod := pDiv.WithModIdx(vec.Range(0, len(outMod))...)
		pDivScaleMod := pDiv.WithModIdx(vec.Range(len(outMod), len(inMod))...)

		scaler := crt.NewScaler(outMod, inMod)

		scInvOutMod := crt.NewScalar(len(outMod)) // scale^{-1} mod outMod
		outModInvSc := crt.NewScalar(len(outMod)) // outMod^{-1} mod scale

		for i := range outMod {
			for j := range scaleMod {
				scInvOutMod.Coeffs[i][0] = num.Mul(scInvOutMod.Coeffs[i][0], num.Inv(scaleMod[j].Value(), outMod[i]), outMod[i])
				outModInvSc.Coeffs[j][0] = num.Mul(outModInvSc.Coeffs[j][0], num.Inv(outMod[i].Value(), scaleMod[j]), scaleMod[j])
			}
		}

		pDivInMod.CopyFrom(ct.Body.Value)
		pDivOutMod.IsNTT = pDivInMod.IsNTT
		pDivScaleMod.IsNTT = pDivInMod.IsNTT

		crtOpOutMod.MulTo(pDivOutMod, pDivOutMod, scInvOutMod)
		crtOpScale.MulTo(pDivScaleMod, pDivScaleMod, outModInvSc)
		if pDivScaleMod.IsNTT {
			crtOpScale.InvNTTTo(pDivScaleMod, pDivScaleMod)
		}

		scaler.ScaleTo(ctOut.Body.Value, pDivScaleMod)
		ctOut.Body.Value.IsNTT = false

		if isNTT {
			crtOpOutMod.FwdNTTTo(ctOut.Body.Value, ctOut.Body.Value)
			if !pDivOutMod.IsNTT {
				crtOpOutMod.FwdNTTTo(pDivOutMod, pDivOutMod)
			}
		} else {
			if pDivOutMod.IsNTT {
				crtOpOutMod.InvNTTTo(pDivOutMod, pDivOutMod)
			}
		}
		crtOpOutMod.AddTo(ctOut.Body.Value, ctOut.Body.Value, pDivOutMod)

		pDivInMod.CopyFrom(ct.Mask.Value)
		pDivOutMod.IsNTT = pDivInMod.IsNTT
		pDivScaleMod.IsNTT = pDivInMod.IsNTT

		crtOpOutMod.MulTo(pDivOutMod, pDivOutMod, scInvOutMod)
		crtOpScale.MulTo(pDivScaleMod, pDivScaleMod, outModInvSc)
		if pDivScaleMod.IsNTT {
			crtOpScale.InvNTTTo(pDivScaleMod, pDivScaleMod)
		}

		scaler.ScaleTo(ctOut.Mask.Value, pDivScaleMod)
		ctOut.Mask.Value.IsNTT = false

		if isNTT {
			crtOpOutMod.FwdNTTTo(ctOut.Mask.Value, ctOut.Mask.Value)
			if !pDivOutMod.IsNTT {
				crtOpOutMod.FwdNTTTo(pDivOutMod, pDivOutMod)
			}
		} else {
			if pDivOutMod.IsNTT {
				crtOpOutMod.InvNTTTo(pDivOutMod, pDivOutMod)
			}
		}
		crtOpOutMod.AddTo(ctOut.Mask.Value, ctOut.Mask.Value, pDivOutMod)

	case len(inMod) < len(outMod):
		scaleMod := outMod[:len(inMod)]
		scale := crt.NewScalar(len(inMod))
		for i := range inMod {
			for j := range scaleMod {
				scale.Coeffs[i][0] = num.Mul(scale.Coeffs[i][0], num.Inv(scaleMod[j].Value(), inMod[i]), inMod[i])
			}
		}

		ctOut.Body.Value.IsNTT = ct.Body.Value.IsNTT
		ctOut.Mask.Value.IsNTT = ct.Mask.Value.IsNTT

		for i := 0; i < len(inMod); i++ {
			vec.MulScalarTo(ctOut.Body.Value.Coeffs[i], ct.Body.Value.Coeffs[i], scale.Coeffs[i][0], outMod[i])
			vec.MulScalarTo(ctOut.Mask.Value.Coeffs[i], ct.Mask.Value.Coeffs[i], scale.Coeffs[i][0], outMod[i])
		}
		for i := len(inMod); i < len(outMod); i++ {
			clear(ctOut.Body.Value.Coeffs[i])
			clear(ctOut.Mask.Value.Coeffs[i])
		}

		if isNTT {
			if !ctOut.IsNTT() {
				op.plainOp.FwdNTTTo(ctOut.Body, ctOut.Body)
				op.plainOp.FwdNTTTo(ctOut.Mask, ctOut.Mask)
			}
		} else {
			if ctOut.IsNTT() {
				op.plainOp.InvNTTTo(ctOut.Body, ctOut.Body)
				op.plainOp.InvNTTTo(ctOut.Mask, ctOut.Mask)
			}
		}

	case ct.ModLen() == l:
		ctOut.CopyFrom(ct)

		if isNTT {
			if !ctOut.IsNTT() {
				op.plainOp.FwdNTTTo(ctOut.Body, ctOut.Body)
				op.plainOp.FwdNTTTo(ctOut.Mask, ctOut.Mask)
			}
		} else {
			if ctOut.IsNTT() {
				op.plainOp.InvNTTTo(ctOut.Body, ctOut.Body)
				op.plainOp.InvNTTTo(ctOut.Mask, ctOut.Mask)
			}
		}
	}
}

// HoistedGadgetProdLazy returns ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) HoistedGadgetProdLazy(pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.Params, op.Params.HasAuxModulus(), true)
	op.HoistedGadgetProdLazyTo(ctOut, pDcmp, ctGadEnc, isNTT)
	return ctOut
}

// HoistedGadgetProdLazyTo computes ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) HoistedGadgetProdLazyTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) {
	if op.Params.HasAuxModulus() {
		if !ctOut.HasAuxModulus() || !pDcmp.HasAuxModulus() || !ctGadEnc.HasAuxModulus() {
			panic("inconsistent input(s)")
		}
	} else {
		if ctOut.HasAuxModulus() || pDcmp.HasAuxModulus() || ctGadEnc.HasAuxModulus() {
			panic("inconsistent input(s)")
		}
	}

	baseLen, auxLen := pDcmp.ModLen()-len(op.Params.auxMod), len(op.Params.auxMod)
	if ctOut.ModLen() != baseLen+auxLen || pDcmp.Len() != op.dcmp.DecomposeLen(baseLen) {
		panic("inconsistent input(s)")
	}

	ctBuf := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(ctBuf)

	ctBuf.Clear()
	ctBuf.Body.hasAux = ctGadEnc.HasAuxModulus()
	ctBuf.Mask.hasAux = ctGadEnc.HasAuxModulus()

	for i := 0; i < pDcmp.Len(); i++ {
		op.plainOp.MulAddTo(ctBuf.Body, pDcmp.Value[i], ctGadEnc.Value[i].Body.WithModIdx(vec.Range(0, baseLen+auxLen)...))
		op.plainOp.MulAddTo(ctBuf.Mask, pDcmp.Value[i], ctGadEnc.Value[i].Mask.WithModIdx(vec.Range(0, baseLen+auxLen)...))
	}

	if !isNTT {
		op.InvNTTTo(ctBuf, ctBuf)
	}

	ctOut.CopyFrom(ctBuf)
}

// HoistedGadgetProd returns ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
func (op *Operator) HoistedGadgetProd(pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.Params, op.Params.HasAuxModulus(), true)
	op.HoistedGadgetProdTo(ctOut, pDcmp, ctGadEnc, isNTT)
	return ctOut
}

// HoistedGadgetProdTo computes ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
func (op *Operator) HoistedGadgetProdTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) {
	ctOutFullMod := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(ctOutFullMod)

	op.HoistedGadgetProdLazyTo(ctOutFullMod, pDcmp, ctGadEnc, true)

	if op.Params.HasAuxModulus() {
		op.DivByAuxModulusTo(ctOut, ctOutFullMod, isNTT)
	} else {
		ctOut.CopyFrom(ctOutFullMod)
		if !isNTT {
			op.InvNTTTo(ctOut, ctOut)
		}
	}
}

// GadgetProdLazy returns ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProdLazy(p *Element, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.Params, op.Params.HasAuxModulus(), true)
	op.GadgetProdLazyTo(ctOut, p, ctGadEnc, isNTT)
	return ctOut
}

// GadgetProdLazyTo computes ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProdLazyTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, isNTT bool) {
	if p.HasAuxModulus() {
		panic("inconsistent input(s)")
	}

	modLen, auxLen := p.ModLen(), len(op.Params.auxMod)

	pInvNTT := op.pPool.Get().(*Element)
	pInvNTT.CopyFrom(p)

	pInvNTT = pInvNTT.WithModIdx(vec.Range(0, modLen)...)
	if p.Value.IsNTT {
		op.plainOp.InvNTTTo(pInvNTT, p)
	} else {
		pInvNTT.CopyFrom(p)
	}

	pDcmp := op.dcmpPool.Get().(*Vector)
	defer op.dcmpPool.Put(pDcmp)
	pDcmp = pDcmp.WithLenModIdx(op.dcmp.DecomposeLen(modLen), vec.Range(0, modLen+auxLen)...)
	op.dcmp.DecomposeTo(pDcmp, pInvNTT)

	for i := range pDcmp.Value {
		op.plainOp.FwdNTTTo(pDcmp.Value[i], pDcmp.Value[i])
	}

	op.HoistedGadgetProdLazyTo(ctOut, pDcmp, ctGadEnc, isNTT)
}

// GadgetProd returns ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProd(p *Element, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.Params, op.Params.HasAuxModulus(), true)
	op.GadgetProdTo(ctOut, p, ctGadEnc, isNTT)
	return ctOut
}

// GadgetProdTo computes ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProdTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, isNTT bool) {
	if p.HasAuxModulus() {
		panic("inconsistent input(s)")
	}

	modLen, auxLen := p.ModLen(), len(op.Params.auxMod)

	pInvNTT := op.pPool.Get().(*Element)
	pInvNTT.CopyFrom(p)

	pInvNTT = pInvNTT.WithModIdx(vec.Range(0, modLen)...)
	if p.Value.IsNTT {
		op.plainOp.InvNTTTo(pInvNTT, p)
	} else {
		pInvNTT.CopyFrom(p)
	}

	pDcmp := op.dcmpPool.Get().(*Vector)
	defer op.dcmpPool.Put(pDcmp)
	pDcmp = pDcmp.WithLenModIdx(op.dcmp.DecomposeLen(modLen), vec.Range(0, modLen+auxLen)...)
	op.dcmp.DecomposeTo(pDcmp, pInvNTT)

	for i := range pDcmp.Value {
		op.plainOp.FwdNTTTo(pDcmp.Value[i], pDcmp.Value[i])
	}

	op.HoistedGadgetProdTo(ctOut, pDcmp, ctGadEnc, isNTT)
}

// // Relin performs a relinearisation and returns the result.
// func (o *Operator) Relin(cIn *Tensor, rlk *RelinKey, isNTT bool) *Ciphertext {
// 	modLen := cIn.ModLen()
// 	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
// 	o.RelinTo(cOut, cIn, rlk, isNTT)
// 	return cOut
// }

// // RelinTo performs a relinearisation and stores the result in cOut.
// // Revise to deal with the case where the input and output share the same elements.
// func (o *Operator) RelinTo(cOut *Ciphertext, cIn *Tensor, rlk *RelinKey, isNTT bool) {
// 	if cIn.HasAux {
// 		panic("Ciphertext must not have auxiliary modulus.")
// 	} else if cIn.Degree() != 3 {
// 		panic("Ciphertext must have degree 3.")
// 	}

// 	modLen := cIn.ModLen()
// 	eval := o.SubEvaluatorAt(false, modLen)
// 	buf := o.buf.pNTT.WithModIdx(vec.Range(0, modLen)...)

// 	o.GadgetProdTo(cOut, &PlainPoly{Value: cIn.Value[2]}, &GadgetEncryption{Value: rlk.Value}, isNTT)

// 	if cIn.Value[0].IsNTT && !isNTT {
// 		eval.InvNTTTo(buf, cIn.Value[0])
// 	} else if !cIn.Value[0].IsNTT && isNTT {
// 		eval.FwdNTTTo(buf, cIn.Value[0])
// 	} else {
// 		buf.CopyFrom(cIn.Value[0])
// 	}

// 	eval.AddTo(cOut.Body, cOut.Body, buf)

// 	if cIn.Value[1].IsNTT && !isNTT {
// 		eval.InvNTTTo(buf, cIn.Value[1])
// 	} else if !cIn.Value[1].IsNTT && isNTT {
// 		eval.FwdNTTTo(buf, cIn.Value[1])
// 	} else {
// 		buf.CopyFrom(cIn.Value[1])
// 	}

// 	eval.AddTo(cOut.Mask, cOut.Mask, buf)
// }

// // KeySwitch performs a key switch and returns the result.
// func (o *Operator) KeySwitch(cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) *Ciphertext {
// 	if cIn.HasAux {
// 		panic("Ciphertext must not have auxiliary modulus.")
// 	}

// 	modLen := cIn.ModLen()
// 	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
// 	o.KeySwitchTo(cOut, cIn, ksk, isNTT)
// 	return cOut
// }

// // KeySwitchTo performs a key switch and stores the result in cOut.
// func (o *Operator) KeySwitchTo(cOut *Ciphertext, cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) {
// 	if cIn.HasAux {
// 		panic("Ciphertext must not have auxiliary modulus.")
// 	}

// 	modLen := cIn.ModLen()
// 	eval := o.SubEvaluatorAt(false, modLen)
// 	buf := o.buf.pKsw.WithModIdx(vec.Range(0, modLen)...)
// 	if cIn.Mask.IsNTT {
// 		eval.InvNTTTo(buf, cIn.Mask)
// 	} else {
// 		buf.CopyFrom(cIn.Mask)
// 	}

// 	var auxLen int
// 	if o.Params.auxModulus != nil {
// 		auxLen = len(o.Params.auxModulus)
// 	}
// 	evalAux := o.SubEvaluatorAt(true, modLen+auxLen)

// 	decmpLen := o.Decmp.DecomposeLen(modLen)
// 	decmp := o.buf.decmp.WithDegreeAndModIdx(decmpLen, vec.Range(0, modLen+auxLen)...)
// 	o.Decmp.DecomposeTo(decmp, buf)
// 	for i := 0; i < decmpLen; i++ {
// 		evalAux.FwdNTTTo(decmp.Value[i], decmp.Value[i])
// 	}

// 	o.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)
// }

// // HoistedKeySwitch performs a hoisted key switch and returns the result.
// func (o *Operator) HoistedKeySwitch(decmp *Tensor, cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) *Ciphertext {
// 	modLen := cIn.ModLen()
// 	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
// 	o.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)
// 	return cOut
// }

// // HoistedKeySwitchTo performs a hoisted key switch and stores the result in cOut.
// func (o *Operator) HoistedKeySwitchTo(cOut *Ciphertext, decmp *Tensor, cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) {
// 	if cIn.HasAux {
// 		panic("Ciphertext must not have auxiliary modulus.")
// 	}

// 	modLen := cIn.ModLen()
// 	eval := o.SubEvaluatorAt(false, modLen)
// 	buf := o.buf.pKsw.WithModIdx(vec.Range(0, modLen)...)

// 	if cIn.Body.IsNTT && !isNTT {
// 		eval.InvNTTTo(buf, cIn.Body)
// 	} else if !cIn.Body.IsNTT && isNTT {
// 		eval.FwdNTTTo(buf, cIn.Body)
// 	} else {
// 		buf.CopyFrom(cIn.Body)
// 	}

// 	kskGad := &GadgetEncryption{Value: ksk.Value}
// 	o.HoistedGadgetProdTo(cOut, decmp, kskGad, isNTT)
// 	eval.AddTo(cOut.Body, cOut.Body, buf)
// }

// // Aut performs an automorphism and returns the result.
// func (o *Operator) Aut(cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) *Ciphertext {
// 	modLen := cIn.ModLen()
// 	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
// 	o.AutTo(cOut, cIn, atk, isNTT)
// 	return cOut
// }

// // AutTo performs an automorphism and stores the result in cOut.
// func (o *Operator) AutTo(cOut *Ciphertext, cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) {
// 	modLen := cIn.ModLen()
// 	ksk := &KeySwitchKey{Value: atk.Value}
// 	o.KeySwitchTo(cOut, cIn, ksk, isNTT)

// 	eval := o.SubEvaluatorAt(false, modLen)

// 	eval.AutTo(cOut.Body, cOut.Body, atk.Idx)
// 	eval.AutTo(cOut.Mask, cOut.Mask, atk.Idx)
// }

// // HoistedAut performs a hoisted automorphism and returns the result.
// func (o *Operator) HoistedAut(decmp *Tensor, cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) *Ciphertext {
// 	modLen := cIn.ModLen()
// 	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
// 	o.HoistedAutTo(cOut, decmp, cIn, atk, isNTT)
// 	return cOut
// }

// // HoistedAutTo performs a hoisted automorphism and stores the result in cOut.
// func (o *Operator) HoistedAutTo(cOut *Ciphertext, decmp *Tensor, cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) {
// 	if cIn.HasAux {
// 		panic("Ciphertext must not have auxiliary modulus.")
// 	}

// 	modLen := cIn.ModLen()
// 	eval := o.SubEvaluatorAt(false, modLen)

// 	ksk := &KeySwitchKey{Value: atk.Value}
// 	o.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)
// 	eval.AutTo(cOut.Body, cOut.Body, atk.Idx)
// 	eval.AutTo(cOut.Mask, cOut.Mask, atk.Idx)
// }

// // ExtProd performs an external product and returns the result.
// func (o *Operator) ExtProd(cIn *Ciphertext, gsw *RGSW, isNTT bool) *Ciphertext {
// 	modLen := cIn.ModLen()
// 	cOut := NewCiphertextCustom(o.Params.ringParams.Rank(), modLen, 0, false)
// 	o.ExtProdTo(cOut, cIn, gsw, isNTT)
// 	return cOut
// }

// // ExtProdTo performs an external product and stores the result in cOut.
// func (o *Operator) ExtProdTo(cOut *Ciphertext, cIn *Ciphertext, gsw *RGSW, isNTT bool) {
// 	if cIn.HasAux {
// 		panic("Ciphertext must not have auxiliary modulus.")
// 	}

// 	modLen := cIn.ModLen()
// 	buf := o.buf.ctExt.WithModIdx(vec.Range(0, modLen)...)

// 	o.GadgetProdTo(buf, &PlainPoly{Value: cIn.Body}, gsw.Body, isNTT)
// 	o.GadgetProdTo(cOut, &PlainPoly{Value: cIn.Mask}, gsw.Mask, isNTT)

// 	o.AddTo(cOut, cOut, buf)
// }

// // HoistedExtProdTo performs a hoisted external product and stores the result in cOut.
// func (o *Operator) HoistedExtProdTo(cOut *Ciphertext, decmpBody *Tensor, decmpMask *Tensor, gsw *RGSW, isNTT bool) {
// 	modLen := decmpBody.ModLen()
// 	if o.Params.auxModulus != nil {
// 		modLen -= len(o.Params.auxModulus)
// 	}
// 	buf := o.buf.ctExt.WithModIdx(vec.Range(0, modLen)...)

// 	o.HoistedGadgetProdLazyTo(buf, decmpBody, gsw.Body, isNTT)
// 	o.HoistedGadgetProdLazyTo(cOut, decmpMask, gsw.Mask, isNTT)

// 	o.AddTo(cOut, cOut, buf)
// }
