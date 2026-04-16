package rlwe

import (
	"math/big"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// PlainOperator evaluates operations over [*Element].
type PlainOperator struct {
	Params Parameters
	crtOp  *crt.Operator

	pPool *pool.Pool[*crt.Element]
}

// NewPlainOperator creates a new [PlainOperator].
func NewPlainOperator(params Parameters) *PlainOperator {
	return &PlainOperator{
		Params: params,
		crtOp:  params.crtOp,

		pPool: pool.NewPool(func() *crt.Element {
			return crt.NewPoly(params.Rank(), len(params.fullMod))
		}),
	}
}

// Add returns e0 + e1.
func (op *PlainOperator) Add(e0, e1 *Element) *Element {
	eOut := NewElement(e0.Rank(), e0.BaseModLen(), e0.AuxModLen(), e0.IsNTT())
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *PlainOperator) AddTo(eOut, e0, e1 *Element) {
	checkOperable(eOut, e0, e1)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).AddTo(eOut.Value, e0.Value, e1.Value)
}

// Sub returns e0 - e1.
func (op *PlainOperator) Sub(e0, e1 *Element) *Element {
	eOut := NewElement(e0.Rank(), e0.BaseModLen(), e0.AuxModLen(), e0.IsNTT())
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *PlainOperator) SubTo(eOut, e0, e1 *Element) {
	checkOperable(eOut, e0, e1)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).SubTo(eOut.Value, e0.Value, e1.Value)
}

// Neg returns -e.
func (op *PlainOperator) Neg(e *Element) *Element {
	eOut := NewElement(e.Rank(), e.BaseModLen(), e.AuxModLen(), e.IsNTT())
	op.NegTo(eOut, e)
	return eOut
}

// NegTo computes eOut = -e.
func (op *PlainOperator) NegTo(eOut, e *Element) {
	checkOperable(eOut, e)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).NegTo(eOut.Value, e.Value)
}

// Mul returns e0 * e1.
func (op *PlainOperator) Mul(e0, e1 *Element) *Element {
	eOut := NewElement(e0.Rank(), e0.BaseModLen(), e0.AuxModLen(), e0.IsNTT())
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
func (op *PlainOperator) MulTo(eOut, e0, e1 *Element) {
	checkOperable(eOut, e0, e1)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).MulTo(eOut.Value, e0.Value, e1.Value)
}

// MulAdd returns eOut += e0 * e1.
func (op *PlainOperator) MulAdd(e0, e1 *Element) *Element {
	eOut := NewElement(e0.Rank(), e0.BaseModLen(), e0.AuxModLen(), e0.IsNTT())
	op.MulAddTo(eOut, e0, e1)
	return eOut
}

// MulAddTo computes eOut += e0 * e1.
func (op *PlainOperator) MulAddTo(eOut, e0, e1 *Element) {
	checkOperable(eOut, e0, e1)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).MulAddTo(eOut.Value, e0.Value, e1.Value)
}

// MulSub returns eOut -= e0 * e1.
func (op *PlainOperator) MulSub(e0, e1 *Element) *Element {
	eOut := NewElement(e0.Rank(), e0.BaseModLen(), e0.AuxModLen(), e0.IsNTT())
	op.MulSubTo(eOut, e0, e1)
	return eOut
}

// MulSubTo computes eOut -= e0 * e1.
func (op *PlainOperator) MulSubTo(eOut, e0, e1 *Element) {
	checkOperable(eOut, e0, e1)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).MulSubTo(eOut.Value, e0.Value, e1.Value)
}

// FwdNTTT returns FwdNTT(e).
func (op *PlainOperator) FwdNTT(e *Element) *Element {
	eOut := NewElement(e.Rank(), e.BaseModLen(), e.AuxModLen(), true)
	op.FwdNTTTo(eOut, e)
	return eOut
}

// FwdNTTTo computes FwdNTT(e) and stores the result in eOut.
func (op *PlainOperator) FwdNTTTo(eOut, e *Element) {
	checkOperable(eOut, e)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).FwdNTTTo(eOut.Value, e.Value)
}

// InvNTT returns InvNTT(e).
func (op *PlainOperator) InvNTT(e *Element) *Element {
	eOut := NewElement(e.Rank(), e.BaseModLen(), e.AuxModLen(), true)
	op.InvNTTTo(eOut, e)
	return eOut
}

// InvNTTTo computes InvNTT(e) and stores the result in eOut.
func (op *PlainOperator) InvNTTTo(eOut, e *Element) {
	checkOperable(eOut, e)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).InvNTTTo(eOut.Value, e.Value)
}

// CanAut returns whether the given automorphism index is valid.
func (op *PlainOperator) CanAut(idx int) bool {
	return op.crtOp.CanAut(idx)
}

// Aut returns aut(e, idx).
func (op *PlainOperator) Aut(e *Element, idx int) *Element {
	eOut := NewElement(e.Rank(), e.BaseModLen(), e.AuxModLen(), e.IsNTT())
	op.AutTo(eOut, e, idx)
	return eOut
}

// AutTo computes aut(e, idx) and stores the result in eOut.
func (op *PlainOperator) AutTo(eOut, e *Element, idx int) {
	checkOperable(eOut, e)
	op.withModIdxCRT(eOut.BaseModLen(), eOut.AuxModLen()).AutTo(eOut.Value, e.Value, idx)
}

// withModIdxCRT returns a [crt.Operator] for modulus up to given modulus length.
func (op *PlainOperator) withModIdxCRT(baseLen, auxLen int) *crt.Operator {
	var idx []int

	if auxLen > 0 {
		idx = vec.Range(len(op.Params.auxMod)-auxLen, len(op.Params.auxMod)+baseLen)
	} else {
		idx = vec.Range(len(op.Params.auxMod), len(op.Params.auxMod)+baseLen)
	}
	return op.crtOp.WithModIdx(idx...)
}

// AsBig returns e as *[big.Int] vector.
func (op *PlainOperator) AsBig(e *Element) []*big.Int {
	baseLen := e.BaseModLen()
	auxLen := e.AuxModLen()

	paramAuxLen := len(op.Params.auxMod)

	lo := paramAuxLen - auxLen
	hi := paramAuxLen + baseLen

	crtOp := op.crtOp.WithModIdx(vec.Range(lo, hi)...)

	return crtOp.AsBig(e.Value)
}

// ModRaise raises e to the modulus with length l.
func (op *PlainOperator) ModRaise(e *Element, l int, isNTT bool) *Element {
	eOut := NewElement(e.Rank(), l, 0, isNTT)
	op.ModRaiseTo(eOut, e, isNTT)
	return eOut
}

// ModRaiseTo raises e to the modulus with length l and stores the result in eOut.
func (op *PlainOperator) ModRaiseTo(eOut, e *Element, isNTT bool) {
	if e.AuxModLen() > 0 || eOut.AuxModLen() > 0 {
		panic("input(s) should not have auxiliary modulus")
	} else if e.BaseModLen() > eOut.BaseModLen() {
		panic("inconsistent output")
	}

	inLen := e.BaseModLen()
	outLen := eOut.BaseModLen()

	switch {
	case inLen < outLen:
		if e.IsNTT() && isNTT {
			emb := crt.NewEmbedder(op.Params.baseMod[inLen:outLen], op.Params.baseMod[:inLen])

			pNTT := op.pPool.Get()
			defer op.pPool.Put(pNTT)
			pNTT = pNTT.WithModIdx(vec.Range(0, inLen)...)

			crtOpIn := op.crtOp.WithModIdx(vec.Range(0, inLen)...)
			crtOpDiff := op.crtOp.WithModIdx(vec.Range(inLen, outLen)...)
			crtOpIn.InvNTTTo(pNTT, e.Value)

			eOutIn := &crt.Element{
				Coeffs: eOut.Value.Coeffs[:inLen],
				IsNTT:  false,
			}
			eOutDiff := &crt.Element{
				Coeffs: eOut.Value.Coeffs[inLen:outLen],
				IsNTT:  false,
			}

			eOutIn.CopyFrom(e.Value)
			emb.EmbedTo(eOutDiff, pNTT)
			crtOpDiff.FwdNTTTo(eOutDiff, eOutDiff)

			eOut.Value.IsNTT = true
		} else {
			emb := crt.NewEmbedder(op.Params.baseMod[:outLen], op.Params.baseMod[:inLen])

			if e.IsNTT() {
				pNTT := op.pPool.Get()
				defer op.pPool.Put(pNTT)
				pNTTIn := pNTT.WithModIdx(vec.Range(0, inLen)...)

				crtOpIn := op.crtOp.WithModIdx(vec.Range(0, inLen)...)
				crtOpIn.InvNTTTo(pNTTIn, e.Value)

				emb.EmbedTo(eOut.Value, pNTTIn)
				eOut.Value.IsNTT = false
			} else {
				emb.EmbedTo(eOut.Value, e.Value)
				eOut.Value.IsNTT = false
			}

			if isNTT {
				op.FwdNTTTo(eOut, eOut)
			}
		}

	case inLen == outLen:
		eOut.CopyFrom(e)
		if isNTT && !eOut.Value.IsNTT {
			op.FwdNTTTo(eOut, eOut)
		} else if !isNTT && eOut.Value.IsNTT {
			op.InvNTTTo(eOut, eOut)
		}
	}

	eOut.auxLen = 0
}

// DivByAuxModulus returns round(e / AuxModulus).
func (op *PlainOperator) DivByAuxModulus(e *Element, isNTT bool) *Element {
	eOut := NewElement(e.Rank(), e.BaseModLen(), e.AuxModLen(), isNTT)
	op.DivByAuxModulusTo(eOut, e, isNTT)
	return eOut
}

// DivByAuxModulusTo computes round(e / AuxModulus) and stores the result in eOut.
func (op *PlainOperator) DivByAuxModulusTo(eOut, e *Element, isNTT bool) {
	if e.AuxModLen() == 0 {
		panic("input(s) must have auxiliary modulus")
	} else if eOut.BaseModLen() != e.BaseModLen() {
		panic("inconsistent output")
	}

	baseLen, auxLen := e.BaseModLen(), e.AuxModLen()
	paramAuxLen := len(op.Params.auxMod)

	opBase := op.crtOp.WithModIdx(vec.Range(paramAuxLen, paramAuxLen+baseLen)...)
	opAux := op.crtOp.WithModIdx(vec.Range(paramAuxLen-auxLen, paramAuxLen)...)

	baseMod := op.Params.baseMod[:baseLen]
	auxMod := op.Params.auxMod[paramAuxLen-auxLen : paramAuxLen]
	scaler := crt.NewScaler(baseMod, auxMod)

	auxInvBase := crt.NewScalarFrom(1, baseMod) // auxMod^{-1} mod baseMod
	baseInvAux := crt.NewScalarFrom(1, auxMod)  // baseMod^{-1} mod auxMod
	for i := range baseMod {
		for j := range auxMod {
			auxInvBase.Coeffs[i][0] = num.Mul(auxInvBase.Coeffs[i][0], num.Inv(auxMod[j].Value(), baseMod[i]), baseMod[i])
			baseInvAux.Coeffs[j][0] = num.Mul(baseInvAux.Coeffs[j][0], num.Inv(baseMod[i].Value(), auxMod[j]), auxMod[j])
		}
	}

	pDiv := op.pPool.Get()
	defer op.pPool.Put(pDiv)
	pDivBase := pDiv.WithModIdx(vec.Range(auxLen, auxLen+baseLen)...)
	pDivAux := pDiv.WithModIdx(vec.Range(0, auxLen)...)
	pDiv = pDiv.WithModIdx(vec.Range(0, baseLen+auxLen)...)

	pDiv.CopyFrom(e.Value)
	pDivBase.IsNTT = pDiv.IsNTT
	pDivAux.IsNTT = pDiv.IsNTT

	opBase.MulTo(pDivBase, pDivBase, auxInvBase)
	opAux.MulTo(pDivAux, pDivAux, baseInvAux)

	if pDivAux.IsNTT {
		opAux.InvNTTTo(pDivAux, pDivAux)
	}

	scaler.ScaleTo(eOut.Value, pDivAux)
	eOut.Value.IsNTT = false

	if isNTT {
		opBase.FwdNTTTo(eOut.Value, eOut.Value)
		if !pDivBase.IsNTT {
			opBase.FwdNTTTo(pDivBase, pDivBase)
		}
	} else {
		if pDivBase.IsNTT {
			opBase.InvNTTTo(pDivBase, pDivBase)
		}
	}
	opBase.AddTo(eOut.Value, eOut.Value, pDivBase)

	eOut.auxLen = 0
}

// Scale scales e to l-th modulus.
//
//   - When l < e.ModLen, it returns round(e / Modulus[l:e.ModLen]).
//   - When l > e.ModLen, it returns e * Modulus[l:e.ModLen].
//   - When l == e.ModLen, it returns a copy of e.
func (op *PlainOperator) Scale(e *Element, l int, isNTT bool) *Element {
	eOut := NewElement(e.Rank(), l, 0, isNTT)
	op.ScaleTo(eOut, e, l, isNTT)
	return eOut
}

// ScaleTo scales e to l-th modulus and stores the result in eOut.
//
//   - When l < e.ModLen, eOut = round(e / Modulus[l:e.ModLen]).
//   - When l > e.ModLen, eOut = e * Modulus[l:e.ModLen].
//   - When l == e.ModLen, eOut = e.
func (op *PlainOperator) ScaleTo(eOut, e *Element, l int, isNTT bool) {
	if e.AuxModLen() > 0 {
		panic("input(s) should not have auxiliary modulus")
	} else if eOut.BaseModLen() != l {
		panic("inconsistent output")
	}

	inLen := e.BaseModLen()
	outLen := l

	inMod := op.Params.baseMod[:inLen]
	outMod := op.Params.baseMod[:outLen]

	switch {
	case inLen > outLen:
		scMod := inMod[outLen:]

		auxLen := len(op.Params.auxMod)
		opOut := op.crtOp.WithModIdx(vec.Range(auxLen, auxLen+outLen)...)
		opScale := op.crtOp.WithModIdx(vec.Range(auxLen+outLen, auxLen+inLen)...)

		p := op.pPool.Get()
		defer op.pPool.Put(p)

		pIn := p.WithModIdx(vec.Range(0, inLen)...)
		pOut := p.WithModIdx(vec.Range(0, outLen)...)
		pScale := p.WithModIdx(vec.Range(outLen, inLen)...)

		scaler := crt.NewScaler(outMod, scMod)

		scInvOut := crt.NewScalarFrom(1, outMod) // scale^{-1} mod outMod
		outInvSc := crt.NewScalarFrom(1, scMod)  // outMod^{-1} mod scale
		for i := range outMod {
			for j := range scMod {
				scInvOut.Coeffs[i][0] = num.Mul(scInvOut.Coeffs[i][0], num.Inv(scMod[j].Value(), outMod[i]), outMod[i])
				outInvSc.Coeffs[j][0] = num.Mul(outInvSc.Coeffs[j][0], num.Inv(outMod[i].Value(), scMod[j]), scMod[j])
			}
		}

		pIn.CopyFrom(e.Value)
		pOut.IsNTT = pIn.IsNTT
		pScale.IsNTT = pIn.IsNTT

		opOut.MulTo(pOut, pOut, scInvOut)
		opScale.MulTo(pScale, pScale, outInvSc)
		if pScale.IsNTT {
			opScale.InvNTTTo(pScale, pScale)
		}

		scaler.ScaleTo(eOut.Value, pScale)
		eOut.Value.IsNTT = false

		if isNTT {
			opOut.FwdNTTTo(eOut.Value, eOut.Value)
			if !pOut.IsNTT {
				opOut.FwdNTTTo(pOut, pOut)
			}
		} else {
			if pOut.IsNTT {
				opOut.InvNTTTo(pOut, pOut)
			}
		}
		opOut.AddTo(eOut.Value, eOut.Value, pOut)

	case inLen < outLen:
		scMod := outMod[inLen:]

		auxLen := len(op.Params.auxMod)
		opOut := op.crtOp.WithModIdx(vec.Range(auxLen, auxLen+outLen)...)

		scale := crt.NewScalarFrom(1, inMod)
		for i := range inMod {
			for j := range scMod {
				scale.Coeffs[i][0] = num.Mul(scale.Coeffs[i][0], scMod[j].Value(), inMod[i])
			}
		}

		eOut.Value.IsNTT = e.Value.IsNTT

		for i := 0; i < inLen; i++ {
			vec.MulScalarTo(eOut.Value.Coeffs[i], e.Value.Coeffs[i], scale.Coeffs[i][0], outMod[i])
		}
		for i := inLen; i < outLen; i++ {
			clear(eOut.Value.Coeffs[i])
		}

		if isNTT {
			if !eOut.Value.IsNTT {
				opOut.FwdNTTTo(eOut.Value, eOut.Value)
			}
		} else {
			if eOut.Value.IsNTT {
				opOut.InvNTTTo(eOut.Value, eOut.Value)
			}
		}

	case inLen == outLen:
		eOut.CopyFrom(e)

		if isNTT {
			if !eOut.Value.IsNTT {
				op.FwdNTTTo(eOut, eOut)
			}
		} else {
			if eOut.Value.IsNTT {
				op.InvNTTTo(eOut, eOut)
			}
		}
	}

	eOut.auxLen = 0
}

// Operator evaluates operations over [*Ciphertext].
type Operator struct {
	Params Parameters

	plainOp *PlainOperator
	dcmp    Decomposer

	pPool    *pool.Pool[*Element]
	ctPool   *pool.Pool[*Ciphertext]
	dcmpPool *pool.Pool[*Vector]
}

// NewOperator creates a new [Operator].
func NewOperator(params Parameters) *Operator {
	return &Operator{
		Params: params,

		plainOp: NewPlainOperator(params),
		dcmp:    NewDecomposer(params),

		pPool: pool.NewPool(func() *Element {
			return NewElement(params.Rank(), len(params.baseMod), len(params.auxMod), true)
		}),
		ctPool: pool.NewPool(func() *Ciphertext {
			return NewCiphertext(params, params.HasAuxModulus(), true)
		}),
		dcmpPool: pool.NewPool(func() *Vector {
			return NewVector(params, params.GadgetLen(), params.HasAuxModulus(), false)
		}),
	}
}

// PlainOperator returns the underlying [PlainOperator].
func (op *Operator) PlainOperator() *PlainOperator {
	return op.plainOp
}

// Decomposer returns the underlying [Decomposer].
func (op *Operator) Decomposer() Decomposer {
	return op.dcmp
}

// AddElement returns ct + pt.
func (op *Operator) AddElement(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), ct.IsNTT())
	op.AddElementTo(ctOut, ct, pt)
	return ctOut
}

// AddElementTo computes ctOut = ct + pt.
func (op *Operator) AddElementTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.AddTo(ctOut.Body, ct.Body, pt)
	ctOut.Mask.CopyFrom(ct.Mask)
}

// SubElement returns ct - pt.
func (op *Operator) SubElement(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), ct.IsNTT())
	op.SubElementTo(ctOut, ct, pt)
	return ctOut
}

// SubElementTo computes ctOut = ct - pt.
func (op *Operator) SubElementTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.SubTo(ctOut.Body, ct.Body, pt)
	ctOut.Mask.CopyFrom(ct.Mask)
}

// MulElement returns ct * pt.
func (op *Operator) MulElement(ct *Ciphertext, pt *Element) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), ct.IsNTT())
	op.MulElementTo(ctOut, ct, pt)
	return ctOut
}

// MulElementTo computes ctOut = ct * pt.
func (op *Operator) MulElementTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.MulTo(ctOut.Body, ct.Body, pt)
	op.plainOp.MulTo(ctOut.Mask, ct.Mask, pt)
}

// MulAddElementTo computes ctOut += ct * pt.
func (op *Operator) MulAddElementTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.MulAddTo(ctOut.Body, ct.Body, pt)
	op.plainOp.MulAddTo(ctOut.Mask, ct.Mask, pt)
}

// MulSubElementTo computes ctOut -= ct * pt.
func (op *Operator) MulSubElementTo(ctOut, ct *Ciphertext, pt *Element) {
	op.plainOp.MulSubTo(ctOut.Body, ct.Body, pt)
	op.plainOp.MulSubTo(ctOut.Mask, ct.Mask, pt)
}

// Add returns ct0 + ct1.
func (op *Operator) Add(ct0, ct1 *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.BaseModLen(), ct0.AuxModLen(), ct0.IsNTT())
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
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.BaseModLen(), ct0.AuxModLen(), ct0.IsNTT())
	op.SubTo(ctOut, ct0, ct1)
	return ctOut
}

// SubTo computes ctOut = ct0 - ct1.
func (op *Operator) SubTo(ctOut, ct0, ct1 *Ciphertext) {
	op.plainOp.SubTo(ctOut.Body, ct0.Body, ct1.Body)
	op.plainOp.SubTo(ctOut.Mask, ct0.Mask, ct1.Mask)
}

func (op *Operator) Neg(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), ct.IsNTT())
	op.NegTo(ctOut, ct)
	return ctOut
}

// NegTo computes ctOut = -ct.
func (op *Operator) NegTo(ctOut, ct *Ciphertext) {
	op.plainOp.NegTo(ctOut.Body, ct.Body)
	op.plainOp.NegTo(ctOut.Mask, ct.Mask)
}

// TODO: FwdNTT/InvNTT of ciphertexts are very common, especially with isNTT flags.
// We should create a seperate structure for this. (something like `rlwe.Transformer`.)

// FwdNTT returns FwdNTT(ct).
func (op *Operator) FwdNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), ct.IsNTT())
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
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), ct.IsNTT())
	op.InvNTTTo(ctOut, ct)
	return ctOut
}

// InvNTTTo computes ctOut = InvNTT(ct).
func (op *Operator) InvNTTTo(ctOut, ct *Ciphertext) {
	op.plainOp.InvNTTTo(ctOut.Body, ct.Body)
	op.plainOp.InvNTTTo(ctOut.Mask, ct.Mask)
}

// ModRaise raises ct to the modulus with length l.
func (op *Operator) ModRaise(ct *Ciphertext, l int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), l, 0, ct.IsNTT())
	op.ModRaiseTo(ctOut, ct, isNTT)
	return ctOut
}

// ModRaiseTo raises ct to the modulus with length l and stores the result in ctOut.
func (op *Operator) ModRaiseTo(ctOut, ct *Ciphertext, isNTT bool) {
	op.plainOp.ModRaiseTo(ctOut.Body, ct.Body, isNTT)
	op.plainOp.ModRaiseTo(ctOut.Mask, ct.Mask, isNTT)
}

// DivByAuxModulus returns round(ct / AuxModulus).
func (op *Operator) DivByAuxModulus(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), 0, ct.IsNTT())
	op.DivByAuxModulusTo(ctOut, ct, isNTT)
	return ctOut
}

// DivByAuxModulusTo computes ctOut = round(ct / AuxModulus).
func (op *Operator) DivByAuxModulusTo(ctOut, ct *Ciphertext, isNTT bool) {
	op.plainOp.DivByAuxModulusTo(ctOut.Body, ct.Body, isNTT)
	op.plainOp.DivByAuxModulusTo(ctOut.Mask, ct.Mask, isNTT)
}

// Scale scales ct to l-th modulus.
//
//   - When l < ct.ModLen, it returns round(ct / Modulus[l:ct.ModLen]).
//   - When l > ct.ModLen, it returns ct * Modulus[l:ct.ModLen].
//   - When l == ct.ModLen, it returns a copy of ct.
func (op *Operator) Scale(ct *Ciphertext, l int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), l, 0, ct.IsNTT())
	op.ScaleTo(ctOut, ct, l, isNTT)
	return ctOut
}

// ScaleTo scales ct to l-th modulus and writes the result to ctOut.
//
//   - When l < ct.ModLen, ctOut = round(ct / Modulus[l:ct.ModLen]).
//   - When l > ct.ModLen, ctOut = ct * Modulus[l:ct.ModLen].
//   - When l == ct.ModLen, ctOut = ct.
func (op *Operator) ScaleTo(ctOut, ct *Ciphertext, l int, isNTT bool) {
	op.plainOp.ScaleTo(ctOut.Body, ct.Body, l, isNTT)
	op.plainOp.ScaleTo(ctOut.Mask, ct.Mask, l, isNTT)
}

// tensorTo tensors two ciphertexts in the ambient modulus into a vector.
//
// Input must be in NTT form, and output is in NTT form.
func (op *Operator) TensorTo(vOut *Vector, ct0, ct1 *Ciphertext) {
	if vOut.BaseModLen() != ct0.BaseModLen() || vOut.BaseModLen() != ct1.BaseModLen() || vOut.AuxModLen() != 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.plainOp
	pOp.MulTo(vOut.Value[0], ct0.Body, ct1.Body)
	pOp.MulTo(vOut.Value[1], ct0.Body, ct1.Mask)
	pOp.MulAddTo(vOut.Value[1], ct0.Mask, ct1.Body)
	pOp.MulTo(vOut.Value[2], ct0.Mask, ct1.Mask)
}

// HoistedGadgetProdLazy returns ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) HoistedGadgetProdLazy(pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(op.Params.Rank(), pDcmp.BaseModLen(), 0, true)
	op.HoistedGadgetProdLazyTo(ctOut, pDcmp, ctGadEnc, isNTT)
	return ctOut
}

// HoistedGadgetProdLazyTo computes ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) HoistedGadgetProdLazyTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) {
	if ctOut.AuxModLen() != pDcmp.AuxModLen() || ctOut.AuxModLen() != op.dcmp.AuxModLen(ctOut.BaseModLen()) {
		panic("inconsistent input(s)")
	}

	baseLen, auxLen := pDcmp.BaseModLen(), pDcmp.AuxModLen()
	if ctOut.BaseModLen() != baseLen || pDcmp.Len() != op.dcmp.DecomposeLen(baseLen) {
		panic("inconsistent input(s)")
	}

	ctBuf := op.ctPool.Get()
	defer op.ctPool.Put(ctBuf)

	ctBuf = ctBuf.WithModLen(baseLen, auxLen)
	ctBuf.Clear()

	for i := 0; i < pDcmp.Len(); i++ {
		op.plainOp.MulAddTo(ctBuf.Body, pDcmp.Value[i], ctGadEnc.Value[i].Body.WithModLen(baseLen, auxLen))
		op.plainOp.MulAddTo(ctBuf.Mask, pDcmp.Value[i], ctGadEnc.Value[i].Mask.WithModLen(baseLen, auxLen))
	}

	if !isNTT {
		op.InvNTTTo(ctBuf, ctBuf)
	}

	ctOut.CopyFrom(ctBuf)
}

// HoistedGadgetProd returns ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
func (op *Operator) HoistedGadgetProd(pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	baseLen, auxLen := pDcmp.BaseModLen(), pDcmp.AuxModLen()
	ctOut := NewCiphertextCustom(op.Params.Rank(), baseLen, auxLen, true)
	op.HoistedGadgetProdTo(ctOut, pDcmp, ctGadEnc, isNTT)
	return ctOut
}

// HoistedGadgetProdTo computes ctOut = p * ctGadEnc, where the decomposition of p is precomputed.
func (op *Operator) HoistedGadgetProdTo(ctOut *Ciphertext, pDcmp *Vector, ctGadEnc *GadgetEncryption, isNTT bool) {
	ctOutAux := op.ctPool.Get()
	defer op.ctPool.Put(ctOutAux)
	ctOutAux = ctOutAux.WithModLen(pDcmp.BaseModLen(), pDcmp.AuxModLen())

	op.HoistedGadgetProdLazyTo(ctOutAux, pDcmp, ctGadEnc, true)

	if op.Params.HasAuxModulus() {
		op.DivByAuxModulusTo(ctOut, ctOutAux, isNTT)
	} else {
		ctOut.CopyFrom(ctOutAux)
		if !isNTT {
			op.InvNTTTo(ctOut, ctOut)
		}
	}
}

// GadgetProdLazy returns ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProdLazy(p *Element, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	baseLen, auxLen := p.BaseModLen(), op.dcmp.AuxModLen(p.BaseModLen())
	ctOut := NewCiphertextCustom(op.Params.Rank(), baseLen, auxLen, true)
	op.GadgetProdLazyTo(ctOut, p, ctGadEnc, isNTT)
	return ctOut
}

// GadgetProdLazyTo computes ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProdLazyTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, isNTT bool) {
	if p.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	} else if ctOut.BaseModLen() != p.BaseModLen() || ctOut.AuxModLen() != op.dcmp.AuxModLen(p.BaseModLen()) {
		panic("inconsistent input(s)")
	}

	baseLen, auxLen := p.BaseModLen(), op.dcmp.AuxModLen(p.BaseModLen())

	pInvNTT := op.pPool.Get()
	defer op.pPool.Put(pInvNTT)
	pInvNTT = pInvNTT.WithModLen(baseLen, 0)
	if p.Value.IsNTT {
		op.plainOp.InvNTTTo(pInvNTT, p)
	} else {
		pInvNTT.CopyFrom(p)
	}

	pDcmp := op.dcmpPool.Get()
	defer op.dcmpPool.Put(pDcmp)
	pDcmp = pDcmp.Slice(vec.Range(0, op.dcmp.DecomposeLen(baseLen))...).WithModLen(baseLen, auxLen)
	op.dcmp.DecomposeTo(pDcmp, pInvNTT)

	for i := range pDcmp.Value {
		op.plainOp.FwdNTTTo(pDcmp.Value[i], pDcmp.Value[i])
	}

	op.HoistedGadgetProdLazyTo(ctOut, pDcmp, ctGadEnc, isNTT)
}

// GadgetProd returns ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProd(p *Element, ctGadEnc *GadgetEncryption, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(op.Params.Rank(), p.BaseModLen(), 0, isNTT)
	op.GadgetProdTo(ctOut, p, ctGadEnc, isNTT)
	return ctOut
}

// GadgetProdTo computes ctOut = p * ctGadEnc.
// The modulus of the output includes the auxillary modulus if present.
func (op *Operator) GadgetProdTo(ctOut *Ciphertext, p *Element, ctGadEnc *GadgetEncryption, isNTT bool) {
	if p.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	baseLen, auxLen := p.BaseModLen(), op.dcmp.AuxModLen(p.BaseModLen())

	pInvNTT := op.pPool.Get()
	defer op.pPool.Put(pInvNTT)
	pInvNTT = pInvNTT.WithModLen(baseLen, 0)
	if p.Value.IsNTT {
		op.plainOp.InvNTTTo(pInvNTT, p)
	} else {
		pInvNTT.CopyFrom(p)
	}

	pDcmp := op.dcmpPool.Get()
	defer op.dcmpPool.Put(pDcmp)
	pDcmp = pDcmp.Slice(vec.Range(0, op.dcmp.DecomposeLen(baseLen))...).WithModLen(baseLen, auxLen)
	op.dcmp.DecomposeTo(pDcmp, pInvNTT)

	for i := range pDcmp.Value {
		op.plainOp.FwdNTTTo(pDcmp.Value[i], pDcmp.Value[i])
	}

	op.HoistedGadgetProdTo(ctOut, pDcmp, ctGadEnc, isNTT)
}

// Relin performs a relinearisation and returns the result.
func (op *Operator) Relin(cIn *Vector, rlk *RelinKey, isNTT bool) *Ciphertext {
	cOut := NewCiphertextCustom(op.Params.Rank(), cIn.BaseModLen(), 0, false)
	op.RelinTo(cOut, cIn, rlk, isNTT)
	return cOut
}

// RelinTo performs a relinearisation and stores the result in cOut.
// Revise to deal with the case where the input and output share the same elements.
func (op *Operator) RelinTo(cOut *Ciphertext, cIn *Vector, rlk *RelinKey, isNTT bool) {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	} else if cIn.Len() != 3 {
		panic("Ciphertext must have degree 3.")
	}

	baseLen := cIn.BaseModLen()
	pOp := op.plainOp

	pNTT := op.pPool.Get()
	defer op.pPool.Put(pNTT)
	pNTT = pNTT.WithModLen(baseLen, 0)

	op.GadgetProdTo(cOut, cIn.Value[2], (*GadgetEncryption)(rlk), isNTT)

	if cIn.Value[0].IsNTT() && !isNTT {
		pOp.InvNTTTo(pNTT, cIn.Value[0])
	} else if !cIn.Value[0].IsNTT() && isNTT {
		pOp.FwdNTTTo(pNTT, cIn.Value[0])
	} else {
		pNTT.CopyFrom(cIn.Value[0])
	}
	pOp.AddTo(cOut.Body, cOut.Body, pNTT)

	if cIn.Value[1].IsNTT() && !isNTT {
		pOp.InvNTTTo(pNTT, cIn.Value[1])
	} else if !cIn.Value[1].IsNTT() && isNTT {
		pOp.FwdNTTTo(pNTT, cIn.Value[1])
	} else {
		pNTT.CopyFrom(cIn.Value[1])
	}
	pOp.AddTo(cOut.Mask, cOut.Mask, pNTT)
}

// KeySwitch performs a key switch and returns the result.
func (op *Operator) KeySwitch(cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) *Ciphertext {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	cOut := NewCiphertextCustom(op.Params.Rank(), cIn.BaseModLen(), 0, false)
	op.KeySwitchTo(cOut, cIn, ksk, isNTT)
	return cOut
}

// KeySwitchTo performs a key switch and stores the result in cOut.
func (op *Operator) KeySwitchTo(cOut *Ciphertext, cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	pOp := op.plainOp

	baseLen := cIn.BaseModLen()
	auxLen := op.dcmp.AuxModLen(baseLen)
	dcmpLen := op.dcmp.DecomposeLen(baseLen)

	pNTT := op.pPool.Get()
	defer op.pPool.Put(pNTT)
	pNTT = pNTT.WithModLen(baseLen, 0)

	if cIn.Mask.IsNTT() {
		pOp.InvNTTTo(pNTT, cIn.Mask)
	} else {
		pNTT.CopyFrom(cIn.Mask)
	}

	pDcmp := op.dcmpPool.Get()
	defer op.dcmpPool.Put(pDcmp)
	pDcmp = pDcmp.Slice(vec.Range(0, dcmpLen)...).WithModLen(baseLen, auxLen)

	op.dcmp.DecomposeTo(pDcmp, pNTT)
	for i := range pDcmp.Value {
		pOp.FwdNTTTo(pDcmp.Value[i], pDcmp.Value[i])
	}

	op.HoistedKeySwitchTo(cOut, pDcmp, cIn, ksk, isNTT)
}

// HoistedKeySwitch performs a hoisted key switch and returns the result.
func (op *Operator) HoistedKeySwitch(decmp *Vector, cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) *Ciphertext {
	cOut := NewCiphertextCustom(op.Params.Rank(), cIn.BaseModLen(), 0, false)
	op.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)
	return cOut
}

// HoistedKeySwitchTo performs a hoisted key switch and stores the result in cOut.
func (op *Operator) HoistedKeySwitchTo(cOut *Ciphertext, decmp *Vector, cIn *Ciphertext, ksk *KeySwitchKey, isNTT bool) {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	pOp := op.plainOp

	pNTT := op.pPool.Get()
	defer op.pPool.Put(pNTT)
	pNTT = pNTT.WithModLen(cIn.BaseModLen(), 0)

	if cIn.Body.IsNTT() && !isNTT {
		pOp.InvNTTTo(pNTT, cIn.Body)
	} else if !cIn.Body.IsNTT() && isNTT {
		pOp.FwdNTTTo(pNTT, cIn.Body)
	} else {
		pNTT.CopyFrom(cIn.Body)
	}

	kskGad := (*GadgetEncryption)(ksk)
	op.HoistedGadgetProdTo(cOut, decmp, kskGad, isNTT)
	pOp.AddTo(cOut.Body, cOut.Body, pNTT)
}

// Aut performs an automorphism and returns the result.
func (op *Operator) Aut(cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) *Ciphertext {
	cOut := NewCiphertextCustom(op.Params.Rank(), cIn.BaseModLen(), 0, false)
	op.AutTo(cOut, cIn, atk, isNTT)
	return cOut
}

// AutTo performs an automorphism and stores the result in cOut.
func (op *Operator) AutTo(cOut *Ciphertext, cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) {
	ksk := &KeySwitchKey{Value: atk.Value}
	op.KeySwitchTo(cOut, cIn, ksk, isNTT)

	pOp := op.plainOp
	pOp.AutTo(cOut.Body, cOut.Body, atk.Idx)
	pOp.AutTo(cOut.Mask, cOut.Mask, atk.Idx)
}

// HoistedAut performs a hoisted automorphism and returns the result.
func (op *Operator) HoistedAut(decmp *Vector, cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) *Ciphertext {
	cOut := NewCiphertextCustom(op.Params.Rank(), cIn.BaseModLen(), 0, false)
	op.HoistedAutTo(cOut, decmp, cIn, atk, isNTT)
	return cOut
}

// HoistedAutTo performs a hoisted automorphism and stores the result in cOut.
func (op *Operator) HoistedAutTo(cOut *Ciphertext, decmp *Vector, cIn *Ciphertext, atk *AutomorphismKey, isNTT bool) {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	ksk := &KeySwitchKey{Value: atk.Value}
	op.HoistedKeySwitchTo(cOut, decmp, cIn, ksk, isNTT)

	pOp := op.plainOp
	pOp.AutTo(cOut.Body, cOut.Body, atk.Idx)
	pOp.AutTo(cOut.Mask, cOut.Mask, atk.Idx)
}

// ExtProd performs an external product and returns the result.
func (op *Operator) ExtProd(cIn *Ciphertext, gsw *RGSW, isNTT bool) *Ciphertext {
	cOut := NewCiphertextCustom(op.Params.Rank(), cIn.BaseModLen(), 0, false)
	op.ExtProdTo(cOut, cIn, gsw, isNTT)
	return cOut
}

// ExtProdTo performs an external product and stores the result in cOut.
func (op *Operator) ExtProdTo(cOut *Ciphertext, cIn *Ciphertext, gsw *RGSW, isNTT bool) {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	ctExt := op.ctPool.Get()
	defer op.ctPool.Put(ctExt)
	ctExt = ctExt.WithModLen(cIn.BaseModLen(), 0)

	op.GadgetProdTo(ctExt, cIn.Body, gsw.Body, isNTT)
	op.GadgetProdTo(cOut, cIn.Mask, gsw.Mask, isNTT)

	op.AddTo(cOut, cOut, ctExt)
}

// ExtProdLazyTo performs an external product and stores the result in cOut.
func (op *Operator) ExtProdLazyTo(cOut *Ciphertext, cIn *Ciphertext, gsw *RGSW, isNTT bool) {
	if cIn.AuxModLen() > 0 {
		panic("Ciphertext must not have auxiliary modulus.")
	}

	ctExt := op.ctPool.Get()
	defer op.ctPool.Put(ctExt)
	ctExt = ctExt.WithModLen(cIn.BaseModLen(), op.dcmp.AuxModLen(cIn.BaseModLen()))

	op.GadgetProdLazyTo(ctExt, cIn.Body, gsw.Body, isNTT)
	op.GadgetProdLazyTo(cOut, cIn.Mask, gsw.Mask, isNTT)

	op.AddTo(cOut, cOut, ctExt)
}

// HoistedExtProdTo performs a hoisted external product and stores the result in cOut.
func (op *Operator) HoistedExtProdTo(cOut *Ciphertext, decmpBody *Vector, decmpMask *Vector, gsw *RGSW, isNTT bool) {
	ctExt := op.ctPool.Get()
	defer op.ctPool.Put(ctExt)
	ctExt = ctExt.WithModLen(decmpBody.BaseModLen(), 0)

	op.HoistedGadgetProdTo(ctExt, decmpBody, gsw.Body, isNTT)
	op.HoistedGadgetProdTo(cOut, decmpMask, gsw.Mask, isNTT)

	op.AddTo(cOut, cOut, ctExt)
}

// HoistedExtProdLazyTo performs a hoisted external product and stores the result in cOut.
func (op *Operator) HoistedExtProdLazyTo(cOut *Ciphertext, decmpBody *Vector, decmpMask *Vector, gsw *RGSW, isNTT bool) {
	ctExt := op.ctPool.Get()
	defer op.ctPool.Put(ctExt)
	ctExt = ctExt.WithModLen(decmpBody.BaseModLen(), op.dcmp.AuxModLen(decmpBody.BaseModLen()))

	op.HoistedGadgetProdLazyTo(ctExt, decmpBody, gsw.Body, isNTT)
	op.HoistedGadgetProdLazyTo(cOut, decmpMask, gsw.Mask, isNTT)

	op.AddTo(cOut, cOut, ctExt)
}
