package ckks

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates operations over [*Ciphertext] and [*Plaintext].
type Operator struct {
	params rlwe.Parameters
	scFac  float64

	rlweOp *rlwe.Operator
	ecd    *Encoder

	fPool   *pool.Pool[*big.Float]
	intPool *pool.Pool[*big.Int]

	ePool  *rlwe.ElementPool
	ctPool *pool.Pool[*Ciphertext]
	vPool  *pool.Pool[*rlwe.Vector]
}

func NewOperator(params rlwe.Parameters, scFac float64) *Operator {
	return &Operator{
		params: params,
		scFac:  scFac,

		rlweOp: rlwe.NewOperator(params),
		ecd:    NewEncoder(params),

		fPool: pool.NewPool(func() *big.Float {
			return big.NewFloat(0).SetPrec(52 + 52 + 64)
		}),
		intPool: pool.NewPool(func() *big.Int {
			return new(big.Int)
		}),
		ePool: rlwe.NewElementPool(params, true, true),
		ctPool: pool.NewPool(func() *Ciphertext {
			return NewCiphertext(params, true)
		}),
		vPool: pool.NewPool(func() *rlwe.Vector {
			return rlwe.NewVector(params, 3, false, true)
		}),
	}
}

// Parameters returns the parameters.
func (op *Operator) Parameters() rlwe.Parameters {
	return op.params
}

// Encoder returns the encoder.
func (op *Operator) Encoder() *Encoder {
	return op.ecd
}

// Rescale rescales the ciphertext to the target modulus.
func (op *Operator) Rescale(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.RescaleTo(ctOut, ct, isNTT)
	return ctOut
}

// RescaleTo rescales the ciphertext to the target modulus.
func (op *Operator) RescaleTo(ctOut, ct *Ciphertext, isNTT bool) {
	resScFac := op.fPool.Get()
	tmpFloat := op.fPool.Get()
	divMod := op.fPool.Get()
	scFac := op.fPool.Get()
	defer op.fPool.Put(resScFac)
	defer op.fPool.Put(tmpFloat)
	defer op.fPool.Put(divMod)
	defer op.fPool.Put(scFac)

	resScFac.SetFloat64(ct.scFac)
	scFac.SetFloat64(op.scFac)
	tarLen := ct.ModLen()
	for tarLen > 0 {
		divMod.SetUint64(op.params.BaseModulus()[tarLen-1].Value())
		if tmpFloat.Quo(resScFac, divMod).Cmp(scFac) > 0 {
			tarLen--
			resScFac.Quo(resScFac, divMod)
		} else {
			break
		}
	}
	if tarLen == 0 {
		panic("cannot rescale ciphertext")
	}

	buf := op.ctPool.Get()
	defer op.ctPool.Put(buf)
	buf = buf.WithModLen(tarLen)

	op.rlweOp.ScaleTo(buf.Value, ct.Value, tarLen, isNTT)
	ctOut.Value.Resize(tarLen, 0)
	ctOut.Value.CopyFrom(buf.Value)
	ctOut.scFac, _ = resScFac.Float64()
}

// FwdNTT returns FwdNTT(ct).
func (op *Operator) FwdNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.FwdNTTTo(ctOut, ct)
	return ctOut
}

// FwdNTTTo computes ctOut = FwdNTT(ct).
func (op *Operator) FwdNTTTo(ctOut, ct *Ciphertext) {
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
	op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
}

// InvNTT returns InvNTT(ct).
func (op *Operator) InvNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.InvNTTTo(ctOut, ct)
	return ctOut
}

// InvNTTTo computes ctOut = InvNTT(ct).
func (op *Operator) InvNTTTo(ctOut, ct *Ciphertext) {
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
	op.rlweOp.InvNTTTo(ctOut.Value, ct.Value)
}

// Neg computes ctOut = -ct.
func (op *Operator) Neg(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.NegTo(ctOut, ct, isNTT)
	return ctOut
}

// NegTo computes ctOut = -ct.
func (op *Operator) NegTo(ctOut, ct *Ciphertext, isNTT bool) {
	ctOut.Resize(ct.ModLen())
	op.rlweOp.NegTo(ctOut.Value, ct.Value)
	if isNTT && !ctOut.IsNTT() {
		op.rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
	} else if !isNTT && ctOut.IsNTT() {
		op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
	}
	ctOut.scFac = ct.scFac
}

// Add computes ctOut = ct0 + ct1.
func (op *Operator) Add(ct0, ct1 *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.AddTo(ctOut, ct0, ct1, isNTT)
	return ctOut
}

// AddTo computes ctOut = ct0 + ct1.
func (op *Operator) AddTo(ctOut, ct0, ct1 *Ciphertext, isNTT bool) {
	// Ensure ct0 has smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.scFac < ct1.scFac) {
		ct0, ct1 = ct1, ct0
	}
	tarLen := ct0.ModLen()
	curLen := ct1.ModLen()

	// Find 'diff', which needs to be multiplied to ct1 to make the scaling factors equal.
	diffFloat := op.fPool.Get()
	tmpFloat := op.fPool.Get()
	diffRound := op.intPool.Get()
	tmpInt := op.intPool.Get()
	defer op.fPool.Put(diffFloat)
	defer op.fPool.Put(tmpFloat)
	defer op.intPool.Put(diffRound)
	defer op.intPool.Put(tmpInt)

	diff := op.ePool.Get(crt.TypeScalar)
	defer op.ePool.Put(diff)
	diff = diff.WithModLen(curLen, 0)

	diffFloat.SetFloat64(ct0.scFac)
	diffFloat.Quo(diffFloat, tmpFloat.SetFloat64(ct1.scFac))
	for i := tarLen; i < curLen; i++ {
		tmpFloat.SetUint64(op.params.BaseModulus()[i].Value())
		diffFloat.Mul(diffFloat, tmpFloat)
	}
	diffFloat.Add(diffFloat, tmpFloat.SetFloat64(0.5))
	diffFloat.Int(diffRound)

	for i := 0; i < curLen; i++ {
		tmpInt.SetUint64(op.params.BaseModulus()[i].Value())
		tmpInt.Rem(diffRound, tmpInt)
		diff.Value.Coeffs[i][0] = tmpInt.Uint64()
	}

	// multiply diff to ct1 and scale to the target modulus.
	c0 := op.ctPool.Get()
	c1 := op.ctPool.Get()
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)
	c0 = c0.WithModLen(tarLen)
	c1 = c1.WithModLen(curLen)
	c1Tar := c1.WithModLen(tarLen)

	if isNTT && !ct0.IsNTT() {
		op.rlweOp.FwdNTTTo(c0.Value, ct0.Value)
	} else if !isNTT && ct0.IsNTT() {
		op.rlweOp.InvNTTTo(c0.Value, ct0.Value)
	} else {
		c0.CopyFrom(ct0)
	}

	op.rlweOp.MulElementTo(c1.Value, ct1.Value, diff)
	op.rlweOp.ScaleTo(c1Tar.Value, c1.Value, tarLen, isNTT)

	ctOut.Value.Resize(tarLen, 0)
	op.rlweOp.AddTo(ctOut.Value, c0.Value, c1Tar.Value)
	ctOut.scFac = ct0.scFac
}

// AddPlain computes ctOut = ct + pt.
func (op *Operator) AddPlain(ct *Ciphertext, pt *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.AddPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// AddPlainTo computes ctOut = ct + pt.
func (op *Operator) AddPlainTo(ctOut, ct *Ciphertext, pt *Plaintext, isNTT bool) {
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()

	var e *rlwe.Element
	switch pt.Rank() {
	case 1:
		e = op.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		e = op.ePool.Get(crt.TypePoly)
	}
	defer op.ePool.Put(e)
	e = e.WithModLen(ct.ModLen(), 0)

	op.ecd.EncodeTo(e, pt, ct.scFac, false)
	ctOut.Resize(ct.ModLen())
	if isNTT && !ct.IsNTT() {
		rlweOp.AddElementTo(ctOut.Value, ct.Value, e)
		rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
	} else if !isNTT && ct.IsNTT() {
		rlweOp.InvNTTTo(ctOut.Value, ct.Value)
		rlweOp.AddElementTo(ctOut.Value, ctOut.Value, e)
	} else {
		if isNTT {
			pOp.FwdNTTTo(e, e)
		}
		rlweOp.AddElementTo(ctOut.Value, ct.Value, e)
	}

	ctOut.scFac = ct.scFac
}

// AddElement computes ctOut = ct + e.
func (op *Operator) AddElement(ct *Ciphertext, e *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.AddElementTo(ctOut, ct, e, isNTT)
	return ctOut
}

// AddElementTo computes ctOut = ct + e.
func (op *Operator) AddElementTo(ctOut, ct *Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.ModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.rlweOp.PlainOperator()

	ctOut.Resize(ct.ModLen())
	if ct.IsNTT() == e.IsNTT() { // First add, then perform NTT/iNTT needed.
		op.rlweOp.AddElementTo(ctOut.Value, ct.Value, e)
		if ct.IsNTT() && !isNTT {
			op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		} else if !ct.IsNTT() && isNTT {
			op.rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		}
	} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then add.
		eNTT := op.ePool.Get(e.Type())
		defer op.ePool.Put(eNTT)
		eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

		if e.IsNTT() {
			pOp.InvNTTTo(eNTT, e)
		} else {
			pOp.FwdNTTTo(eNTT, e)
		}
		op.rlweOp.AddElementTo(ctOut.Value, ct.Value, eNTT)
	} else { // First perform NTT/iNTT on ct, then add.
		if ct.IsNTT() {
			op.rlweOp.InvNTTTo(ctOut.Value, ct.Value)
		} else {
			op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
		}
		op.rlweOp.AddElementTo(ctOut.Value, ctOut.Value, e)
	}

	ctOut.scFac = ct.scFac
}

// Sub computes ctOut = ct0 - ct1.
func (op *Operator) Sub(ct0, ct1 *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.SubTo(ctOut, ct0, ct1, isNTT)
	return ctOut
}

// SubTo computes ctOut = ct0 - ct1.
func (op *Operator) SubTo(ctOut, ct0, ct1 *Ciphertext, isNTT bool) {
	// Ensure ct0 has smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.scFac < ct1.scFac) {
		ct0, ct1 = ct1, ct0
	}
	tarLen := ct0.ModLen()
	curLen := ct1.ModLen()

	// Find 'diff', which needs to be multiplied to ct1 to make the scaling factors equal.
	diffFloat := op.fPool.Get()
	tmpFloat := op.fPool.Get()
	diffRound := op.intPool.Get()
	tmpInt := op.intPool.Get()
	defer op.fPool.Put(diffFloat)
	defer op.fPool.Put(tmpFloat)
	defer op.intPool.Put(diffRound)
	defer op.intPool.Put(tmpInt)

	diff := op.ePool.Get(crt.TypeScalar)
	defer op.ePool.Put(diff)
	diff = diff.WithModLen(curLen, 0)

	diffFloat.SetFloat64(ct0.scFac)
	diffFloat.Quo(diffFloat, tmpFloat.SetFloat64(ct1.scFac))
	for i := tarLen; i < curLen; i++ {
		tmpFloat.SetUint64(op.params.BaseModulus()[i].Value())
		diffFloat.Mul(diffFloat, tmpFloat)
	}
	diffFloat.Add(diffFloat, tmpFloat.SetFloat64(0.5))
	diffFloat.Int(diffRound)

	for i := 0; i < curLen; i++ {
		tmpInt.SetUint64(op.params.BaseModulus()[i].Value())
		tmpInt.Rem(diffRound, tmpInt)
		diff.Value.Coeffs[i][0] = tmpInt.Uint64()
	}

	// multiply diff to ct1 and scale to the target modulus.
	c0 := op.ctPool.Get()
	c1 := op.ctPool.Get()
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)
	c0 = c0.WithModLen(tarLen)
	c1 = c1.WithModLen(curLen)
	c1Tar := c1.WithModLen(tarLen)

	if isNTT && !ct0.IsNTT() {
		op.rlweOp.FwdNTTTo(c0.Value, ct0.Value)
	} else if !isNTT && ct0.IsNTT() {
		op.rlweOp.InvNTTTo(c0.Value, ct0.Value)
	} else {
		c0.CopyFrom(ct0)
	}

	op.rlweOp.MulElementTo(c1.Value, ct1.Value, diff)
	op.rlweOp.ScaleTo(c1Tar.Value, c1.Value, tarLen, isNTT)

	ctOut.Value.Resize(tarLen, 0)
	op.rlweOp.SubTo(ctOut.Value, c0.Value, c1Tar.Value)
	ctOut.scFac = ct0.scFac
}

// SubPlain computes ctOut = ct - pt.
func (op *Operator) SubPlain(ct *Ciphertext, pt *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.SubPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// SubPlainTo computes ctOut = ct - pt.
func (op *Operator) SubPlainTo(ctOut, ct *Ciphertext, pt *Plaintext, isNTT bool) {
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()

	var e *rlwe.Element
	switch pt.Rank() {
	case 1:
		e = op.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		e = op.ePool.Get(crt.TypePoly)
	}
	defer op.ePool.Put(e)
	e = e.WithModLen(ct.ModLen(), 0)

	op.ecd.EncodeTo(e, pt, ct.scFac, false)
	ctOut.Resize(ct.ModLen())
	if isNTT && !ct.IsNTT() {
		rlweOp.SubElementTo(ctOut.Value, ct.Value, e)
		rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
	} else if !isNTT && ct.IsNTT() {
		rlweOp.InvNTTTo(ctOut.Value, ct.Value)
		rlweOp.SubElementTo(ctOut.Value, ctOut.Value, e)
	} else {
		if isNTT {
			pOp.FwdNTTTo(e, e)
		}
		rlweOp.SubElementTo(ctOut.Value, ct.Value, e)
	}

	ctOut.scFac = ct.scFac
}

// SubElement computes ctOut = ct - e.
func (op *Operator) SubElement(ct *Ciphertext, e *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.SubElementTo(ctOut, ct, e, isNTT)
	return ctOut
}

// SubElementTo computes ctOut = ct - e.
func (op *Operator) SubElementTo(ctOut, ct *Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.ModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.rlweOp.PlainOperator()

	ctOut.Resize(ct.ModLen())
	if ct.IsNTT() == e.IsNTT() { // First sub, then perform NTT/iNTT needed.
		op.rlweOp.SubElementTo(ctOut.Value, ct.Value, e)
		if ct.IsNTT() && !isNTT {
			op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		} else if !ct.IsNTT() && isNTT {
			op.rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		}
	} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then add.
		eNTT := op.ePool.Get(e.Type())
		defer op.ePool.Put(eNTT)
		eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

		if e.IsNTT() {
			pOp.InvNTTTo(eNTT, e)
		} else {
			pOp.FwdNTTTo(eNTT, e)
		}
		op.rlweOp.SubElementTo(ctOut.Value, ct.Value, eNTT)
	} else { // First perform NTT/iNTT on ct, then add.
		if ct.IsNTT() {
			op.rlweOp.InvNTTTo(ctOut.Value, ct.Value)
		} else {
			op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
		}
		op.rlweOp.SubElementTo(ctOut.Value, ctOut.Value, e)
	}

	ctOut.scFac = ct.scFac
}

// Mul computes ctOut = ct0 * ct1.
func (op *Operator) Mul(ct0, ct1 *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.MulTo(ctOut, ct0, ct1, rlk, isNTT)
	return ctOut
}

// MulTo computes ctOut = ct0 * ct1.
func (op *Operator) MulTo(ctOut, ct0, ct1 *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) {
	tarLen, auxMod := op.getAuxMod(ct0, ct1)

	// Tensoring the ciphertexts.
	v := op.vPool.Get()
	defer op.vPool.Put(v)
	v = v.WithModLen(tarLen, 0)

	outScFac := op.tensorTo(v, ct0, ct1, tarLen, auxMod, isNTT)

	// Relinearise the result.
	ctOut.Value.Resize(tarLen, 0)
	op.rlweOp.RelinTo(ctOut.Value, v, rlk, isNTT)

	// Rescale if needed.
	ctOut.scFac = outScFac
	op.RescaleTo(ctOut, ctOut, isNTT)
}

// getAuxMod gets the auxiliary modulus for the multiplication.
func (op *Operator) getAuxMod(ct0, ct1 *Ciphertext) (int, *num.Modulus) {
	// Force ct0 to have the smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.scFac < ct1.scFac) {
		ct0, ct1 = ct1, ct0
	}

	curLen := ct0.ModLen()
	tarLen := curLen

	// compute the difference.
	diffFloat := op.fPool.Get()
	tmpFloat := op.fPool.Get()
	diffRound := op.intPool.Get()
	tmpInt := op.intPool.Get()
	defer op.fPool.Put(diffFloat)
	defer op.fPool.Put(tmpFloat)
	defer op.intPool.Put(diffRound)
	defer op.intPool.Put(tmpInt)

	diffFloat.SetFloat64(op.scFac)
	diffFloat.Quo(diffFloat, tmpFloat.SetFloat64(ct0.scFac))
	for tarLen > 0 {
		diffFloat.Mul(diffFloat, tmpFloat.SetUint64(op.params.BaseModulus()[tarLen-1].Value()))

		tmpFloat.Add(diffFloat, tmpFloat.SetFloat64(0.5))
		tmpFloat.Int(diffRound)
		if diffRound.Cmp(tmpInt.SetUint64(1<<num.MaxModulusBits)) <= 0 {
			break
		}

		tarLen--
	}
	if tarLen == 0 {
		panic("ciphertext scaling factor is too large to perform multiplication")
	}

	auxMod := diffRound.Uint64()

	switch {
	case auxMod == op.params.BaseModulus()[tarLen-1].Value():
		return tarLen, nil
	case auxMod == 1:
		return tarLen - 1, nil
	default:
		return tarLen, num.NewModulus(auxMod)
	}
}

// scaleToMulModTo scales the ciphertext to the multiplication modulus.
// Output is in the NTT form.
func (op *Operator) scaleToMulModTo(ctOut *Ciphertext, ctIn *Ciphertext, auxMod *num.Modulus) {
	baseMod := op.params.BaseModulus()

	inLen := ctIn.ModLen()
	outLen := ctOut.ModLen()

	inMod := baseMod[:inLen]
	outMod := make([]*num.Modulus, outLen)
	if auxMod == nil {
		copy(outMod, baseMod[:outLen])
	} else {
		copy(outMod[:outLen-1], baseMod[:outLen-1])
		outMod[outLen-1] = auxMod
	}

	// Find 'diff', which needs to be multiplied to ctIn to make the scaling factors equal.
	diffFloat := op.fPool.Get()
	tmpFloat := op.fPool.Get()
	diffRound := op.intPool.Get()
	tmpInt := op.intPool.Get()
	defer op.fPool.Put(diffFloat)
	defer op.fPool.Put(tmpFloat)
	defer op.intPool.Put(diffRound)
	defer op.intPool.Put(tmpInt)

	diff := op.ePool.Get(crt.TypeScalar)
	defer op.ePool.Put(diff)
	diff = diff.WithModLen(inLen, 0)

	diffFloat.SetFloat64(op.scFac)
	diffFloat.Quo(diffFloat, tmpFloat.SetFloat64(ctIn.scFac))
	for i := 0; i < inLen; i++ {
		diffFloat.Mul(diffFloat, tmpFloat.SetUint64(inMod[i].Value()))
	}
	for i := 0; i < outLen; i++ {
		diffFloat.Quo(diffFloat, tmpFloat.SetUint64(outMod[i].Value()))
	}
	diffFloat.Add(diffFloat, tmpFloat.SetFloat64(0.5))
	diffFloat.Int(diffRound)

	for i := 0; i < inLen; i++ {
		tmpInt.SetUint64(op.params.BaseModulus()[i].Value())
		tmpInt.Rem(diffRound, tmpInt)
		diff.Value.Coeffs[i][0] = tmpInt.Uint64()
	}

	// TODO: optimise later.
	ctBuf := op.ctPool.Get()
	defer op.ctPool.Put(ctBuf)
	ctBuf = ctBuf.WithModLen(inLen)

	op.rlweOp.MulElementTo(ctBuf.Value, ctIn.Value, diff)

	sc := crt.NewScaler(outMod, inMod)
	opOut := crt.NewOperator(op.params.RingParams(), outMod)
	if ctIn.Value.IsNTT() {
		op.InvNTTTo(ctBuf, ctBuf)
	}
	sc.ScaleTo(ctOut.Value.Body.Value, ctBuf.Value.Body.Value)
	sc.ScaleTo(ctOut.Value.Mask.Value, ctBuf.Value.Mask.Value)

	opOut.FwdNTTTo(ctOut.Value.Body.Value, ctOut.Value.Body.Value)
	opOut.FwdNTTTo(ctOut.Value.Mask.Value, ctOut.Value.Mask.Value)

	// Set output scaling factor.
	scale := op.fPool.Get()
	defer op.fPool.Put(scale)

	scale.SetInt(diffRound)
	scale.Mul(scale, tmpFloat.SetFloat64(ctIn.scFac))
	for i := 0; i < inLen; i++ {
		scale.Quo(scale, tmpFloat.SetUint64(inMod[i].Value()))
	}
	for i := 0; i < outLen; i++ {
		scale.Mul(scale, tmpFloat.SetUint64(outMod[i].Value()))
	}

	ctOut.scFac, _ = scale.Float64()
}

// TensorTo tensors two ciphertexts in the multiplication modulus into a vector.
// Input and output are in the NTT form.
func (op *Operator) tensorTo(v *rlwe.Vector, ct0, ct1 *Ciphertext, tarLen int, auxMod *num.Modulus, isNTT bool) float64 {
	var c0, c1 *Ciphertext
	c0 = op.ctPool.Get()
	defer op.ctPool.Put(c0)
	c0 = c0.WithModLen(tarLen)

	// scale to the multiplication modulus.
	op.scaleToMulModTo(c0, ct0, auxMod)
	if ct0 == ct1 {
		c1 = c0
	} else {
		c1 = op.ctPool.Get()
		defer op.ctPool.Put(c1)
		c1 = c1.WithModLen(tarLen)

		op.scaleToMulModTo(c1, ct1, auxMod)
	}

	if v.BaseModLen() != c0.ModLen() || v.BaseModLen() != c1.ModLen() {
		panic("inconsistent input(s)")
	}

	// TODO: optimise later.
	baseMod := op.params.BaseModulus()
	mulMod := make([]*num.Modulus, tarLen)
	if auxMod == nil {
		copy(mulMod, baseMod[:tarLen])
	} else {
		copy(mulMod[:tarLen-1], baseMod[:tarLen-1])
		mulMod[tarLen-1] = auxMod
	}
	opAux := crt.NewOperator(op.params.RingParams(), mulMod)

	opAux.MulTo(v.Value[0].Value, c0.Value.Body.Value, c1.Value.Body.Value)
	opAux.MulTo(v.Value[1].Value, c0.Value.Body.Value, c1.Value.Mask.Value)
	opAux.MulAddTo(v.Value[1].Value, c0.Value.Mask.Value, c1.Value.Body.Value)
	opAux.MulTo(v.Value[2].Value, c0.Value.Mask.Value, c1.Value.Mask.Value)

	// ScaleEmbed to the target modulus.
	op.scaleFromMulModTo(v, v, auxMod, isNTT)

	// compute the scaling factor.
	outScFac := c0.scFac * c1.scFac
	if auxMod != nil {
		outScFac *= float64(baseMod[tarLen-1].Value()) / float64(auxMod.Value())
	}

	return outScFac
}

// scaleFromMulModTo scales and embeds the ciphertext to the target modulus.
// input is in the NTT form.
func (op *Operator) scaleFromMulModTo(vOut *rlwe.Vector, vIn *rlwe.Vector, auxMod *num.Modulus, isNTT bool) {
	if vOut.BaseModLen() != vIn.BaseModLen() {
		panic("inconsistent input(s)")
	}

	// TODO: optimise later.
	if auxMod == nil {
		vOut.Value[0].CopyFrom(vIn.Value[0])
		vOut.Value[1].CopyFrom(vIn.Value[1])
		vOut.Value[2].CopyFrom(vIn.Value[2])
	} else {
		baseMod := op.params.BaseModulus()

		modLen := vIn.BaseModLen()
		inMod := make([]*num.Modulus, modLen)
		copy(inMod[:modLen-1], baseMod[:modLen-1])
		inMod[modLen-1] = auxMod

		outMod := baseMod[:modLen]
		opAux := crt.NewOperator(op.params.RingParams(), inMod)

		sc := crt.NewScaler(outMod, inMod)

		opAux.InvNTTTo(vOut.Value[0].Value, vOut.Value[0].Value)
		opAux.InvNTTTo(vOut.Value[1].Value, vOut.Value[1].Value)
		opAux.InvNTTTo(vOut.Value[2].Value, vOut.Value[2].Value)

		sc.ScaleTo(vOut.Value[0].Value, vOut.Value[0].Value)
		sc.ScaleTo(vOut.Value[1].Value, vOut.Value[1].Value)
		sc.ScaleTo(vOut.Value[2].Value, vOut.Value[2].Value)
	}

	pOp := op.rlweOp.PlainOperator()
	for i := 0; i < 3; i++ {
		if isNTT && !vOut.Value[i].IsNTT() {
			pOp.FwdNTTTo(vOut.Value[i], vOut.Value[i])
		} else if !isNTT && vOut.Value[i].IsNTT() {
			pOp.InvNTTTo(vOut.Value[i], vOut.Value[i])
		}
	}
}

// MulPlain computes ctOut = ct * pt.
func (op *Operator) MulPlain(ct *Ciphertext, pt *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.MulPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// MulPlainTo computes ctOut = ct * pt.
func (op *Operator) MulPlainTo(ctOut, ct *Ciphertext, pt *Plaintext, isNTT bool) {
	var e *rlwe.Element
	switch pt.Rank() {
	case 1:
		e = op.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		e = op.ePool.Get(crt.TypePoly)
	}
	defer op.ePool.Put(e)
	e = e.WithModLen(ct.ModLen(), 0)

	switch pt.valueType {
	case TypeReal:
		op.ecd.EncodeTo(e, pt, op.scFac, true)
	case TypeInt:
		op.ecd.EncodeTo(e, pt, 1, true)
	}

	ctOut.Resize(ct.ModLen())
	if ct.IsNTT() {
		op.rlweOp.MulElementTo(ctOut.Value, ct.Value, e)
	} else {
		op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
		op.rlweOp.MulElementTo(ctOut.Value, ctOut.Value, e)
	}

	if !isNTT {
		op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
	}

	switch pt.valueType {
	case TypeReal:
		ctOut.scFac = ct.scFac * op.scFac
	case TypeInt:
		ctOut.scFac = ct.scFac
	}

	op.RescaleTo(ctOut, ctOut, isNTT)
}

// MulElement computes ctOut = ct * e.
func (op *Operator) MulElement(ct *Ciphertext, e *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(op.params, true)
	op.MulElementTo(ctOut, ct, e, isNTT)
	return ctOut
}

// MulElementTo computes ctOut = ct * e.
func (op *Operator) MulElementTo(ctOut, ct *Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.ModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.rlweOp.PlainOperator()
	ctOut.Resize(ct.ModLen())
	if e.Type() == crt.TypeScalar || (ct.IsNTT() && e.IsNTT()) {
		op.rlweOp.MulElementTo(ctOut.Value, ct.Value, e)
	} else {
		if ct.IsNTT() {
			eNTT := op.ePool.Get(e.Type())
			defer op.ePool.Put(eNTT)
			eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

			pOp.FwdNTTTo(eNTT, e)
			op.rlweOp.MulElementTo(ctOut.Value, ct.Value, eNTT)
		} else {
			op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)

			if e.IsNTT() {
				op.rlweOp.MulElementTo(ctOut.Value, ctOut.Value, e)
			} else {
				eNTT := op.ePool.Get(e.Type())
				defer op.ePool.Put(eNTT)
				eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

				pOp.FwdNTTTo(eNTT, e)
				op.rlweOp.MulElementTo(ctOut.Value, ctOut.Value, eNTT)
			}
		}
	}

	if ctOut.IsNTT() && !isNTT {
		op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
	} else if !ctOut.IsNTT() && isNTT {
		op.rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
	}

	ctOut.scFac = ct.scFac * op.scFac
	op.RescaleTo(ctOut, ctOut, isNTT)
}
