package bgv

import (
	"math"
	"sync"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates operations over [*Ciphertext] and [*Plaintext].
type Operator struct {
	Params rlwe.Parameters
	msgMod *num.Modulus

	intOp  *heint.Operator
	rlweOp *rlwe.Operator

	noise *NoiseEstimator

	ePool  *rlwe.ElementPool
	ctPool *sync.Pool
	vPool  *sync.Pool
}

// NewOperator creates a new [Operator].
func NewOperator(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *Operator {
	return &Operator{
		Params: params,
		msgMod: msgMod,

		intOp:  heint.NewOperator(params, msgMod),
		rlweOp: rlwe.NewOperator(params),

		noise: NewNoiseEstimator(params, msgMod, estimType),

		ePool: rlwe.NewElementPool(params, false, true),
		ctPool: &sync.Pool{
			New: func() any {
				return NewCiphertext(params, true)
			},
		},
		vPool: &sync.Pool{
			New: func() any {
				return rlwe.NewVector(params, 3, false, true)
			},
		},
	}
}

// Parameters returns the parameters.
func (op *Operator) Parameters() rlwe.Parameters {
	return op.Params
}

// ModSwitch switches the modulus of the ciphertext to the given length.
func (op *Operator) ModSwitch(ct *Ciphertext, l int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), l, true)
	op.ModSwitchTo(ctOut, ct, l, isNTT)
	return ctOut
}

// ModSwitchTo switches the modulus of the ciphertext to the given length.
func (op *Operator) ModSwitchTo(ctOut, ct *Ciphertext, l int, isNTT bool) {
	op.rlweOp.ScaleTo(ctOut.Value, ct.Value, l, isNTT)
	op.noise.ModSwitchTo(ctOut, ct, l)
}

// FwdNTT returns FwdNTT(ct).
func (op *Operator) FwdNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.FwdNTTTo(ctOut, ct)
	return ctOut
}

// FwdNTTTo computes ctOut = FwdNTT(ct).
func (op *Operator) FwdNTTTo(ctOut, ct *Ciphertext) {
	op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
	op.noise.FwdNTTTo(ctOut, ct)
}

// InvNTT returns InvNTT(ct).
func (op *Operator) InvNTT(ct *Ciphertext) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.InvNTTTo(ctOut, ct)
	return ctOut
}

// InvNTTTo computes ctOut = InvNTT(ct).
func (op *Operator) InvNTTTo(ctOut, ct *Ciphertext) {
	op.rlweOp.InvNTTTo(ctOut.Value, ct.Value)
	op.noise.InvNTTTo(ctOut, ct)
}

// Neg computes ctOut = -ct.
func (op *Operator) Neg(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.NegTo(ctOut, ct, isNTT)
	return ctOut
}

// NegTo computes ctOut = -ct.
func (op *Operator) NegTo(ctOut, ct *Ciphertext, isNTT bool) {
	op.intOp.NegTo(ctOut.Value, ct.Value, isNTT)
	op.noise.NegTo(ctOut, ct)
}

// Add computes ctOut = ct0 + ct1.
func (op *Operator) Add(ct0, ct1 *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), true)
	op.AddTo(ctOut, ct0, ct1, isNTT)
	return ctOut
}

// AddTo computes ctOut = ct0 + ct1.
func (op *Operator) AddTo(ctOut, ct0, ct1 *Ciphertext, isNTT bool) {
	op.intOp.AddTo(ctOut.Value, ct0.Value, ct1.Value, isNTT)
	op.noise.AddTo(ctOut, ct0, ct1)
}

// AddPlain computes ctOut = ct + pt.
func (op *Operator) AddPlain(ct *Ciphertext, pt Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.AddPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// AddPlainTo computes ctOut = ct + pt.
func (op *Operator) AddPlainTo(ctOut, ct *Ciphertext, pt Plaintext, isNTT bool) {
	op.intOp.AddPlainTo(ctOut.Value, ct.Value, pt, isNTT)
	op.noise.AddPlainTo(ctOut, ct, pt)
}

// AddElement computes ctOut = ct + e.
func (op *Operator) AddElement(ct *Ciphertext, e *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.AddElementTo(ctOut, ct, e, isNTT)
	return ctOut
}

// AddElementTo computes ctOut = ct + e.
func (op *Operator) AddElementTo(ctOut, ct *Ciphertext, e *rlwe.Element, isNTT bool) {
	op.intOp.AddElementTo(ctOut.Value, ct.Value, e, isNTT)
	op.noise.AddElementTo(ctOut, ct, e)
}

// Sub computes ctOut = ct0 - ct1.
func (op *Operator) Sub(ct0, ct1 *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), true)
	op.SubTo(ctOut, ct0, ct1, isNTT)
	return ctOut
}

// SubTo computes ctOut = ct0 - ct1.
func (op *Operator) SubTo(ctOut, ct0, ct1 *Ciphertext, isNTT bool) {
	op.intOp.SubTo(ctOut.Value, ct0.Value, ct1.Value, isNTT)
	op.noise.SubTo(ctOut, ct0, ct1)
}

// SubPlain computes ctOut = ct - pt.
func (op *Operator) SubPlain(ct *Ciphertext, pt Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.SubPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// SubPlainTo computes ctOut = ct - pt.
func (op *Operator) SubPlainTo(ctOut, ct *Ciphertext, pt Plaintext, isNTT bool) {
	op.intOp.SubPlainTo(ctOut.Value, ct.Value, pt, isNTT)
	op.noise.SubPlainTo(ctOut, ct, pt)
}

// SubElement computes ctOut = ct - e.
func (op *Operator) SubElement(ct *Ciphertext, e *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.SubElementTo(ctOut, ct, e, isNTT)
	return ctOut
}

// SubElementTo computes ctOut = ct - e.
func (op *Operator) SubElementTo(ctOut, ct *Ciphertext, e *rlwe.Element, isNTT bool) {
	op.intOp.SubElementTo(ctOut.Value, ct.Value, e, isNTT)
	op.noise.SubElementTo(ctOut, ct, e)
}

// Mul computes ctOut = ct0 * ct1.
func (op *Operator) Mul(ct0, ct1 *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), true)
	op.MulTo(ctOut, ct0, ct1, rlk, isNTT)
	return ctOut
}

// MulTo computes ctOut = ct0 * ct1.
func (op *Operator) MulTo(ctOut, ct0, ct1 *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) {
	tarLen, auxIdx, auxMod := op.getAuxMod(ct0, ct1)

	c0 := op.ctPool.Get().(*Ciphertext)
	c1 := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)
	c0 = c0.WithModLen(tarLen)
	c1 = c1.WithModLen(tarLen)

	// scale to the multiplication modulus.
	op.scaleToMulModTo(c0, ct0, auxIdx, auxMod)
	op.scaleToMulModTo(c1, ct1, auxIdx, auxMod)

	// Tensoring the ciphertexts.
	v := op.vPool.Get().(*rlwe.Vector)
	defer op.vPool.Put(v)
	v = v.WithModLen(tarLen, 0)

	op.tensorTo(v, c0, c1, auxIdx, auxMod)

	// ScaleEmbed to the target modulus.
	op.scaleFromMulModTo(v, v, auxIdx, auxMod, isNTT)

	// Relinearise the result.
	op.rlweOp.RelinTo(ctOut.Value, v, rlk, isNTT)

	// Estimate the noise.
	op.noise.MulTo(ctOut, ct0, ct1)
}

// setMulMod sets the modulus for the multiplication.
func (op *Operator) getAuxMod(ct0, ct1 *Ciphertext) (int, int, *num.Modulus) {
	// Force ct0 to have the smaller noise.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise > ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	tarLen := ct0.ModLen()
	auxIdx := 0
	for auxIdx < tarLen-1 {
		if num.GCD(op.Params.BaseModulus()[auxIdx].Value(), op.msgMod.Value()) != 1 {
			break
		}
		auxIdx++
	}

	remInv := uint64(1)
	for i := 0; i < tarLen; i++ {
		if i != auxIdx {
			remInv = num.Mul(remInv, op.Params.BaseModulus()[i].Value(), op.msgMod)
		}
	}
	rem := num.Inv(remInv, op.msgMod)

	// TODO: auxMod should be a coprime with the moduli.
	var scale uint64
	switch op.noise.estimType {
	case heint.VarianceType:
		scale = uint64(math.Ceil(math.Sqrt(ct0.noise / op.noise.noise.RoundNoise())))
	case heint.WorstCaseType:
		scale = uint64(math.Ceil(ct0.noise / op.noise.noise.RoundNoise()))
	}
	auxMod := uint64(math.Floor(float64(op.Params.BaseModulus()[auxIdx].Value())/float64(scale*op.msgMod.Value())))*scale*op.msgMod.Value() + rem

	return tarLen, auxIdx, num.NewModulus(auxMod)
}

// scaleToMulModTo scales the ciphertext to the multiplication modulus.
// Output is in the NTT form.
func (op *Operator) scaleToMulModTo(ctOut *Ciphertext, ctIn *Ciphertext, auxIdx int, auxMod *num.Modulus) {
	inLen := ctIn.ModLen()
	outLen := ctOut.ModLen()

	baseMod := op.Params.BaseModulus()

	inMod := baseMod[:inLen]
	outMod := make([]*num.Modulus, outLen)
	for i := 0; i < outLen; i++ {
		if i == auxIdx {
			outMod[i] = auxMod
		} else {
			outMod[i] = baseMod[i]
		}
	}

	// TODO: optimise later.
	sc := crt.NewScaler(outMod, inMod)
	opOut := crt.NewOperator(op.Params.RingParams(), outMod)
	if ctIn.Value.IsNTT() {
		ctBuf := op.ctPool.Get().(*Ciphertext)
		defer op.ctPool.Put(ctBuf)
		ctBuf = ctBuf.WithModLen(inLen)

		op.InvNTTTo(ctBuf, ctIn)
		sc.ScaleTo(ctOut.Value.Body.Value, ctBuf.Value.Body.Value)
		sc.ScaleTo(ctOut.Value.Mask.Value, ctBuf.Value.Mask.Value)
	} else {
		sc.ScaleTo(ctOut.Value.Body.Value, ctIn.Value.Body.Value)
		sc.ScaleTo(ctOut.Value.Mask.Value, ctIn.Value.Mask.Value)
	}

	opOut.FwdNTTTo(ctOut.Value.Body.Value, ctOut.Value.Body.Value)
	opOut.FwdNTTTo(ctOut.Value.Mask.Value, ctOut.Value.Mask.Value)
}

// TensorTo tensors two ciphertexts in the multiplication modulus into a vector.
// Input and output are in the NTT form.
func (op *Operator) tensorTo(v *rlwe.Vector, c0, c1 *Ciphertext, auxIdx int, auxMod *num.Modulus) {
	if v.BaseModLen() != c0.ModLen() || v.BaseModLen() != c1.ModLen() {
		panic("inconsistent input(s)")
	}

	// TODO: optimise later.
	mulLen := c0.ModLen()
	baseMod := op.Params.BaseModulus()
	mulMod := make([]*num.Modulus, mulLen)
	for i := 0; i < mulLen; i++ {
		if i == auxIdx {
			mulMod[i] = auxMod
		} else {
			mulMod[i] = baseMod[i]
		}
	}
	opAux := crt.NewOperator(op.Params.RingParams(), mulMod)

	msgMod := op.ePool.Get(crt.TypeScalar)
	defer op.ePool.Put(msgMod)
	msgMod = msgMod.WithModLen(mulLen, 0)
	for i := 0; i < mulLen; i++ {
		msgMod.Value.Coeffs[i][0] = num.Reduce(mulMod[i].Value()-op.msgMod.Value(), mulMod[i])
	}

	opAux.MulTo(v.Value[0].Value, c0.Value.Body.Value, c1.Value.Body.Value)
	opAux.MulTo(v.Value[1].Value, c0.Value.Body.Value, c1.Value.Mask.Value)
	opAux.MulAddTo(v.Value[1].Value, c0.Value.Mask.Value, c1.Value.Body.Value)
	opAux.MulTo(v.Value[2].Value, c0.Value.Mask.Value, c1.Value.Mask.Value)

	opAux.MulTo(v.Value[0].Value, v.Value[0].Value, msgMod.Value)
	opAux.MulTo(v.Value[1].Value, v.Value[1].Value, msgMod.Value)
	opAux.MulTo(v.Value[2].Value, v.Value[2].Value, msgMod.Value)
}

// scaleFromMulModTo scales and embeds the ciphertext to the target modulus.
// input is in the NTT form.
func (op *Operator) scaleFromMulModTo(vOut *rlwe.Vector, vIn *rlwe.Vector, auxIdx int, auxMod *num.Modulus, isNTT bool) {
	if vOut.BaseModLen() != vIn.BaseModLen() {
		panic("inconsistent input(s)")
	}

	baseMod := op.Params.BaseModulus()

	// TODO: optimise later.
	modLen := vIn.BaseModLen()
	inMod := make([]*num.Modulus, modLen)
	for i := 0; i < modLen; i++ {
		if i == auxIdx {
			inMod[i] = auxMod
		} else {
			inMod[i] = baseMod[i]
		}
	}
	outMod := baseMod[:modLen]
	opAux := crt.NewOperator(op.Params.RingParams(), inMod)

	sc := crt.NewScaler(outMod, inMod)

	opAux.InvNTTTo(vOut.Value[0].Value, vOut.Value[0].Value)
	opAux.InvNTTTo(vOut.Value[1].Value, vOut.Value[1].Value)
	opAux.InvNTTTo(vOut.Value[2].Value, vOut.Value[2].Value)

	sc.ScaleTo(vOut.Value[0].Value, vOut.Value[0].Value)
	sc.ScaleTo(vOut.Value[1].Value, vOut.Value[1].Value)
	sc.ScaleTo(vOut.Value[2].Value, vOut.Value[2].Value)

	if isNTT {
		pOp := op.rlweOp.PlainOperator()
		pOp.FwdNTTTo(vOut.Value[0], vOut.Value[0])
		pOp.FwdNTTTo(vOut.Value[1], vOut.Value[1])
		pOp.FwdNTTTo(vOut.Value[2], vOut.Value[2])
	}
}

// MulPlain computes ctOut = ct * pt.
func (op *Operator) MulPlain(ct *Ciphertext, pt Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.MulPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// MulPlainTo computes ctOut = ct * pt.
func (op *Operator) MulPlainTo(ctOut, ct *Ciphertext, pt Plaintext, isNTT bool) {
	op.intOp.MulPlainTo(ctOut.Value, ct.Value, pt, isNTT)
	op.noise.MulPlainTo(ctOut, ct, pt)
}

// MulElement computes ctOut = ct * e.
func (op *Operator) MulElement(ct *Ciphertext, e *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.MulElementTo(ctOut, ct, e, isNTT)
	return ctOut
}

// MulElementTo computes ctOut = ct * e.
func (op *Operator) MulElementTo(ctOut, ct *Ciphertext, e *rlwe.Element, isNTT bool) {
	op.intOp.MulElementTo(ctOut.Value, ct.Value, e, isNTT)
	op.noise.MulElementTo(ctOut, ct, e)
}

// KeySwitch performs a key switch and returns the result.
func (op *Operator) KeySwitch(ct *Ciphertext, ksk *rlwe.KeySwitchKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.KeySwitchTo(ctOut, ct, ksk, isNTT)
	return ctOut
}

// KeySwitchTo performs a key switch and stores the result in ctOut.
func (op *Operator) KeySwitchTo(ctOut *Ciphertext, ct *Ciphertext, ksk *rlwe.KeySwitchKey, isNTT bool) {
	op.rlweOp.KeySwitchTo(ctOut.Value, ct.Value, ksk, isNTT)
	op.noise.KeySwitchTo(ctOut, ct)
}

// HoistedKeySwitch performs a hoisted key switch and returns the result.
func (op *Operator) HoistedKeySwitch(decmp *rlwe.Vector, ct *Ciphertext, ksk *rlwe.KeySwitchKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.HoistedKeySwitchTo(ctOut, decmp, ct, ksk, isNTT)
	return ctOut
}

// HoistedKeySwitchTo performs a hoisted key switch and stores the result in ctOut.
func (op *Operator) HoistedKeySwitchTo(ctOut *Ciphertext, decmp *rlwe.Vector, ct *Ciphertext, ksk *rlwe.KeySwitchKey, isNTT bool) {
	op.rlweOp.HoistedKeySwitchTo(ctOut.Value, decmp, ct.Value, ksk, isNTT)
	op.noise.KeySwitchTo(ctOut, ct)
}

// Rotate performs a rotation and returns the result.
func (op *Operator) Rotate(ct *Ciphertext, atk *rlwe.AutomorphismKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.RotateTo(ctOut, ct, atk, isNTT)
	return ctOut
}

// RotateTo performs a rotation and stores the result in ctOut.
func (op *Operator) RotateTo(ctOut *Ciphertext, ct *Ciphertext, rtk *rlwe.AutomorphismKey, isNTT bool) {
	op.AutTo(ctOut, ct, rtk, isNTT)
	op.noise.AutTo(ctOut, ct)
}

// Aut performs an automorphism and returns the result.
func (op *Operator) Aut(ct *Ciphertext, atk *rlwe.AutomorphismKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.AutTo(ctOut, ct, atk, isNTT)
	return ctOut
}

// AutTo performs an automorphism and stores the result in ctOut.
func (op *Operator) AutTo(ctOut *Ciphertext, ct *Ciphertext, atk *rlwe.AutomorphismKey, isNTT bool) {
	op.rlweOp.AutTo(ctOut.Value, ct.Value, atk, isNTT)
	op.noise.AutTo(ctOut, ct)
}

// HoistedAut performs a hoisted automorphism and returns the result.
func (op *Operator) HoistedAut(decmp *rlwe.Vector, ct *Ciphertext, atk *rlwe.AutomorphismKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.HoistedAutTo(ctOut, decmp, ct, atk, isNTT)
	return ctOut
}

// HoistedAutTo performs a hoisted automorphism and stores the result in ctOut.
func (op *Operator) HoistedAutTo(ctOut *Ciphertext, decmp *rlwe.Vector, ct *Ciphertext, atk *rlwe.AutomorphismKey, isNTT bool) {
	op.rlweOp.HoistedAutTo(ctOut.Value, decmp, ct.Value, atk, isNTT)
	op.noise.AutTo(ctOut, ct)
}

// TODO: estimate the noise.
// ExtProd performs an external product and returns the result.
func (op *Operator) ExtProd(ct *Ciphertext, gsw *rlwe.RGSW, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.ExtProdTo(ctOut, ct, gsw, isNTT)
	return ctOut
}

// ExtProdTo performs an external product and stores the result in ctOut.
func (op *Operator) ExtProdTo(ctOut *Ciphertext, ct *Ciphertext, gsw *rlwe.RGSW, isNTT bool) {
	op.rlweOp.ExtProdTo(ctOut.Value, ct.Value, gsw, isNTT)
}

// HoistedExtProd performs a hoisted external product and returns the result.
func (op *Operator) HoistedExtProd(decmpBody *rlwe.Vector, decmpMask *rlwe.Vector, ct *Ciphertext, gsw *rlwe.RGSW, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.HoistedExtProdTo(ctOut, decmpBody, decmpMask, ct, gsw, isNTT)
	return ctOut
}

// HoistedExtProdTo performs a hoisted external product and stores the result in ctOut.
func (op *Operator) HoistedExtProdTo(ctOut *Ciphertext, decmpBody *rlwe.Vector, decmpMask *rlwe.Vector, ct *Ciphertext, gsw *rlwe.RGSW, isNTT bool) {
	op.rlweOp.HoistedExtProdTo(ctOut.Value, decmpBody, decmpMask, gsw, isNTT)
}
