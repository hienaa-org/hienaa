package bgv

import (
	"math"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Operator evaluates operations over [*Ciphertext] and [*Plaintext].
type Operator struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	intOp  *heint.Operator
	rlweOp *rlwe.Operator

	noise *NoiseEstimator

	ePool   *rlwe.ElementPool
	ctPool  *pool.Pool[*Ciphertext]
	vPool   *pool.Pool[*rlwe.Vector]
	embPool *pool.Pool[*[]uint64]
}

// NewOperator creates a new [Operator].
func NewOperator(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *Operator {
	return &Operator{
		params: params,
		msgMod: msgMod,

		intOp:  heint.NewOperator(params, msgMod),
		rlweOp: rlwe.NewOperator(params),

		noise: NewNoiseEstimator(params, msgMod, estimType),

		ePool: rlwe.NewElementPool(params, false, true),
		ctPool: pool.NewPool(func() *Ciphertext {
			return NewCiphertext(params, true)
		}),
		vPool: pool.NewPool(func() *rlwe.Vector {
			return rlwe.NewVector(params, 3, false, true)
		}),
		embPool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.Rank())
			return &v
		}),
	}
}

// Parameters returns the parameters.
func (op *Operator) Parameters() rlwe.Parameters {
	return op.params
}

// Rescale rescales the ciphertext to the target modulus.
func (op *Operator) Rescale(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.RescaleTo(ctOut, ct, isNTT)
	return ctOut
}

// RescaleTo rescales the ciphertext to the target modulus.
func (op *Operator) RescaleTo(ctOut, ct *Ciphertext, isNTT bool) {
	var scale float64
	switch op.noise.estimType {
	case heint.VarianceType:
		scale = math.Sqrt(ct.noise / op.noise.noise.RoundNoise())
	case heint.WorstCaseType:
		scale = ct.noise / op.noise.noise.RoundNoise()
	}

	tarLen := ct.ModLen()
	for tarLen > 0 {
		divMod := float64(op.params.BaseModulus()[tarLen-1].Value())
		if scale < divMod {
			break
		}

		scale /= divMod
		tarLen--
	}
	if tarLen == 0 {
		panic("ciphertext noise is too large")
	}
	op.noise.ModSwitchTo(ctOut, ct, tarLen)

	buf := op.ctPool.Get()
	defer op.ctPool.Put(buf)
	buf = buf.WithModLen(tarLen)

	op.rlweOp.ScaleTo(buf.Value, ct.Value, tarLen, isNTT)
	ctOut.Value.Resize(tarLen, 0)
	ctOut.Value.CopyFrom(buf.Value)
}

// ModSwitch switches the modulus of the ciphertext to the given length.
func (op *Operator) ModSwitch(ct *Ciphertext, l int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), l, true)
	op.ModSwitchTo(ctOut, ct, l, isNTT)
	return ctOut
}

// ModSwitchTo switches the modulus of the ciphertext to the given length.
func (op *Operator) ModSwitchTo(ctOut, ct *Ciphertext, l int, isNTT bool) {
	ctOutTmp := ctOut.WithModLen(l)
	op.rlweOp.ScaleTo(ctOutTmp.Value, ct.Value, l, isNTT)
	ctOut.Value.Resize(l, 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	// Compute the auxiliary modulus for the multiplication.
	tarLen, auxIdx, auxMod := op.noise.getAuxMod(ct0, ct1)

	// Estimate the noise.
	op.noise.MulTo(ctOut, ct0, ct1)

	// Tensoring the ciphertexts.
	v := op.vPool.Get()
	defer op.vPool.Put(v)
	v = v.WithModLen(tarLen, 0)

	op.tensorTo(v, ct0, ct1, tarLen, auxIdx, auxMod, true)

	// Relinearise the result.
	ctOut.Value.Resize(tarLen, 0)
	op.rlweOp.RelinTo(ctOut.Value, v, rlk, true)

	// Rescale if needed.
	op.RescaleTo(ctOut, ctOut, isNTT)
}

// scaleToMulModTo scales the ciphertext to the multiplication modulus.
// Output is in the NTT form.
func (op *Operator) scaleToMulModTo(ctOut *Ciphertext, ctIn *Ciphertext, auxIdx int, auxMod *num.Modulus) {
	inLen := ctIn.ModLen()
	outLen := ctOut.ModLen()

	auxLen := len(op.rlweOp.PlainOperator().Params.AuxModulus())
	opIn := op.rlweOp.PlainOperator().Params.Operator().WithModIdx(vec.Range(auxLen, auxLen+inLen)...)

	var opOut *crt.Operator
	if auxMod == nil {
		opOut = opIn.WithModIdx(vec.Range(0, outLen)...)
	} else {
		opOut = opIn.WithModIdx(vec.Range(0, auxIdx)...).AppendTmpModulus(auxMod).Append(opIn.WithModIdx(vec.Range(auxIdx+1, outLen)...))
	}

	sc := crt.NewScaler(opOut, opIn).WithPool(op.embPool)
	sc.ScaleTo(ctOut.Value.Body.Value, ctIn.Value.Body.Value, true)
	sc.ScaleTo(ctOut.Value.Mask.Value, ctIn.Value.Mask.Value, true)
}

// TensorTo tensors two ciphertexts in the multiplication modulus into a vector.
//
// Input and output are in the NTT form.
// The algorithm works for any modulus, however, the noise is minimized only when Q mod t = +-1.
func (op *Operator) tensorTo(v *rlwe.Vector, ct0, ct1 *Ciphertext, tarLen int, auxIdx int, auxMod *num.Modulus, isNTT bool) {
	// scale to the multiplication modulus.
	var c0, c1 *Ciphertext
	c0 = op.ctPool.Get()
	defer op.ctPool.Put(c0)
	c0 = c0.WithModLen(tarLen)

	op.scaleToMulModTo(c0, ct0, auxIdx, auxMod)
	if ct0 == ct1 {
		c1 = c0
	} else {
		c1 = op.ctPool.Get()
		defer op.ctPool.Put(c1)
		c1 = c1.WithModLen(tarLen)
		op.scaleToMulModTo(c1, ct1, auxIdx, auxMod)
	}

	if v.BaseModLen() != c0.ModLen() || v.BaseModLen() != c1.ModLen() {
		panic("inconsistent input(s)")
	}

	// Compute required constants.
	auxLen := len(op.rlweOp.PlainOperator().Params.AuxModulus())
	baseOp := op.rlweOp.PlainOperator().Params.Operator().WithModIdx(vec.Range(auxLen, auxLen+tarLen)...)
	var opAux *crt.Operator
	if auxMod == nil {
		opAux = baseOp
	} else {
		opAux = baseOp.WithModIdx(vec.Range(0, auxIdx)...).AppendTmpModulus(auxMod).Append(baseOp.WithModIdx(vec.Range(auxIdx+1, tarLen)...))
	}

	// Tensoring the ciphertexts.
	opAux.MulTo(v.Value[0].Value, c0.Value.Body.Value, c1.Value.Body.Value)
	opAux.MulTo(v.Value[1].Value, c0.Value.Body.Value, c1.Value.Mask.Value)
	opAux.MulAddTo(v.Value[1].Value, c0.Value.Mask.Value, c1.Value.Body.Value)
	opAux.MulTo(v.Value[2].Value, c0.Value.Mask.Value, c1.Value.Mask.Value)

	// Compute mulConst = t * [-Q^{-1}]_t mod Q.
	mulConst := op.ePool.Get(crt.TypeScalar)
	defer op.ePool.Put(mulConst)
	mulConst = mulConst.WithModLen(tarLen, 0)

	mulMod := opAux.Modulus()
	rem := uint64(1)
	for i := 0; i < tarLen; i++ {
		rem = num.Mul(rem, mulMod[i].Value(), op.msgMod)
	}
	remInv := num.Inv(rem, op.msgMod)

	if remInv <= op.msgMod.Value()>>1 {
		for i := 0; i < tarLen; i++ {
			mulConst.Value.Coeffs[i][0] = num.Mul(remInv, op.msgMod.Value(), mulMod[i])
		}
		opAux.NegTo(mulConst.Value, mulConst.Value)
	} else {
		for i := 0; i < tarLen; i++ {
			mulConst.Value.Coeffs[i][0] = num.Mul(op.msgMod.Value()-remInv, op.msgMod.Value(), mulMod[i])
		}
	}

	// Multiply by mulConst.
	opAux.MulTo(v.Value[0].Value, v.Value[0].Value, mulConst.Value)
	opAux.MulTo(v.Value[1].Value, v.Value[1].Value, mulConst.Value)
	opAux.MulTo(v.Value[2].Value, v.Value[2].Value, mulConst.Value)

	// ScaleEmbed to the target modulus.
	op.scaleFromMulModTo(v, v, auxIdx, auxMod, isNTT)
}

// scaleFromMulModTo scales and embeds the ciphertext to the target modulus.
// input is in the NTT form.
func (op *Operator) scaleFromMulModTo(vOut *rlwe.Vector, vIn *rlwe.Vector, auxIdx int, auxMod *num.Modulus, isNTT bool) {
	if vOut.BaseModLen() != vIn.BaseModLen() {
		panic("inconsistent input(s)")
	}

	inLen := vIn.BaseModLen()
	outLen := vOut.BaseModLen()
	auxLen := len(op.rlweOp.PlainOperator().Params.AuxModulus())
	opOut := op.rlweOp.PlainOperator().Params.Operator().WithModIdx(vec.Range(auxLen, auxLen+outLen)...)

	var opIn *crt.Operator
	if auxMod == nil {
		opIn = opOut
	} else {
		opIn = opOut.WithModIdx(vec.Range(0, auxIdx)...).AppendTmpModulus(auxMod).Append(opOut.WithModIdx(vec.Range(auxIdx+1, inLen)...))
	}

	sc := crt.NewScaler(opOut, opIn).WithPool(op.embPool)
	sc.ScaleTo(vOut.Value[0].Value, vIn.Value[0].Value, isNTT)
	sc.ScaleTo(vOut.Value[1].Value, vIn.Value[1].Value, isNTT)
	sc.ScaleTo(vOut.Value[2].Value, vIn.Value[2].Value, isNTT)
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
	op.RescaleTo(ctOut, ctOut, isNTT)
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
	op.RescaleTo(ctOut, ctOut, isNTT)
}

// KeySwitch performs a key switch and returns the result.
func (op *Operator) KeySwitch(ct *Ciphertext, ksk *rlwe.KeySwitchKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.KeySwitchTo(ctOut, ct, ksk, isNTT)
	return ctOut
}

// KeySwitchTo performs a key switch and stores the result in ctOut.
func (op *Operator) KeySwitchTo(ctOut *Ciphertext, ct *Ciphertext, ksk *rlwe.KeySwitchKey, isNTT bool) {
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
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
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
	op.rlweOp.HoistedExtProdTo(ctOut.Value, decmpBody, decmpMask, gsw, isNTT)
}

// MulPlainMatrix computes ctOut = mat * ct.
func (op *Operator) MulPlainMatrix(mat *rlwe.PlainMatrix, ct *Ciphertext, atk map[int]*rlwe.AutomorphismKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.MulPlainMatrixTo(ctOut, ct, mat, atk, isNTT)
	return ctOut
}

// MulPlainMatrixTo computes ctOut = mat * ct.
func (op *Operator) MulPlainMatrixTo(ctOut, ct *Ciphertext, mat *rlwe.PlainMatrix, atk map[int]*rlwe.AutomorphismKey, isNTT bool) {
	ctOut.Value.Resize(ct.Value.BaseModLen(), 0)
	op.noise.MulPlainMatrixTo(ctOut, ct, mat)
	op.rlweOp.MulPlainMatrixTo(ctOut.Value, ct.Value, mat, atk, isNTT)
	op.RescaleTo(ctOut, ctOut, isNTT)
}
