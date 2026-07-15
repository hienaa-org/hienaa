package bfv

import (
	"math"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates operations over [*Ciphertext] and [Plaintext].
type Operator struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	intOp  *heint.Operator
	rlweOp *rlwe.Operator
	ambOp  *rlwe.Operator

	noise *NoiseEstimator

	ePool   *rlwe.ElementPool
	ctPool  *pool.Pool[*Ciphertext]
	vPool   *pool.Pool[*rlwe.Vector]
	embPool *pool.Pool[*[]uint64]
}

// NewOperator creates a new [Operator].
func NewOperator(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *Operator {
	ne := NewNoiseEstimator(params, msgMod, estimType)

	ambParams := rlwe.ParametersLiteral{
		RingParams:  params.RingParams(),
		BaseModulus: ne.ambMod,

		SecretKeyParams: params.SecretKeyParams(),
		NoiseParams:     params.NoiseParams(),
	}.Compile()
	ambOp := rlwe.NewOperator(ambParams)

	return &Operator{
		params: params,
		msgMod: msgMod,

		intOp:  heint.NewOperator(params, msgMod),
		rlweOp: rlwe.NewOperator(params),
		ambOp:  ambOp,

		noise: ne,

		ePool: rlwe.NewElementPool(ambParams, false, true),
		ctPool: pool.NewPool(func() *Ciphertext {
			return NewCiphertext(ambParams, false)
		}),
		vPool: pool.NewPool(func() *rlwe.Vector {
			return rlwe.NewVector(ambParams, 3, false, true)
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
		scale = math.Ceil(math.Sqrt(ct.noise / op.noise.noise.RoundNoise()))
	case heint.WorstCaseType:
		scale = math.Ceil(ct.noise / op.noise.noise.RoundNoise())
	}

	tarLen := ct.ModLen()
	if float64(op.params.BaseModulus()[tarLen-1].Value()) < scale {
		tarLen--
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
	ctOut.Value.Resize(l, 0)
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
	// Force ct0 to have the smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	// Parameters.
	tarLen := ct0.ModLen()
	ambLen := op.noise.getAmbLen(ct0, ct1)

	// Lift and tensor the ciphertexts.
	vAmb := op.vPool.Get()
	defer op.vPool.Put(vAmb)
	vAmb = vAmb.WithModLen(tarLen+ambLen, 0)
	op.liftAndTensorTo(vAmb, ct0, ct1, ambLen, true)

	// Switch the modulus from ambient modulus to the target modulus.
	vBase := vAmb.WithModLen(tarLen, 0)
	op.divRoundTo(vBase, vAmb, true)

	// Relinearise the result.
	ctOut.Value.Resize(vBase.BaseModLen(), 0)
	op.rlweOp.RelinTo(ctOut.Value, vBase, rlk, isNTT)

	// Estimate the noise.
	op.noise.MulTo(ctOut, ct0, ct1)
}

// liftAndTensorTo lifts the ciphertexts to the ambient modulus and tensors them.
func (op *Operator) liftAndTensorTo(vAmb *rlwe.Vector, ct0, ct1 *Ciphertext, ambLen int, isNTT bool) {
	// Force ct0 to have the smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	// Parameters.
	tarLen := ct0.ModLen()

	// Allocate temporary ciphertexts.
	cAmb0 := op.ctPool.Get()
	defer op.ctPool.Put(cAmb0)
	cAmb0 = cAmb0.WithModLen(tarLen + ambLen)

	cAmb1 := op.ctPool.Get()
	defer op.ctPool.Put(cAmb1)
	cAmb1 = cAmb1.WithModLen(tarLen + ambLen)
	cAux1 := cAmb1.WithModLen(ambLen)

	// Modulus switch ct1.
	ctOp := op.ambOp.Params.Operator().Slice(0, ct1.ModLen())
	auxOp := op.ambOp.Params.Operator().Slice(tarLen, tarLen+ambLen)

	sc := crt.NewScaler(auxOp, ctOp).WithPool(op.embPool)
	sc.ScaleTo(cAux1.Value.Body.Value, ct1.Value.Body.Value, true)
	sc.ScaleTo(cAux1.Value.Mask.Value, ct1.Value.Mask.Value, true)

	// Lift to the ambient modulus, in the NTT form.
	op.ambOp.ModRaiseTo(cAmb0.Value, ct0.Value, true)

	ambOp := op.ambOp.Params.Operator().Slice(0, tarLen+ambLen)
	emb := crt.NewEmbedder(ambOp, auxOp).WithPool(op.embPool)
	emb.EmbedTo(cAmb1.Value.Body.Value, cAux1.Value.Body.Value, true)
	emb.EmbedTo(cAmb1.Value.Mask.Value, cAux1.Value.Mask.Value, true)

	// Tensoring the ciphertexts.
	op.ambOp.TensorTo(vAmb, cAmb0.Value, cAmb1.Value)

	// Multiply by the message modulus.
	msgMod := op.ePool.Get(crt.TypeScalar)
	defer op.ePool.Put(msgMod)
	msgMod = msgMod.WithModLen(tarLen+ambLen, 0)
	for i := range msgMod.Value.Coeffs {
		msgMod.Value.Coeffs[i][0] = num.Reduce(op.msgMod.Value(), op.ambOp.Params.FullModulus()[i])
	}

	pAmbOp := op.ambOp.PlainOperator()
	for i := range vAmb.Value {
		pAmbOp.MulTo(vAmb.Value[i], vAmb.Value[i], msgMod)
		if !isNTT {
			pAmbOp.InvNTTTo(vAmb.Value[i], vAmb.Value[i])
		}
	}
}

// divRoundTo divides the vector by the auxMod/msgMod and rounds the result.
// Assumes that vAmb is in the NTT form.
func (op *Operator) divRoundTo(vBase, vAmb *rlwe.Vector, isNTT bool) {
	for i := range vBase.Value {
		op.ambOp.PlainOperator().ScaleTo(vBase.Value[i], vAmb.Value[i], vBase.BaseModLen(), isNTT)
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
}
