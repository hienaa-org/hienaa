package bgv

import (
	"cmp"
	"math"
	"math/big"
	"slices"
	"sync"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates operations over [*Ciphertext] and [*Plaintext].
type Operator struct {
	Params rlwe.Parameters
	msgMod *num.Modulus

	intOp  *heint.Operator
	rlweOp *rlwe.Operator
	ambOp  *rlwe.Operator

	noise *NoiseEstimator

	bigPool *sync.Pool
	ptPool  *sync.Pool
	ctPool  *sync.Pool
	vPool   *sync.Pool
}

// NewOperator creates a new [Operator].
func NewOperator(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *Operator {
	ringParams := params.RingParams()
	baseMod := params.BaseModulus()

	needModBig := big.NewInt(int64(msgMod.Value()))
	modi := big.NewInt(int64(msgMod.Value()))
	for _, mod := range baseMod {
		modi.SetUint64(mod.Value())
		needModBig.Mul(needModBig, modi)
		needModBig.Mul(needModBig, modi.SetInt64(int64(ringParams.ExpFactor())))
	}
	needModBig.Mul(needModBig, modi.SetInt64(int64(ringParams.ExpFactor())))

	extraBits := float64(0)
	for _, mod := range params.BaseModulus() {
		extraBits += num.Log2(mod.Value())
	}
	extraBits = extraBits + num.Log2(ringParams.ExpFactor())
	extraLen := int(math.Ceil(extraBits / num.MaxModulusBits))

	var extraMod []*num.Modulus
	var extraModBig *big.Int
	for {
		extraMod = make([]*num.Modulus, extraLen)
		extraModBig = big.NewInt(1)

		gap, err := dft.NTTPrimeGap(ringParams)
		if err != nil {
			panic(err)
		}

		start := (uint64(math.Floor(num.MaxModulus/float64(gap))))*gap + 1
		if start > num.MaxModulus {
			for start > num.MaxModulus {
				start -= gap
			}
		}
		prime := num.MustPrevPrime(start, gap)

		cnt := 0
		for cnt < extraLen {
			primemod := num.NewModulus(prime)
			for {
				_, ok := slices.BinarySearchFunc(baseMod, primemod, func(a, b *num.Modulus) int {
					return cmp.Compare(a.Value(), b.Value())
				})

				if !ok {
					break
				} else {
					prime = num.MustPrevPrime(prime, gap)
					primemod = num.NewModulus(prime)
				}
			}
			extraMod[cnt] = primemod
			extraModBig.Mul(extraModBig, modi.SetUint64(primemod.Value()))
			prime = num.MustPrevPrime(prime, gap)
			cnt++
		}

		if extraModBig.Cmp(needModBig) < 0 {
			extraLen++
		} else {
			break
		}
	}

	ambBaseMod := append(baseMod, extraMod...)
	ambOp := rlwe.NewOperator(rlwe.ParametersLiteral{
		RingParams:  ringParams,
		BaseModulus: ambBaseMod,

		SecretKeyParams: params.SecretKeyParams(),
		NoiseParams:     params.NoiseParams(),
	}.Compile())

	return &Operator{
		Params: params,
		msgMod: msgMod,

		intOp:  heint.NewOperator(params, msgMod),
		rlweOp: rlwe.NewOperator(params),
		ambOp:  ambOp,

		noise: NewNoiseEstimator(params, msgMod, estimType),

		bigPool: &sync.Pool{
			New: func() any {
				return new(big.Int)
			},
		},
		ptPool: &sync.Pool{
			New: func() any {
				return rlwe.NewElement(params.Rank(), len(ambBaseMod), len(params.AuxModulus()), true)
			},
		},
		ctPool: &sync.Pool{
			New: func() any {
				return NewCiphertextCustom(params.Rank(), len(ambBaseMod), true)
			},
		},
		vPool: &sync.Pool{
			New: func() any {
				return rlwe.NewVectorCustom(params.Rank(), len(ambBaseMod), len(params.AuxModulus()), 3, true)
			},
		},
	}
}

// Parameters returns the parameters.
func (op *Operator) Parameters() rlwe.Parameters {
	return op.Params
}

// AmbientOperator returns the ambient operator.
func (op *Operator) AmbientOperator() *rlwe.Operator {
	return op.ambOp
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
func (op *Operator) AddPlain(ct *Ciphertext, pt *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.AddPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// AddPlainTo computes ctOut = ct + pt.
func (op *Operator) AddPlainTo(ctOut, ct *Ciphertext, pt *Plaintext, isNTT bool) {
	op.intOp.AddPlainTo(ctOut.Value, ct.Value, (*crt.Element)(pt), isNTT)
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
func (op *Operator) SubPlain(ct *Ciphertext, pt *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.SubPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// SubPlainTo computes ctOut = ct - pt.
func (op *Operator) SubPlainTo(ctOut, ct *Ciphertext, pt *Plaintext, isNTT bool) {
	op.intOp.SubPlainTo(ctOut.Value, ct.Value, (*crt.Element)(pt), isNTT)
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
	tarLen := min(ct0.ModLen(), ct1.ModLen())

	mulMod := op.bigPool.Get().(*big.Int)
	defer op.bigPool.Put(mulMod)

	ambLen := op.getMulModAmbLen(ct0, ct1, mulMod)

	cAmb0 := op.ctPool.Get().(*Ciphertext)
	cAmb1 := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(cAmb0)
	defer op.ctPool.Put(cAmb1)

	cAmb0 = cAmb0.WithModLen(ambLen)
	cAmb1 = cAmb1.WithModLen(ambLen)

	// ScaleEmbed to the ambient modulus.
	op.scaleEmbedToAmbTo(cAmb0, ct0, mulMod)
	op.scaleEmbedToAmbTo(cAmb1, ct1, mulMod)

	// Tensoring the ciphertexts.
	vAmb := op.vPool.Get().(*rlwe.Vector)
	defer op.vPool.Put(vAmb)
	vAmb = vAmb.WithModLen(ambLen, 0)
	vBase := vAmb.WithModLen(tarLen, 0)

	op.ambOp.TensorTo(vAmb, cAmb0.Value, cAmb1.Value)

	// ScaleEmbed to the target modulus.
	op.scaleEmbedToTarTo(vBase, vAmb, mulMod)

	// Relinearise the result.
	op.rlweOp.RelinTo(ctOut.Value, vBase, rlk, isNTT)

	// Estimate the noise.
	op.noise.MulTo(ctOut, ct0, ct1)
}

// setMulMod sets the modulus for the multiplication.
func (op *Operator) getMulModAmbLen(ct0, ct1 *Ciphertext, mulMod *big.Int) int {
	scale0 := op.bigPool.Get().(*big.Int)
	scale1 := op.bigPool.Get().(*big.Int)
	defer op.bigPool.Put(scale0)
	defer op.bigPool.Put(scale1)

	noise := new(big.Float)
	switch op.noise.estimType {
	case heint.VarianceType:
		noise.SetFloat64(math.Ceil(math.Sqrt(ct0.noise / op.noise.noise.RoundNoise())))
		noise.Int(scale0)
		noise.SetFloat64(math.Ceil(math.Sqrt(ct1.noise / op.noise.noise.RoundNoise())))
		noise.Int(scale1)
	case heint.WorstCaseType:
		noise.SetFloat64(math.Ceil(ct0.noise / op.noise.noise.RoundNoise()))
		noise.Int(scale0)
		noise.SetFloat64(math.Ceil(ct1.noise / op.noise.noise.RoundNoise()))
		noise.Int(scale1)
	}

	mod0 := op.bigPool.Get().(*big.Int)
	mod1 := op.bigPool.Get().(*big.Int)
	modi := op.bigPool.Get().(*big.Int)
	defer op.bigPool.Put(mod0)
	defer op.bigPool.Put(mod1)
	defer op.bigPool.Put(modi)

	mod0.SetInt64(1)
	for i := 0; i < ct0.ModLen(); i++ {
		modi.SetUint64(op.Params.BaseModulus()[i].Value())
		mod0.Mul(mod0, modi)
	}
	mod0.Quo(mod0, scale0)

	mod1.SetInt64(1)
	for i := 0; i < ct1.ModLen(); i++ {
		modi.SetUint64(op.Params.BaseModulus()[i].Value())
		mod1.Mul(mod1, modi)
	}
	mod1.Quo(mod1, scale1)

	// reuse bigInt.
	msgMod := scale0
	tarMod := scale1
	msgMod.SetUint64(op.msgMod.Value())

	if mod0.Cmp(mod1) > 0 {
		mulMod.Set(mod0)
	} else {
		mulMod.Set(mod1)
	}
	mulMod.Quo(mulMod, msgMod)
	mulMod.Mul(mulMod, msgMod)
	mulMod.Add(mulMod, modi.SetInt64(1))
	tarMod.Mul(mulMod, mulMod)
	tarMod.Mul(tarMod, modi.SetInt64(int64(op.Params.RingParams().ExpFactor())))

	ambLen := 0
	ambMod := mod0 // reuse bigInt.
	ambMod.SetInt64(1)
	for ambMod.Cmp(tarMod) < 0 {
		ambLen += 1
		modi.SetUint64(op.ambOp.Params.FullModulus()[ambLen-1].Value())
		ambMod.Mul(ambMod, modi)
	}

	return ambLen
}

// scaleEmbedForwardTo scales and embeds the ciphertext to the ambient modulus.
func (op *Operator) scaleEmbedToAmbTo(ctOut *Ciphertext, ctIn *Ciphertext, mulMod *big.Int) {
	inLen := ctIn.ModLen()
	outLen := ctOut.ModLen()

	baseMod := op.bigPool.Get().(*big.Int)
	modi := op.bigPool.Get().(*big.Int)
	defer op.bigPool.Put(baseMod)
	defer op.bigPool.Put(modi)
	baseMod.SetInt64(1)

	for i := 0; i < inLen; i++ {
		modi.SetUint64(op.Params.BaseModulus()[i].Value())
		baseMod.Mul(baseMod, modi)
	}

	inMod := op.Params.BaseModulus()[:inLen]
	outMod := op.ambOp.Params.FullModulus()[:outLen]
	scEmb := crt.NewScaleEmbedder(outMod, inMod, new(big.Rat).SetFrac(mulMod, baseMod))

	if ctIn.IsNTT() {
		ctNTT := op.ctPool.Get().(*Ciphertext)
		defer op.ctPool.Put(ctNTT)
		ctNTT = ctNTT.WithModLen(inLen)

		op.rlweOp.InvNTTTo(ctNTT.Value, ctIn.Value)

		scEmb.ScaleEmbedTo(ctOut.Value.Body.Value, ctNTT.Value.Body.Value)
		scEmb.ScaleEmbedTo(ctOut.Value.Mask.Value, ctNTT.Value.Mask.Value)
	} else {
		scEmb.ScaleEmbedTo(ctOut.Value.Body.Value, ctIn.Value.Body.Value)
		scEmb.ScaleEmbedTo(ctOut.Value.Mask.Value, ctIn.Value.Mask.Value)
	}

	ctOut.Value.Body.Value.IsNTT = false
	ctOut.Value.Mask.Value.IsNTT = false
	op.ambOp.FwdNTTTo(ctOut.Value, ctOut.Value)
}

// scaleEmbedBackwardTo scales and embeds the ciphertext to the target modulus.
func (op *Operator) scaleEmbedToTarTo(vBase *rlwe.Vector, vAmb *rlwe.Vector, mulMod *big.Int) {
	inLen := vAmb.BaseModLen()
	outLen := vBase.BaseModLen()

	baseMod := op.bigPool.Get().(*big.Int)
	modi := op.bigPool.Get().(*big.Int)
	defer op.bigPool.Put(baseMod)
	defer op.bigPool.Put(modi)
	baseMod.SetInt64(1)
	for i := 0; i < outLen; i++ {
		modi.SetUint64(op.Params.BaseModulus()[i].Value())
		baseMod.Mul(baseMod, modi)
	}

	inMod := op.ambOp.Params.BaseModulus()[:inLen]
	outMod := op.Params.BaseModulus()[:outLen]
	scEmb := crt.NewScaleEmbedder(outMod, inMod, new(big.Rat).SetFrac(baseMod, mulMod))

	pNTT := op.ptPool.Get().(*rlwe.Element)
	defer op.ptPool.Put(pNTT)
	pNTT = pNTT.WithModLen(inLen, 0)

	msgMod := &rlwe.Element{Value: crt.NewScalarFrom(op.msgMod.Value(), op.ambOp.Params.FullModulus()[:inLen])}
	for i := 0; i < 3; i++ {
		if vAmb.Value[i].IsNTT() {
			op.ambOp.PlainOperator().InvNTTTo(pNTT, vAmb.Value[i])
			op.ambOp.PlainOperator().MulTo(pNTT, pNTT, msgMod)
			scEmb.ScaleEmbedTo(vBase.Value[i].Value, pNTT.Value)
		} else {
			op.ambOp.PlainOperator().MulTo(pNTT, vAmb.Value[i], msgMod)
		}
		scEmb.ScaleEmbedTo(vBase.Value[i].Value, pNTT.Value)
		vBase.Value[i].Value.IsNTT = false
		op.rlweOp.PlainOperator().NegTo(vBase.Value[i], vBase.Value[i])
	}
}

// MulPlain computes ctOut = ct * pt.
func (op *Operator) MulPlain(ct *Ciphertext, pt *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.MulPlainTo(ctOut, ct, pt, isNTT)
	return ctOut
}

// MulPlainTo computes ctOut = ct * pt.
func (op *Operator) MulPlainTo(ctOut, ct *Ciphertext, pt *Plaintext, isNTT bool) {
	op.intOp.MulPlainTo(ctOut.Value, ct.Value, (*crt.Element)(pt), isNTT)
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
