package bfv

import (
	"cmp"
	"math"
	"math/big"
	"slices"
	"sync"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates operations over [*Ciphertext] and [*Plaintext].
type Operator struct {
	params Parameters

	rlweOp *rlwe.Operator
	ambOp  *rlwe.Operator

	packer  *Packer
	encoder *Encoder

	scFacs []*rlwe.Element

	noise *NoiseEstimator

	bigPool *sync.Pool
	sPool   *sync.Pool
	ptPool  *sync.Pool
	ctPool  *sync.Pool
	vPool   *sync.Pool
}

// NewOperator creates a new [Operator].
func NewOperator(params Parameters) *Operator {
	rlweParams := params.RLWEParams()
	rlweOp := rlwe.NewOperator(rlweParams)

	ringParams := rlweParams.RingParams()
	extraBits := float64(0)
	for _, mod := range rlweParams.BaseModulus() {
		extraBits += num.Log2(mod.Value())
	}
	extraBits = extraBits + num.Log2(ringParams.ExpFactor())
	extraLen := int(math.Ceil(extraBits / num.MaxModulusBits))
	baseMod := rlweParams.BaseModulus()

	var extraMod []*num.Modulus
	var bitlen float64
	for {
		extraMod = make([]*num.Modulus, extraLen)

		gap, err := dft.NTTPrimeGap(ringParams)
		if err != nil {
			panic(err)
		}

		bitlen = extraBits
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
			bitlen -= num.Log2(primemod.Value())
			prime = num.MustPrevPrime(prime, gap)
			cnt++
		}

		if bitlen > 0 {
			extraLen++
		} else {
			break
		}
	}

	ambBaseMod := append(baseMod, extraMod...)
	ambOp := rlwe.NewOperator(rlwe.ParametersLiteral{
		RingParams:  ringParams,
		BaseModulus: ambBaseMod,

		SecretKeyParams: rlweParams.SecretKeyParams(),
		NoiseParams:     rlweParams.NoiseParams(),
	}.Compile())

	encoder := NewEncoder(params)

	return &Operator{
		params: params,

		rlweOp: rlweOp,
		ambOp:  ambOp,

		packer:  NewPacker(params),
		encoder: encoder,

		scFacs: computeScalingFactor(rlweParams.BaseModulus(), params.MessageModulus()),

		noise: NewNoiseEstimator(params),

		bigPool: &sync.Pool{
			New: func() any {
				return new(big.Int)
			},
		},
		sPool: &sync.Pool{
			New: func() any {
				return rlwe.NewScalar(rlweParams, true)
			},
		},
		ptPool: &sync.Pool{
			New: func() any {
				return rlwe.NewElement(rlweParams.Rank(), len(ambBaseMod), len(rlweParams.AuxModulus()), true)
			},
		},
		ctPool: &sync.Pool{
			New: func() any {
				return NewCiphertextCustom(rlweParams.Rank(), len(ambBaseMod), true)
			},
		},
		vPool: &sync.Pool{
			New: func() any {
				return rlwe.NewVectorCustom(rlweParams.Rank(), len(ambBaseMod), len(rlweParams.AuxModulus()), 3, true)
			},
		},
	}
}

// Parameters returns the parameters.
func (op *Operator) Parameters() Parameters {
	return op.params
}

// AmbientOperator returns the ambient operator.
func (op *Operator) AmbientOperator() *rlwe.Operator {
	return op.ambOp
}

// Packer returns the packer.
func (op *Operator) Packer() *Packer {
	return op.packer
}

// Encoder returns the encoder.
func (op *Operator) Encoder() *Encoder {
	return op.encoder
}

// NoiseEstimator returns the noise estimator.
func (op *Operator) NoiseEstimator() *NoiseEstimator {
	return op.noise
}

// Pack packs a vector of uint64 into a [*Plaintext].
func (op *Operator) Pack(v []uint64) *Plaintext {
	return op.packer.Pack(v)
}

// PackTo packs a vector of uint64 to eOut.
func (op *Operator) PackTo(eOut *Plaintext, v []uint64) {
	op.packer.PackTo(eOut, v)
}

// UnPack unpacks a [*Plaintext] into a vector of uint64.
func (op *Operator) UnPack(e *Plaintext) []uint64 {
	return op.packer.UnPack(e)
}

// UnPackTo unpacks a [*Plaintext] to vOut.
func (op *Operator) UnPackTo(vOut []uint64, e *Plaintext) {
	op.packer.UnPackTo(vOut, e)
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
}

// Neg computes ctOut = -ct.
func (op *Operator) Neg(ct *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.NegTo(ctOut, ct, isNTT)
	return ctOut
}

// NegTo computes ctOut = -ct.
func (op *Operator) NegTo(ctOut, ct *Ciphertext, isNTT bool) {
	op.rlweOp.NegTo(ctOut.Value, ct.Value)
	if isNTT && !ctOut.IsNTT() {
		op.rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
	} else if !isNTT && ctOut.IsNTT() {
		op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
	}

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
	tarLen := min(ct0.ModLen(), ct1.ModLen())

	c0 := op.ctPool.Get().(*Ciphertext)
	c1 := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)

	c0 = c0.WithModLen(tarLen)
	c1 = c1.WithModLen(tarLen)

	op.ModSwitchTo(c0, ct0, tarLen, isNTT)
	op.ModSwitchTo(c1, ct1, tarLen, isNTT)

	op.rlweOp.AddTo(ctOut.Value, c0.Value, c1.Value)

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
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()

	if pt.Rank() == 1 { // Input is a scalar
		e := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(e)
		e = e.WithModLen(ct.ModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.ModLen()-1], e)
		rlweOp.AddElementTo(ctOut.Value, ct.Value, e)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		}
	} else { // Input is a polynomial
		e := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(e)
		e = e.WithModLen(ct.ModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.ModLen()-1], e)
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
	}

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
	if ct.ModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if e.Rank() == 1 { // Input is a scalar
		eEcd := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.ModLen()-1], e)
		rlweOp.AddElementTo(ctOut.Value, ct.Value, eEcd)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		}
	} else { // Input is a polynomial
		eEcd := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.ModLen()-1], e)

		if ct.IsNTT() == e.IsNTT() { // First add, then perform NTT/iNTT needed.
			rlweOp.AddElementTo(ctOut.Value, ct.Value, eEcd)
			if ct.IsNTT() && !isNTT {
				rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
			} else if !ct.IsNTT() && isNTT {
				rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
			}
		} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then add.
			if e.IsNTT() {
				pOp.InvNTTTo(eEcd, eEcd)
			} else {
				pOp.FwdNTTTo(eEcd, eEcd)
			}
			rlweOp.AddElementTo(ctOut.Value, ct.Value, eEcd)
		} else { // First perform NTT/iNTT on ct, then add.
			if ct.IsNTT() {
				rlweOp.InvNTTTo(ctOut.Value, ct.Value)
			} else {
				rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
			}
			rlweOp.AddElementTo(ctOut.Value, ctOut.Value, eEcd)
		}
	}
}

// Sub computes ctOut = ct0 - ct1.
func (op *Operator) Sub(ct0, ct1 *Ciphertext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct0.Rank(), ct0.ModLen(), true)
	op.SubTo(ctOut, ct0, ct1, isNTT)
	return ctOut
}

// SubTo computes ctOut = ct0 - ct1.
func (op *Operator) SubTo(ctOut, ct0, ct1 *Ciphertext, isNTT bool) {
	tarLen := min(ct0.ModLen(), ct1.ModLen())

	c0 := op.ctPool.Get().(*Ciphertext)
	c1 := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)

	c0 = c0.WithModLen(tarLen)
	c1 = c1.WithModLen(tarLen)

	op.ModSwitchTo(c0, ct0, tarLen, isNTT)
	op.ModSwitchTo(c1, ct1, tarLen, isNTT)

	op.rlweOp.SubTo(ctOut.Value, c0.Value, c1.Value)

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
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if pt.Rank() == 1 { // Input is a scalar
		e := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(e)
		e = e.WithModLen(ct.ModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.ModLen()-1], e)
		rlweOp.SubElementTo(ctOut.Value, ct.Value, e)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		}
	} else { // Input is a polynomial
		e := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(e)
		e = e.WithModLen(ct.ModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.ModLen()-1], e)
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
	}

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
	if ct.ModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if e.Rank() == 1 { // Input is a scalar
		eEcd := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.ModLen()-1], e)
		rlweOp.SubElementTo(ctOut.Value, ct.Value, eEcd)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		}
	} else { // Input is a polynomial
		eEcd := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.ModLen()-1], e)

		if ct.IsNTT() == e.IsNTT() { // First sub, then perform NTT/iNTT needed.
			rlweOp.SubElementTo(ctOut.Value, ct.Value, eEcd)
			if ct.IsNTT() && !isNTT {
				rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
			} else if !ct.IsNTT() && isNTT {
				rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
			}
		} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then sub.
			if e.IsNTT() {
				pOp.InvNTTTo(eEcd, eEcd)
			} else {
				pOp.FwdNTTTo(eEcd, eEcd)
			}
			rlweOp.SubElementTo(ctOut.Value, ct.Value, eEcd)
		} else { // First perform NTT/iNTT on ct, then sub.
			if ct.IsNTT() {
				rlweOp.InvNTTTo(ctOut.Value, ct.Value)
			} else {
				rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
			}
			rlweOp.SubElementTo(ctOut.Value, ctOut.Value, eEcd)
		}
	}
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
	ambLen := op.getMulAmbLen(tarLen)

	cAmb0 := op.ctPool.Get().(*Ciphertext)
	cAmb1 := op.ctPool.Get().(*Ciphertext)
	defer op.ctPool.Put(cAmb0)
	defer op.ctPool.Put(cAmb1)

	cAmb0 = cAmb0.WithModLen(ambLen)
	cAmb1 = cAmb1.WithModLen(ambLen)
	c0 := cAmb0.WithModLen(tarLen)
	c1 := cAmb1.WithModLen(tarLen)

	// Modulus switch to the target length.
	op.ModSwitchTo(c0, ct0, tarLen, c0.IsNTT())
	op.ModSwitchTo(c1, ct1, tarLen, c1.IsNTT())

	// Lift to the ambient modulus.
	op.ambOp.ModRaiseTo(cAmb0.Value, c0.Value, true)
	op.ambOp.ModRaiseTo(cAmb1.Value, c1.Value, true)

	// Tensoring the ciphertexts.
	vAmb := op.vPool.Get().(*rlwe.Vector)
	defer op.vPool.Put(vAmb)
	vAmb = vAmb.WithModLen(ambLen, 0)
	vBase := vAmb.WithModLen(tarLen, 0)

	op.tensorTo(vAmb, cAmb0, cAmb1)

	// Switch the modulus from ambient modulus to the target modulus.
	op.divRoundTo(vBase, vAmb)

	// Relinearise the result.
	op.rlweOp.RelinTo(ctOut.Value, vBase, rlk, isNTT)

	// Estimate the noise.
	op.noise.MulTo(ctOut, ct0, ct1)
}

// getMulAmbLen returns the number of ambient moduli required for tensoring.
func (op *Operator) getMulAmbLen(tarLen int) int {
	tarBits := float64(0)
	for i := 0; i < tarLen; i++ {
		tarBits += num.Log2(op.params.rlweParams.BaseModulus()[i].Value())
	}
	tarBits = 2*tarBits + num.Log2(op.params.rlweParams.RingParams().ExpFactor())

	ambModLen := 1
	for {
		tarBits -= num.Log2(op.ambOp.Params.FullModulus()[ambModLen-1].Value())
		if tarBits <= 0 {
			break
		}
		ambModLen += 1
	}

	return ambModLen
}

// tensorTo tensors two ciphertexts in the ambient modulus into a vector.
//
// Input must be in NTT form, and output is in NTT form.
func (op *Operator) tensorTo(vOut *rlwe.Vector, ct0, ct1 *Ciphertext) {
	if vOut.BaseModLen() != ct0.ModLen() || vOut.BaseModLen() != ct1.ModLen() || vOut.AuxModLen() != 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.ambOp.PlainOperator()
	pOp.MulTo(vOut.Value[0], ct0.Value.Body, ct1.Value.Body)
	pOp.MulTo(vOut.Value[1], ct0.Value.Body, ct1.Value.Mask)
	pOp.MulAddTo(vOut.Value[1], ct0.Value.Mask, ct1.Value.Body)
	pOp.MulTo(vOut.Value[2], ct0.Value.Mask, ct1.Value.Mask)
}

// divRoundTo divides a vector by the base modulus and rounds the result.
//
// Input must be in NTT form, and output is in coefficient form.
func (op *Operator) divRoundTo(vOut *rlwe.Vector, vIn *rlwe.Vector) {
	if vOut.AuxModLen() != 0 || vIn.AuxModLen() != 0 {
		panic("input(s) must not have auxiliary modulus")
	} else if vOut.Len() != 3 || vIn.Len() != 3 {
		panic("input(s) must have 3 elements")
	}

	baseLen := vOut.BaseModLen()
	ambLen := vIn.BaseModLen()

	baseMod := op.params.rlweParams.BaseModulus()[:baseLen]
	ambMod := op.ambOp.Params.FullModulus()[:ambLen]

	scale := big.NewRat(int64(op.params.MessageModulus().Value()), 1)

	modi := op.bigPool.Get().(*big.Int)
	defer op.bigPool.Put(modi)
	for i := 0; i < baseLen; i++ {
		modi.SetUint64(baseMod[i].Value())
		scale.Denom().Mul(scale.Denom(), modi)
	}

	pAmb := op.ptPool.Get().(*rlwe.Element)
	defer op.ptPool.Put(pAmb)
	pAmb = pAmb.WithModLen(ambLen, 0)

	opAmb := op.ambOp.PlainOperator()
	scEmb := crt.NewScaleEmbedder(baseMod, ambMod, scale)
	for i := 0; i < 3; i++ {
		opAmb.InvNTTTo(pAmb, vIn.Value[i])
		scEmb.ScaleEmbedTo(vOut.Value[i].Value, pAmb.Value)
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
	rlweOp := op.rlweOp
	if pt.Rank() == 1 { // Input is a scalar
		e := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(e)
		e = e.WithModLen(ct.ModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		rlweOp.MulElementTo(ctOut.Value, ct.Value, e)
		if ct.IsNTT() && !isNTT {
			rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		} else if !ct.IsNTT() && isNTT {
			rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
		}
	} else { // Input is a polynomial
		e := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(e)
		e = e.WithModLen(ct.ModLen(), 0)

		op.encoder.EncodeTo(e, pt, true)
		if ct.IsNTT() {
			op.rlweOp.MulElementTo(ctOut.Value, ct.Value, e)
		} else {
			op.rlweOp.FwdNTTTo(ctOut.Value, ct.Value)
			op.rlweOp.MulElementTo(ctOut.Value, ctOut.Value, e)
		}

		if !isNTT {
			op.rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
		}
	}

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
	if ct.ModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if e.Rank() == 1 || (ct.IsNTT() && e.IsNTT()) {
		rlweOp.MulElementTo(ctOut.Value, ct.Value, e)
	} else {
		if ct.IsNTT() {
			eNTT := op.ptPool.Get().(*rlwe.Element)
			defer op.ptPool.Put(eNTT)
			eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

			pOp.FwdNTTTo(eNTT, e)
			rlweOp.MulElementTo(ctOut.Value, ct.Value, eNTT)
		} else {
			rlweOp.FwdNTTTo(ctOut.Value, ct.Value)

			if e.IsNTT() {
				rlweOp.MulElementTo(ctOut.Value, ctOut.Value, e)
			} else {
				eNTT := op.ptPool.Get().(*rlwe.Element)
				defer op.ptPool.Put(eNTT)
				eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

				pOp.FwdNTTTo(eNTT, e)
				rlweOp.MulElementTo(ctOut.Value, ctOut.Value, eNTT)
			}
		}
	}

	if ctOut.IsNTT() && !isNTT {
		rlweOp.InvNTTTo(ctOut.Value, ctOut.Value)
	} else if !ctOut.IsNTT() && isNTT {
		rlweOp.FwdNTTTo(ctOut.Value, ctOut.Value)
	}

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

// ExtProd performs an external product and returns the result.
func (op *Operator) ExtProd(ct *Ciphertext, gsw *rlwe.RGSW, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(ct.Rank(), ct.ModLen(), true)
	op.ExtProdTo(ctOut, ct, gsw, isNTT)
	return ctOut
}

// ExtProdTo performs an external product and stores the result in ctOut.
func (op *Operator) ExtProdTo(ctOut *Ciphertext, ct *Ciphertext, gsw *rlwe.RGSW, isNTT bool) {
	op.rlweOp.ExtProdTo(ctOut.Value, ct.Value, gsw, isNTT)
	op.noise.ExtProdTo(ctOut, ct)
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
	op.noise.ExtProdTo(ctOut, ct)
}
