package heint

import (
	"sync"

	"github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encryptor encrypts/decrypts [*Ciphertext] and [*Plaintext].
type Encryptor struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	pOp     *rlwe.PlainOperator
	rlweEnc *rlwe.Encryptor
	packer  pack.IntPacker
	encoder *Encoder

	scFacs []*rlwe.Element
	scaler []*crt.Scaler

	sPool  *sync.Pool
	ptPool *sync.Pool
}

// NewEncryptor creates a new [Encryptor].
func NewEncryptor(params rlwe.Parameters, msgMod *num.Modulus) *Encryptor {
	skSampler := params.SecretKeyParams().Sampler()
	skValue := skSampler.Sample(params.Rank(), params.FullModulus())
	sk := (*rlwe.SecretKey)(rlwe.NewElementFrom(skValue, len(params.AuxModulus())))
	params.Operator().FwdNTTTo(skValue, skValue)

	return NewEncryptorWithKey(params, msgMod, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] with a secret key.
func NewEncryptorWithKey(params rlwe.Parameters, msgMod *num.Modulus, skNTT *rlwe.SecretKey) *Encryptor {
	baseMod := params.BaseModulus()
	scFacs := computeScalingFactor(baseMod, msgMod)

	scaler := make([]*crt.Scaler, len(baseMod))
	for i := range scaler {
		scaler[i] = crt.NewScaler([]*num.Modulus{msgMod}, baseMod[:i+1])
	}

	return &Encryptor{
		params: params,
		msgMod: msgMod,

		pOp:     rlwe.NewPlainOperator(params),
		rlweEnc: rlwe.NewEncryptorWithKey(params, skNTT),
		packer:  pack.NewIntPacker(params.RingParams(), msgMod),
		encoder: NewEncoder(params, msgMod),

		scFacs: scFacs,
		scaler: scaler,

		sPool: &sync.Pool{
			New: func() any {
				return rlwe.NewScalar(params, true)
			},
		},
		ptPool: &sync.Pool{
			New: func() any {
				return rlwe.NewPoly(params, params.HasAuxModulus(), true)
			},
		},
	}
}

// Parameters returns the parameters.
func (e *Encryptor) Parameters() rlwe.Parameters {
	return e.params
}

// SecretKey returns the secret key.
func (e *Encryptor) SecretKey() *rlwe.SecretKey {
	return e.rlweEnc.SecretKey()
}

// NewRelinKey creates a new relinearisation key.
func (e *Encryptor) NewRelinKey() *rlwe.RelinKey {
	return e.rlweEnc.NewRelinKey()
}

// NewKeySwitchKey creates a new key switch key.
func (e *Encryptor) NewKeySwitchKey(skNew *rlwe.SecretKey) *rlwe.KeySwitchKey {
	return e.rlweEnc.NewKeySwitchKey(skNew)
}

// NewAutomorphismKey creates a new automorphism key.
func (e *Encryptor) NewAutomorphismKey(idx int) *rlwe.AutomorphismKey {
	return e.rlweEnc.NewAutomorphismKey(idx)
}

// NewRotationKey creates a new rotation key.
func (e *Encryptor) NewRotationKey(idx []int) *rlwe.AutomorphismKey {
	autIdx := e.packer.RotIdxToAutIdx(idx)
	return e.rlweEnc.NewAutomorphismKey(autIdx)
}

// Encrypt encrypts the message v.
func (e *Encryptor) Encrypt(v *crt.Element, isNTT bool) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertext(e.params, false, isNTT)
	e.EncryptTo(ctOut, v, isNTT)
	return ctOut
}

// EncryptCustom encrypts the message v with custom parameters.
func (e *Encryptor) EncryptCustom(v *crt.Element, baseLen int, isNTT bool) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertextCustom(e.params.Rank(), baseLen, 0, isNTT)
	e.EncryptTo(ctOut, v, isNTT)
	return ctOut
}

// EncryptTo encrypts the message v to ctOut.
func (e *Encryptor) EncryptTo(ctOut *rlwe.Ciphertext, v *crt.Element, isNTT bool) {
	if ctOut.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	} else if v.ModLen() != 1 {
		panic("message modulus length should be one")
	}

	baseLen := ctOut.BaseModLen()

	var pt *rlwe.Element
	if v.Rank() == 1 {
		pt = e.sPool.Get().(*rlwe.Element)
		defer e.sPool.Put(pt)
		pt = pt.WithModLen(baseLen, 0)
	} else {
		pt = e.ptPool.Get().(*rlwe.Element)
		defer e.ptPool.Put(pt)
		pt = pt.WithModLen(baseLen, 0)
	}

	e.encoder.EncodeTo(pt, v, isNTT)
	e.pOp.MulTo(pt, pt, e.scFacs[baseLen-1])
	e.rlweEnc.EncryptTo(ctOut, pt, isNTT)
}

// EncryptElement encrypts the element e.
func (e *Encryptor) EncryptElement(eIn *rlwe.Element, isNTT bool) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertextCustom(e.params.Rank(), eIn.BaseModLen(), 0, isNTT)
	e.EncryptElementTo(ctOut, eIn, isNTT)
	return ctOut
}

// EncryptElementTo encrypts the element e to ctOut.
func (e *Encryptor) EncryptElementTo(ctOut *rlwe.Ciphertext, eIn *rlwe.Element, isNTT bool) {
	if eIn.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	}

	baseLen := eIn.BaseModLen()

	var pt *rlwe.Element
	if eIn.Rank() == 1 {
		pt = e.sPool.Get().(*rlwe.Element)
		defer e.sPool.Put(pt)
		pt = pt.WithModLen(baseLen, 0)
	} else {
		pt = e.ptPool.Get().(*rlwe.Element)
		defer e.ptPool.Put(pt)
		pt = pt.WithModLen(baseLen, 0)
	}

	pt.CopyFrom(eIn)
	e.pOp.MulTo((*rlwe.Element)(pt), (*rlwe.Element)(pt), e.scFacs[baseLen-1])
	e.rlweEnc.EncryptTo(ctOut, (*rlwe.Element)(pt), isNTT)
}

// Decrypt decrypts the ciphertext ct.
func (e *Encryptor) Decrypt(ct *rlwe.Ciphertext) *crt.Element {
	vOut := crt.NewPoly(ct.Rank(), 1)
	e.DecryptTo(vOut, ct)
	return vOut
}

// DecryptTo decrypts the ciphertext ct to vOut.
func (e *Encryptor) DecryptTo(vOut *crt.Element, ct *rlwe.Ciphertext) {
	if vOut.Rank() > 1 && vOut.Rank() != ct.Rank() {
		panic("inconsistent output")
	}

	baseLen := ct.BaseModLen()

	var pt *rlwe.Element
	if vOut.Rank() == 1 {
		pt = e.sPool.Get().(*rlwe.Element)
		defer e.sPool.Put(pt)
		pt = pt.WithModLen(baseLen, 0)
	} else {
		pt = e.ptPool.Get().(*rlwe.Element)
		defer e.ptPool.Put(pt)
		pt = pt.WithModLen(baseLen, 0)
	}
	ptScale := pt.WithModLen(1, 0)

	e.PhaseTo(pt, ct)

	e.scaler[baseLen-1].ScaleTo(ptScale.Value, pt.Value)
	vOut.CopyFrom(ptScale.Value)
}

// Phase performs Phase(ct).
func (e *Encryptor) Phase(ct *rlwe.Ciphertext) *rlwe.Element {
	return e.rlweEnc.Phase(ct)
}

// PhaseTo performs Phase(ct) and stores the result in pt.
func (e *Encryptor) PhaseTo(eOut *rlwe.Element, ct *rlwe.Ciphertext) {
	e.rlweEnc.PhaseTo((*rlwe.Element)(eOut), ct)
}

// Noise returns the noise of the ciphertext.
func (e *Encryptor) Noise(ct *rlwe.Ciphertext) *rlwe.Element {
	noise := rlwe.NewElement(ct.Rank(), ct.BaseModLen(), 0, false)
	e.NoiseTo(noise, ct)
	return noise
}

// NoiseTo stores the noise of the ciphertext in eOut.
func (e *Encryptor) NoiseTo(eOut *rlwe.Element, ct *rlwe.Ciphertext) {
	if !eOut.IsConsistent(ct.Body) {
		panic("inconsistent output")
	}

	baseLen := ct.BaseModLen()

	pt := e.ptPool.Get().(*rlwe.Element)
	vEcd := e.ptPool.Get().(*rlwe.Element)
	defer e.ptPool.Put(pt)
	defer e.ptPool.Put(vEcd)
	pt = pt.WithModLen(baseLen, 0)
	vEcd = vEcd.WithModLen(baseLen, 0)
	v := vEcd.Value.WithModIdx(1)

	e.PhaseTo(pt, ct)
	e.scaler[baseLen-1].ScaleTo(v, pt.Value)
	e.encoder.EncodeTo(vEcd, v, false)
	e.pOp.MulTo(vEcd, vEcd, e.scFacs[baseLen-1])
	e.pOp.SubTo(eOut, pt, vEcd)
}
