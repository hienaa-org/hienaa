package bfv

import (
	"sync"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

type Encryptor struct {
	params Parameters
	noise  *NoiseEstimator

	pOp *rlwe.PlainOperator

	rlweEnc *rlwe.Encryptor
	encoder *Encoder

	scFacs []*rlwe.Element
	scaler []*crt.Scaler

	sPool  *sync.Pool
	ptPool *sync.Pool
}

// NewEncryptor creates a new [Encryptor].
func NewEncryptor(params Parameters) *Encryptor {
	rlweParams := params.RLWEParams()
	skSampler := rlweParams.SecretKeyParams().Sampler()
	skValue := skSampler.Sample(rlweParams.Rank(), rlweParams.FullModulus())
	sk := (*rlwe.SecretKey)(rlwe.NewElementFrom(skValue, len(rlweParams.AuxModulus())))
	rlweParams.Operator().FwdNTTTo(skValue, skValue)

	return NewEncryptorWithKey(params, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] with a secret key.
func NewEncryptorWithKey(params Parameters, skNTT *rlwe.SecretKey) *Encryptor {
	rlweParams := params.RLWEParams()

	baseMod := rlweParams.BaseModulus()
	msgMod := []*num.Modulus{params.MessageModulus()}
	scFacs := computeScalingFactor(baseMod, msgMod[0])

	scaler := make([]*crt.Scaler, len(baseMod))
	for i := range scaler {
		scaler[i] = crt.NewScaler(msgMod, baseMod[:i+1])
	}

	return &Encryptor{
		params: params,
		noise:  NewNoiseEstimator(params),

		pOp: rlwe.NewPlainOperator(rlweParams),

		rlweEnc: rlwe.NewEncryptorWithKey(rlweParams, skNTT),
		encoder: NewEncoder(params),

		scFacs: scFacs,
		scaler: scaler,

		sPool: &sync.Pool{
			New: func() any {
				return rlwe.NewScalar(rlweParams, true)
			},
		},
		ptPool: &sync.Pool{
			New: func() any {
				return rlwe.NewPoly(rlweParams, rlweParams.HasAuxModulus(), true)
			},
		},
	}
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

// TODO: Rotation

// Encrypt encrypts the message v.
func (e *Encryptor) Encrypt(v *Plaintext, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(e.params, isNTT)
	e.EncryptTo(ctOut, v, isNTT)
	return ctOut
}

// EncryptCustom encrypts the message v with custom parameters.
func (e *Encryptor) EncryptCustom(v *Plaintext, baseLen int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(e.params.rlweParams.Rank(), baseLen, isNTT)
	e.EncryptTo(ctOut, v, isNTT)
	return ctOut
}

// EncryptTo encrypts the message v to ctOut.
func (e *Encryptor) EncryptTo(ctOut *Ciphertext, v *Plaintext, isNTT bool) {
	baseLen := ctOut.ModLen()

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
	e.rlweEnc.EncryptTo(ctOut.Value, pt, isNTT)

	e.noise.EncryptTo(ctOut)
}

// EncryptElement encrypts the element e.
func (e *Encryptor) EncryptElement(eIn *rlwe.Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(e.params.rlweParams.Rank(), eIn.BaseModLen(), isNTT)
	e.EncryptElementTo(ctOut, eIn, isNTT)
	return ctOut
}

// EncryptElementTo encrypts the element e to ctOut.
func (e *Encryptor) EncryptElementTo(ctOut *Ciphertext, eIn *rlwe.Element, isNTT bool) {
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
	e.rlweEnc.EncryptTo(ctOut.Value, (*rlwe.Element)(pt), isNTT)

	e.noise.EncryptTo(ctOut)
}

// Decrypt decrypts the ciphertext ct.
func (e *Encryptor) Decrypt(ct *Ciphertext) *Plaintext {
	vOut := NewPoly(ct.Rank())
	e.DecryptTo(vOut, ct)
	return vOut
}

// DecryptTo decrypts the ciphertext ct to vOut.
func (e *Encryptor) DecryptTo(vOut *Plaintext, ct *Ciphertext) {
	if vOut.Rank() > 1 && vOut.Rank() != ct.Rank() {
		panic("inconsistent output")
	}

	baseLen := ct.ModLen()

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
	vOut.CopyFrom((*Plaintext)(ptScale.Value))
}

// Phase performs Phase(ct).
func (e *Encryptor) Phase(ct *Ciphertext) *rlwe.Element {
	return e.rlweEnc.Phase(ct.Value)
}

// PhaseTo performs Phase(ct) and stores the result in pt.
func (e *Encryptor) PhaseTo(eOut *rlwe.Element, ct *Ciphertext) {
	e.rlweEnc.PhaseTo((*rlwe.Element)(eOut), ct.Value)
}

// Noise returns the noise of the ciphertext.
func (e *Encryptor) Noise(ct *Ciphertext) *rlwe.Element {
	noise := rlwe.NewElement(ct.Rank(), ct.ModLen(), 0, false)
	e.NoiseTo(noise, ct)
	return noise
}

// NoiseTo stores the noise of the ciphertext in eOut.
func (e *Encryptor) NoiseTo(eOut *rlwe.Element, ct *Ciphertext) {
	if !eOut.IsConsistent(ct.Value.Body) {
		panic("inconsistent output")
	}

	baseLen := ct.ModLen()

	pt := e.ptPool.Get().(*rlwe.Element)
	vEcd := e.ptPool.Get().(*rlwe.Element)
	defer e.ptPool.Put(pt)
	defer e.ptPool.Put(vEcd)
	pt = pt.WithModLen(baseLen, 0)
	vEcd = vEcd.WithModLen(baseLen, 0)
	v := vEcd.Value.WithModIdx(1)

	e.PhaseTo(pt, ct)
	e.scaler[baseLen-1].ScaleTo(v, pt.Value)
	e.encoder.EncodeTo(vEcd, (*Plaintext)(v), false)
	e.pOp.MulTo(vEcd, vEcd, e.scFacs[baseLen-1])
	e.pOp.SubTo(eOut, pt, vEcd)
}
