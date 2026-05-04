package ckks

import (
	"github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
)

type Encryptor struct {
	params rlwe.Parameters

	pOp     *rlwe.PlainOperator
	rlweEnc *rlwe.Encryptor

	pack pack.ComplexPacker
	ecd  *Encoder

	ePool *rlwe.ElementPool
}

// NewEncryptor creates a new [Encryptor].
func NewEncryptor(params rlwe.Parameters) *Encryptor {
	skSampler := params.SecretKeyParams().Sampler()
	skValue := skSampler.Sample(params.Rank(), params.FullModulus())
	sk := (*rlwe.SecretKey)(rlwe.NewElementFrom(skValue, len(params.AuxModulus())))
	params.Operator().FwdNTTTo(skValue, skValue)

	return NewEncryptorWithKey(params, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] with a secret key.
func NewEncryptorWithKey(params rlwe.Parameters, skNTT *rlwe.SecretKey) *Encryptor {
	return &Encryptor{
		params: params,

		pOp:     rlwe.NewPlainOperator(params),
		rlweEnc: rlwe.NewEncryptorWithKey(params, skNTT),

		pack: pack.NewComplexPacker(params.RingParams()),
		ecd:  NewEncoder(params),

		ePool: rlwe.NewElementPool(params, params.HasAuxModulus(), true),
	}
}

// Parameters returns the parameters.
func (enc *Encryptor) Parameters() rlwe.Parameters {
	return enc.params
}

// SecretKey returns the secret key.
func (enc *Encryptor) SecretKey() *rlwe.SecretKey {
	return enc.rlweEnc.SecretKey()
}

// NewRelinKey creates a new relinearisation key.
func (enc *Encryptor) NewRelinKey() *rlwe.RelinKey {
	return enc.rlweEnc.NewRelinKey()
}

// NewKeySwitchKey creates a new key switch key.
func (enc *Encryptor) NewKeySwitchKey(skNew *rlwe.SecretKey) *rlwe.KeySwitchKey {
	return enc.rlweEnc.NewKeySwitchKey(skNew)
}

// NewAutomorphismKey creates a new automorphism key.
func (enc *Encryptor) NewAutomorphismKey(idx int) *rlwe.AutomorphismKey {
	return enc.rlweEnc.NewAutomorphismKey(idx)
}

// NewRotationKey creates a new rotation key.
func (enc *Encryptor) NewRotationKey(idx []int) *rlwe.AutomorphismKey {
	auxIdx := enc.pack.RotIdxToAutIdx(idx)
	return enc.rlweEnc.NewAutomorphismKey(auxIdx)
}

// Encrypt encrypts the message v.
func (enc *Encryptor) Encrypt(v *Plaintext, scFac float64, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(enc.params, isNTT)
	enc.EncryptTo(ctOut, v, scFac, isNTT)
	return ctOut
}

// EncryptCustom encrypts the message v with custom parameters.
func (enc *Encryptor) EncryptCustom(v *Plaintext, scFac float64, baseLen int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(enc.params.Rank(), baseLen, isNTT)
	enc.EncryptTo(ctOut, v, scFac, isNTT)
	return ctOut
}

// EncryptTo encrypts the message v to ctOut.
func (enc *Encryptor) EncryptTo(ctOut *Ciphertext, v *Plaintext, scFac float64, isNTT bool) {
	if ctOut.Value.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	}

	baseLen := ctOut.Value.BaseModLen()
	var pt *rlwe.Element
	switch v.Rank() {
	case 1:
		pt = enc.ePool.Get(crt.TypeScalar)
	case ctOut.Rank():
		pt = enc.ePool.Get(crt.TypePoly)
	default:
		panic("invalid input length")
	}
	defer enc.ePool.Put(pt)
	pt = pt.WithModLen(baseLen, 0)

	enc.ecd.EncodeTo(pt, v, scFac, isNTT)
	enc.rlweEnc.EncryptTo(ctOut.Value, pt, isNTT)
	ctOut.scFac = scFac
}

// EncryptElement encrypts the element e.
// TODO: does it make sense to have scFac here?
func (enc *Encryptor) EncryptElement(eIn *rlwe.Element, scFac float64, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(enc.params.Rank(), eIn.BaseModLen(), isNTT)
	enc.EncryptElementTo(ctOut, eIn, scFac, isNTT)
	return ctOut
}

// EncryptElementTo encrypts the element e to ctOut.
func (enc *Encryptor) EncryptElementTo(ctOut *Ciphertext, eIn *rlwe.Element, scFac float64, isNTT bool) {
	if eIn.AuxModLen() > 0 || ctOut.Value.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	} else if eIn.BaseModLen() != ctOut.Value.BaseModLen() {
		panic("base modulus length mismatch")
	}

	enc.rlweEnc.EncryptTo(ctOut.Value, eIn, isNTT)
	ctOut.scFac = scFac
}

// Decrypt decrypts the ciphertext ct.
func (enc *Encryptor) Decrypt(ct *Ciphertext) *Plaintext {
	vOut := NewPoly(ct.Rank(), TypeReal)
	enc.DecryptTo(vOut, ct)
	return vOut
}

// DecryptTo decrypts the ciphertext ct to vOut.
func (enc *Encryptor) DecryptTo(vOut *Plaintext, ct *Ciphertext) {
	if vOut.Rank() != 1 && vOut.Rank() != ct.Rank() {
		panic("inconsistent output")
	}

	baseLen := ct.Value.BaseModLen()
	var pt *rlwe.Element
	switch vOut.Rank() {
	case 1:
		pt = enc.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		pt = enc.ePool.Get(crt.TypePoly)
	default:
		panic("invalid input length")
	}
	defer enc.ePool.Put(pt)
	pt = pt.WithModLen(baseLen, 0)

	// TODO: refrain from heavy bigInt vector allocation.
	enc.PhaseTo(pt, ct)
	vBig := enc.pOp.AsBig(pt)
	for i := range vOut.Value {
		val, _ := vBig[i].Float64()
		vOut.Value[i] = val / ct.scFac
	}
}

// Phase performs Phase(ct).
func (enc *Encryptor) Phase(ct *Ciphertext) *rlwe.Element {
	return enc.rlweEnc.Phase(ct.Value)
}

// PhaseTo performs Phase(ct) and stores the result in eOut.
func (enc *Encryptor) PhaseTo(eOut *rlwe.Element, ct *Ciphertext) {
	enc.rlweEnc.PhaseTo(eOut, ct.Value)
}
