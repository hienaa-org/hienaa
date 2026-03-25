package bgv

import (
	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

type Encryptor struct {
	intEnc *heint.Encryptor
	noise  *NoiseEstimator
}

// NewEncryptor creates a new [Encryptor].
func NewEncryptor(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *Encryptor {
	return &Encryptor{
		intEnc: heint.NewEncryptor(params, msgMod),
		noise:  NewNoiseEstimator(params, msgMod, estimType),
	}
}

// NewEncryptorWithKey creates a new [Encryptor] with a secret key.
func NewEncryptorWithKey(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType, skNTT *rlwe.SecretKey) *Encryptor {
	return &Encryptor{
		intEnc: heint.NewEncryptorWithKey(params, msgMod, skNTT),
		noise:  NewNoiseEstimator(params, msgMod, estimType),
	}
}

// Parameters returns the parameters.
func (enc *Encryptor) Parameters() rlwe.Parameters {
	return enc.intEnc.Parameters()
}

// SecretKey returns the secret key.
func (enc *Encryptor) SecretKey() *rlwe.SecretKey {
	return enc.intEnc.SecretKey()
}

// NewRelinKey creates a new relinearisation key.
func (enc *Encryptor) NewRelinKey() *rlwe.RelinKey {
	return enc.intEnc.NewRelinKey()
}

// NewKeySwitchKey creates a new key switch key.
func (enc *Encryptor) NewKeySwitchKey(skNew *rlwe.SecretKey) *rlwe.KeySwitchKey {
	return enc.intEnc.NewKeySwitchKey(skNew)
}

// NewAutomorphismKey creates a new automorphism key.
func (enc *Encryptor) NewAutomorphismKey(idx int) *rlwe.AutomorphismKey {
	return enc.intEnc.NewAutomorphismKey(idx)
}

// NewRotationKey creates a new rotation key.
func (enc *Encryptor) NewRotationKey(idx []int) *rlwe.AutomorphismKey {
	return enc.intEnc.NewRotationKey(idx)
}

// Encrypt encrypts the message v.
func (enc *Encryptor) Encrypt(v Plaintext, isNTT bool) *Ciphertext {
	ct := NewCiphertext(enc.Parameters(), isNTT)
	enc.EncryptTo(ct, v, isNTT)
	return ct
}

// EncryptCustom encrypts the message v with custom parameters.
func (enc *Encryptor) EncryptCustom(v Plaintext, baseLen int, isNTT bool) *Ciphertext {
	ct := NewCiphertextCustom(enc.Parameters().Rank(), baseLen, isNTT)
	enc.EncryptTo(ct, v, isNTT)
	return ct
}

// EncryptTo encrypts the message v to ctOut.
func (enc *Encryptor) EncryptTo(ctOut *Ciphertext, v Plaintext, isNTT bool) {
	enc.intEnc.EncryptTo(ctOut.Value, v, isNTT)
	enc.noise.EncryptTo(ctOut)
}

// EncryptElement encrypts the element e.
func (enc *Encryptor) EncryptElement(eIn *rlwe.Element, isNTT bool) *Ciphertext {
	ct := NewCiphertext(enc.Parameters(), isNTT)
	enc.EncryptElementTo(ct, eIn, isNTT)
	return ct
}

// EncryptElementTo encrypts the element e to ctOut.
func (enc *Encryptor) EncryptElementTo(ctOut *Ciphertext, eIn *rlwe.Element, isNTT bool) {
	enc.intEnc.EncryptElementTo(ctOut.Value, eIn, isNTT)
	enc.noise.EncryptTo(ctOut)
}

// Decrypt decrypts the ciphertext ct.
func (enc *Encryptor) Decrypt(ct *Ciphertext) Plaintext {
	return enc.intEnc.Decrypt(ct.Value)
}

// DecryptTo decrypts the ciphertext ct to vOut.
func (enc *Encryptor) DecryptTo(vOut Plaintext, ct *Ciphertext) {
	enc.intEnc.DecryptTo(vOut, ct.Value)
}

// Phase performs Phase(ct).
func (enc *Encryptor) Phase(ct *Ciphertext) *rlwe.Element {
	return enc.intEnc.Phase(ct.Value)
}

// PhaseTo performs Phase(ct) and stores the result in pt.
func (enc *Encryptor) PhaseTo(eOut *rlwe.Element, ct *Ciphertext) {
	enc.intEnc.PhaseTo((*rlwe.Element)(eOut), ct.Value)
}

// Noise returns the noise of the ciphertext.
func (enc *Encryptor) Noise(ct *Ciphertext) *rlwe.Element {
	return enc.intEnc.Noise(ct.Value)
}

// NoiseTo stores the noise of the ciphertext in eOut.
func (enc *Encryptor) NoiseTo(eOut *rlwe.Element, ct *Ciphertext) {
	enc.intEnc.NoiseTo(eOut, ct.Value)
}
