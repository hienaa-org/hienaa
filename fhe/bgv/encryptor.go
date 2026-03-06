package bgv

import (
	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

type Encryptor struct {
	heintEnc *heint.Encryptor
	noise    *NoiseEstimator
}

// NewEncryptor creates a new [Encryptor].
func NewEncryptor(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *Encryptor {
	return &Encryptor{
		heintEnc: heint.NewEncryptor(params, msgMod),
		noise:    NewNoiseEstimator(params, msgMod, estimType),
	}
}

// NewEncryptorWithKey creates a new [Encryptor] with a secret key.
func NewEncryptorWithKey(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType, skNTT *rlwe.SecretKey) *Encryptor {
	return &Encryptor{
		heintEnc: heint.NewEncryptorWithKey(params, msgMod, skNTT),
		noise:    NewNoiseEstimator(params, msgMod, estimType),
	}
}

// Parameters returns the parameters.
func (e *Encryptor) Parameters() rlwe.Parameters {
	return e.heintEnc.Parameters()
}

// SecretKey returns the secret key.
func (e *Encryptor) SecretKey() *rlwe.SecretKey {
	return e.heintEnc.SecretKey()
}

// NewRelinKey creates a new relinearisation key.
func (e *Encryptor) NewRelinKey() *rlwe.RelinKey {
	return e.heintEnc.NewRelinKey()
}

// NewKeySwitchKey creates a new key switch key.
func (e *Encryptor) NewKeySwitchKey(skNew *rlwe.SecretKey) *rlwe.KeySwitchKey {
	return e.heintEnc.NewKeySwitchKey(skNew)
}

// NewAutomorphismKey creates a new automorphism key.
func (e *Encryptor) NewAutomorphismKey(idx int) *rlwe.AutomorphismKey {
	return e.heintEnc.NewAutomorphismKey(idx)
}

// NewRotationKey creates a new rotation key.
func (e *Encryptor) NewRotationKey(idx []int) *rlwe.AutomorphismKey {
	return e.heintEnc.NewRotationKey(idx)
}

// Encrypt encrypts the message v.
func (e *Encryptor) Encrypt(v *Plaintext, isNTT bool) *Ciphertext {
	ct := NewCiphertext(e.Parameters(), isNTT)
	e.EncryptTo(ct, v, isNTT)
	return ct
}

// EncryptCustom encrypts the message v with custom parameters.
func (e *Encryptor) EncryptCustom(v *Plaintext, baseLen int, isNTT bool) *Ciphertext {
	ct := NewCiphertextCustom(e.Parameters().Rank(), baseLen, isNTT)
	e.EncryptTo(ct, v, isNTT)
	return ct
}

// EncryptTo encrypts the message v to ctOut.
func (e *Encryptor) EncryptTo(ctOut *Ciphertext, v *Plaintext, isNTT bool) {
	e.heintEnc.EncryptTo(ctOut.Value, (*crt.Element)(v), isNTT)
	e.noise.EncryptTo(ctOut)
}

// EncryptElement encrypts the element e.
func (e *Encryptor) EncryptElement(eIn *rlwe.Element, isNTT bool) *Ciphertext {
	ct := NewCiphertext(e.Parameters(), isNTT)
	e.EncryptElementTo(ct, eIn, isNTT)
	return ct
}

// EncryptElementTo encrypts the element e to ctOut.
func (e *Encryptor) EncryptElementTo(ctOut *Ciphertext, eIn *rlwe.Element, isNTT bool) {
	e.heintEnc.EncryptElementTo(ctOut.Value, eIn, isNTT)
	e.noise.EncryptTo(ctOut)
}

// Decrypt decrypts the ciphertext ct.
func (e *Encryptor) Decrypt(ct *Ciphertext) *Plaintext {
	return (*Plaintext)(e.heintEnc.Decrypt(ct.Value))
}

// DecryptTo decrypts the ciphertext ct to vOut.
func (e *Encryptor) DecryptTo(vOut *Plaintext, ct *Ciphertext) {
	e.heintEnc.DecryptTo((*crt.Element)(vOut), ct.Value)
}

// Phase performs Phase(ct).
func (e *Encryptor) Phase(ct *Ciphertext) *rlwe.Element {
	return e.heintEnc.Phase(ct.Value)
}

// PhaseTo performs Phase(ct) and stores the result in pt.
func (e *Encryptor) PhaseTo(eOut *rlwe.Element, ct *Ciphertext) {
	e.heintEnc.PhaseTo((*rlwe.Element)(eOut), ct.Value)
}

// Noise returns the noise of the ciphertext.
func (e *Encryptor) Noise(ct *Ciphertext) *rlwe.Element {
	return e.heintEnc.Noise(ct.Value)
}

// NoiseTo stores the noise of the ciphertext in eOut.
func (e *Encryptor) NoiseTo(eOut *rlwe.Element, ct *Ciphertext) {
	e.heintEnc.NoiseTo(eOut, ct.Value)
}
