package bfv

import (
	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

const (
	VarianceType heint.EstimType = iota
	WorstCaseType
)

// NoiseEstimator estimates the noise of the ciphertext.
type NoiseEstimator struct {
	params    rlwe.Parameters
	msgMod    *num.Modulus
	estimType heint.EstimType

	noise *heint.NoiseEstimator
}

// NewNoiseEstimator creates a new [NoiseEstimator].
func NewNoiseEstimator(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *NoiseEstimator {
	return &NoiseEstimator{
		params:    params,
		msgMod:    msgMod,
		estimType: estimType,

		noise: heint.NewNoiseEstimator(params, msgMod, estimType),
	}
}

// TODO: consider the noise from encoding when the plaintext modulus does not divide the ciphertext modulus.

// Encrypt returns the noise of the fresh ciphertext.
func (ne *NoiseEstimator) EncryptTo(cOut *Ciphertext) {
	cOut.noise = ne.noise.Encrypt()
}

// ModSwitchTo switches the modulus of the ciphertext to the given length.
func (ne *NoiseEstimator) ModSwitchTo(cOut, ct *Ciphertext, l int) {
	cOut.noise = ne.noise.ModSwitch(cOut.noise, ct.ModLen(), l)
}

// FwdNTTTo computes ctOut = FwdNTT(ct).
func (ne *NoiseEstimator) FwdNTTTo(cOut, ct *Ciphertext) {
	cOut.noise = ct.noise
}

// InvNTTTo computes ctOut = InvNTT(ct).
func (ne *NoiseEstimator) InvNTTTo(cOut, ct *Ciphertext) {
	cOut.noise = ct.noise
}

// NegTo returns the noise of the negation of a ciphertext.
func (ne *NoiseEstimator) NegTo(cOut, ct *Ciphertext) {
	cOut.noise = ne.noise.Neg(ct.noise)
}

// Add returns the noise of the sum of two ciphertexts.
func (ne *NoiseEstimator) AddTo(cOut, ct0, ct1 *Ciphertext) {
	cOut.noise = ne.noise.Add(ct0.noise, ct1.noise, ct0.ModLen(), ct1.ModLen())
}

// AddPlainTo returns the noise of the sum of a ciphertext and a plaintext.
func (ne *NoiseEstimator) AddPlainTo(cOut, ct *Ciphertext, pt Plaintext) {
	cOut.noise = ne.noise.AddPlain(ct.noise, pt)
}

// AddElementTo returns the noise of the sum of a ciphertext and an element.
func (ne *NoiseEstimator) AddElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	cOut.noise = ne.noise.AddElement(ct.noise, e)
}

// SubTo returns the noise of the difference of two ciphertexts.
func (ne *NoiseEstimator) SubTo(cOut, ct0, ct1 *Ciphertext) {
	cOut.noise = ne.noise.Sub(ct0.noise, ct1.noise, ct0.ModLen(), ct1.ModLen())
}

// SubPlainTo returns the noise of the difference of a ciphertext and a plaintext.
func (ne *NoiseEstimator) SubPlainTo(cOut, ct *Ciphertext, pt Plaintext) {
	cOut.noise = ne.noise.SubPlain(ct.noise, pt)
}

// SubElementTo returns the noise of the difference of a ciphertext and an element.
func (ne *NoiseEstimator) SubElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	cOut.noise = ne.noise.SubElement(ct.noise, e)
}

// Mul returns the noise of the product of two ciphertexts.
func (ne *NoiseEstimator) MulTo(ctOut, ct0, ct1 *Ciphertext) {
	// TODO: noise estimation for the improved BFV multiplication.

	tarLen := min(ct0.ModLen(), ct1.ModLen())

	noise0 := ne.noise.ModSwitch(ct0.noise, ct0.ModLen(), tarLen)
	noise1 := ne.noise.ModSwitch(ct1.noise, ct1.ModLen(), tarLen)

	// Given phases q/t*m + e + qI and q/t*m' + e' + qI' of lifted ciphertexts,
	// the output noise is given by q/t*mm' + (me'+m'e) + t*(eI' + e'I) + t/q * ee' + e_rnd + keyswitching noise.

	// e_rnd
	keyExpFac := ne.noise.KeyExpansionFactor()
	if ne.estimType == VarianceType {
		ctOut.noise = (1 + keyExpFac + keyExpFac*keyExpFac) / 12
	} else {
		ctOut.noise = (1 + keyExpFac + keyExpFac*keyExpFac) / 2
	}

	// key-switching noise.
	ctOut.noise += ne.noise.KeySwitch(ctOut.noise, tarLen)

	// t/q * ee'
	msgMod := float64(ne.msgMod.Value())
	expFac := float64(ne.params.RingParams().ExpFactor())
	baseMod := float64(1)
	for i := 0; i < tarLen; i++ {
		baseMod *= float64(ne.params.BaseModulus()[i].Value())
	}
	ctOut.noise += msgMod / baseMod * noise0 * noise1 * expFac

	// t*(eI' + e'I)
	ctOut.noise += msgMod * expFac * (noise0 + noise1) * ne.noise.RoundNoise()

	// me' + m'e
	ctOut.noise += msgMod / 2 * (noise0 + noise1) * expFac
}

// MulPlainTo returns the noise of the product of a ciphertext and a plaintext.
func (ne *NoiseEstimator) MulPlainTo(cOut, ct *Ciphertext, pt Plaintext) {
	cOut.noise = ne.noise.MulPlain(ct.noise, pt)
}

// MulElementTo returns the noise of the product of a ciphertext and an element.
func (ne *NoiseEstimator) MulElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	// When we multiply a ciphertext and an element, we cannot compute the tight bound of the noise
	// without expensive operations, such as basis embedding or inverse NTT.
	// Therefore, we use the half of the message modulus as the noise bound.
	cOut.noise = ne.noise.MulElement(ct.noise, e)
}

// KeySwitchTo returns the noise of the key switched ciphertext.
func (ne *NoiseEstimator) KeySwitchTo(cOut, ct *Ciphertext) {
	cOut.noise = ne.noise.KeySwitch(ct.noise, ct.ModLen())
}

// AutTo returns the noise of the automorphism of a ciphertext.
func (ne *NoiseEstimator) AutTo(cOut, ct *Ciphertext) {
	cOut.noise = ne.noise.Aut(ct.noise, ct.ModLen())
}
