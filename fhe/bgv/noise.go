package bgv

import (
	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
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

// EncryptTo returns the noise of the fresh ciphertext.
func (ne *NoiseEstimator) EncryptTo(cOut *Ciphertext) {
	cOut.noise = ne.noise.Encrypt()
}

// ModSwitchTo returns the noise of the modulus switched ciphertext.
func (ne *NoiseEstimator) ModSwitchTo(cOut, ct *Ciphertext, l int) {
	cOut.noise = ne.noise.ModSwitch(cOut.noise, ct.ModLen(), l)
}

// FwdNTTTo returns the noise of the forward NTT transformed ciphertext.
func (ne *NoiseEstimator) FwdNTTTo(cOut, ct *Ciphertext) {
	cOut.noise = ct.noise
}

// InvNTTTo returns the noise of the inverse NTT transformed ciphertext.
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
func (ne *NoiseEstimator) AddPlainTo(cOut, ct *Ciphertext, pt *Plaintext) {
	cOut.noise = ne.noise.AddPlain(ct.noise, (*crt.Element)(pt))
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
func (ne *NoiseEstimator) SubPlainTo(cOut, ct *Ciphertext, pt *Plaintext) {
	cOut.noise = ne.noise.SubPlain(ct.noise, (*crt.Element)(pt))
}

// SubElementTo returns the noise of the difference of a ciphertext and an element.
func (ne *NoiseEstimator) SubElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	cOut.noise = ne.noise.SubElement(ct.noise, e)
}

// MulTo returns the noise of the product of two ciphertexts.
func (ne *NoiseEstimator) MulTo(cOut, ct0, ct1 *Ciphertext) {

}

// MulPlainTo returns the noise of the product of a ciphertext and a plaintext.
func (ne *NoiseEstimator) MulPlainTo(cOut, ct *Ciphertext, pt *Plaintext) {
	cOut.noise = ne.noise.MulPlain(ct.noise, (*crt.Element)(pt))
}

// MulElementTo returns the noise of the product of a ciphertext and an element.
func (ne *NoiseEstimator) MulElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
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
