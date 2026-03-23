package bgv

import (
	"math"

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

// EncryptTo returns the noise of the fresh ciphertext.
func (ne *NoiseEstimator) EncryptTo(ctOut *Ciphertext) {
	ctOut.noise = ne.noise.Encrypt()
}

// RescaleTo rescales the noise of the ciphertext to the target modulus.
func (ne *NoiseEstimator) RescaleTo(ctOut, ct *Ciphertext) {
	ctOut.noise = ne.noise.RescaleTo(ct.noise, ct.ModLen())
}

// ModSwitchTo returns the noise of the modulus switched ciphertext.
func (ne *NoiseEstimator) ModSwitchTo(ctOut, ct *Ciphertext, l int) {
	ctOut.noise = ne.noise.ModSwitch(ctOut.noise, ct.ModLen(), l)
}

// FwdNTTTo returns the noise of the forward NTT transformed ciphertext.
func (ne *NoiseEstimator) FwdNTTTo(ctOut, ct *Ciphertext) {
	ctOut.noise = ct.noise
}

// InvNTTTo returns the noise of the inverse NTT transformed ciphertext.
func (ne *NoiseEstimator) InvNTTTo(ctOut, ct *Ciphertext) {
	ctOut.noise = ct.noise
}

// NegTo returns the noise of the negation of a ciphertext.
func (ne *NoiseEstimator) NegTo(ctOut, ct *Ciphertext) {
	ctOut.noise = ne.noise.Neg(ct.noise)
}

// Add returns the noise of the sum of two ciphertexts.
func (ne *NoiseEstimator) AddTo(ctOut, ct0, ct1 *Ciphertext) {
	ctOut.noise = ne.noise.Add(ct0.noise, ct1.noise, ct0.ModLen(), ct1.ModLen())
}

// AddPlainTo returns the noise of the sum of a ciphertext and a plaintext.
func (ne *NoiseEstimator) AddPlainTo(ctOut, ct *Ciphertext, pt Plaintext) {
	ctOut.noise = ne.noise.AddPlain(ct.noise, pt)
}

// AddElementTo returns the noise of the sum of a ciphertext and an element.
func (ne *NoiseEstimator) AddElementTo(ctOut, ct *Ciphertext, e *rlwe.Element) {
	ctOut.noise = ne.noise.AddElement(ct.noise, e)
}

// SubTo returns the noise of the difference of two ciphertexts.
func (ne *NoiseEstimator) SubTo(ctOut, ct0, ct1 *Ciphertext) {
	ctOut.noise = ne.noise.Sub(ct0.noise, ct1.noise, ct0.ModLen(), ct1.ModLen())
}

// SubPlainTo returns the noise of the difference of a ciphertext and a plaintext.
func (ne *NoiseEstimator) SubPlainTo(ctOut, ct *Ciphertext, pt Plaintext) {
	ctOut.noise = ne.noise.SubPlain(ct.noise, pt)
}

// SubElementTo returns the noise of the difference of a ciphertext and an element.
func (ne *NoiseEstimator) SubElementTo(ctOut, ct *Ciphertext, e *rlwe.Element) {
	ctOut.noise = ne.noise.SubElement(ct.noise, e)
}

// MulTo returns the noise of the product of two ciphertexts.
func (ne *NoiseEstimator) MulTo(ctOut, ct0, ct1 *Ciphertext) {
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	tarLen := ct0.ModLen()
	auxIdx := 0
	for auxIdx < tarLen-1 {
		if num.GCD(ne.params.BaseModulus()[auxIdx].Value(), ne.msgMod.Value()) != 1 {
			break
		}
		auxIdx++
	}

	remInv := uint64(1)
	for i := 0; i < tarLen; i++ {
		if i != auxIdx {
			remInv = num.Mul(remInv, ne.params.BaseModulus()[i].Value(), ne.msgMod)
		}
	}
	rem := float64(num.Inv(remInv, ne.msgMod))

	var scale float64
	switch ne.estimType {
	case heint.VarianceType:
		scale = math.Ceil(math.Sqrt(ct0.noise / ne.noise.RoundNoise()))
	case heint.WorstCaseType:
		scale = math.Ceil(ct0.noise / ne.noise.RoundNoise())
	}
	divMod := float64(ne.params.BaseModulus()[auxIdx].Value())
	msgMod := float64(ne.msgMod.Value())
	auxMod := math.Floor(divMod/scale*msgMod)*scale*msgMod + rem

	switch ne.estimType {
	case heint.VarianceType:
		noise0 := ct0.noise/(divMod/auxMod)/(divMod/auxMod) + ne.noise.RoundNoise()
		noise1 := ct1.noise/(divMod/auxMod)/(divMod/auxMod) + ne.noise.RoundNoise()
		msgMod := float64(ne.msgMod.Value())
		expFac := float64(ne.params.RingParams().ExpFactor())

		mulNoise := noise0 * noise1 * msgMod * msgMod * expFac
		ctOut.noise = mulNoise*(divMod/auxMod)*(divMod/auxMod) + ne.noise.RoundNoise()
	case heint.WorstCaseType:
		noise0 := ct0.noise/(divMod/auxMod) + ne.noise.RoundNoise()
		noise1 := ct1.noise/(divMod/auxMod) + ne.noise.RoundNoise()
		msgMod := float64(ne.msgMod.Value())
		expFac := float64(ne.params.RingParams().ExpFactor())

		mulNoise := noise0 * noise1 * msgMod * expFac
		ctOut.noise = mulNoise*(divMod/auxMod) + ne.noise.RoundNoise()
	}
}

// MulPlainTo returns the noise of the product of a ciphertext and a plaintext.
func (ne *NoiseEstimator) MulPlainTo(ctOut, ct *Ciphertext, pt Plaintext) {
	ctOut.noise = ne.noise.MulPlain(ct.noise, pt)
}

// MulElementTo returns the noise of the product of a ciphertext and an element.
func (ne *NoiseEstimator) MulElementTo(ctOut, ct *Ciphertext, e *rlwe.Element) {
	ctOut.noise = ne.noise.MulElement(ct.noise, e)
}

// KeySwitchTo returns the noise of the key switched ciphertext.
func (ne *NoiseEstimator) KeySwitchTo(ctOut, ct *Ciphertext) {
	ctOut.noise = ne.noise.KeySwitch(ct.noise, ct.ModLen())
}

// AutTo returns the noise of the automorphism of a ciphertext.
func (ne *NoiseEstimator) AutTo(ctOut, ct *Ciphertext) {
	ctOut.noise = ne.noise.Aut(ct.noise, ct.ModLen())
}
