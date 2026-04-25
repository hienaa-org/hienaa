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

// ModSwitchTo returns the noise of the modulus switched ciphertext.
func (ne *NoiseEstimator) ModSwitchTo(ctOut, ct *Ciphertext, l int) {
	ctOut.noise = ne.noise.ModSwitch(ct.noise, ct.ModLen(), l)
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

// getAuxMod gets the auxiliary modulus for the multiplication.
func (ne *NoiseEstimator) getAuxMod(ct0, ct1 *Ciphertext) (int, int, *num.Modulus) {
	// Force ct0 to have the larger noise.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	curLen := ct0.ModLen()
	tarLen := curLen

	var scale float64
	switch ne.estimType {
	case heint.VarianceType:
		scale = math.Max(1, math.Sqrt(ct0.noise/ne.noise.RoundNoise()))
	case heint.WorstCaseType:
		scale = math.Max(1, ct0.noise/ne.noise.RoundNoise())
	}

	for scale > math.Exp2(num.MaxModulusBits) && tarLen > 0 {
		tarLen--
		if tarLen == 0 {
			panic("ciphertext noise is too large to perform multiplication")
		}
		scale /= float64(ne.params.BaseModulus()[tarLen].Value())
	}

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
	rem := num.Inv(remInv, ne.msgMod)

	divMod := ne.params.BaseModulus()[auxIdx].Value()
	msgMod := ne.msgMod.Value()
	auxMod := uint64(math.Floor(float64(divMod)/scale/float64(msgMod)))*msgMod + rem

	if auxMod == 1 {
		return tarLen - 1, tarLen - 1, nil
	} else {
		return tarLen, auxIdx, num.NewModulus(auxMod)
	}
}

// MulTo returns the noise of the product of two ciphertexts.
func (ne *NoiseEstimator) MulTo(ctOut, ct0, ct1 *Ciphertext) {
	_, auxIdx, auxMod := ne.getAuxMod(ct0, ct1)

	var scale float64
	if auxMod != nil {
		scale = float64(ne.params.BaseModulus()[auxIdx].Value()) / float64(auxMod.Value())
	} else {
		scale = 1
	}

	switch ne.estimType {
	case heint.VarianceType:
		noise0 := ct0.noise/scale/scale + ne.noise.RoundNoise()
		noise1 := ct1.noise/scale/scale + ne.noise.RoundNoise()
		msgMod := float64(ne.msgMod.Value())
		expFac := float64(ne.params.RingParams().ExpandFactor())

		mulNoise := noise0 * noise1 * msgMod * msgMod * expFac
		ctOut.noise = mulNoise*scale*scale + ne.noise.RoundNoise()
	case heint.WorstCaseType:
		noise0 := ct0.noise/scale + ne.noise.RoundNoise()
		noise1 := ct1.noise/scale + ne.noise.RoundNoise()
		msgMod := float64(ne.msgMod.Value())
		expFac := float64(ne.params.RingParams().ExpandFactor())

		mulNoise := noise0 * noise1 * msgMod * expFac
		ctOut.noise = mulNoise*scale + ne.noise.RoundNoise()
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
