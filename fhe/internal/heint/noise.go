package heint

import (
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

type EstimType int

const (
	VarianceType EstimType = iota
	WorstCaseType
)

// NoiseEstimator estimates the noise of the ciphertext.
type NoiseEstimator struct {
	params    rlwe.Parameters
	msgMod    *num.Modulus
	estimType EstimType

	// TODO: do we need 'heavy' decomposer just to compute the decomposition length?
	dcmp      rlwe.Decomposer
	keyExpFac float64
	dcmpBound float64
}

// NewNoiseEstimator creates a new [NoiseEstimator].
func NewNoiseEstimator(params rlwe.Parameters, msgMod *num.Modulus, estimType EstimType) *NoiseEstimator {
	// Compute the expansion factor of the noise.
	keyExpFac := float64(params.RingParams().ExpFactor())
	switch keyType := params.SecretKeyParams().(type) {
	case crt.TernarySamplerParameters:
		// TODO: Currently we simply use the expansion factor of a ternary vector with Hamming weight
		// as expFac * Hamming weight / rank. Is this the tight/correct bound?
		if keyType.HammingWeight > 0 {
			keyExpFac *= float64(keyType.HammingWeight) / float64(params.Rank())
		}
	default:
		switch estimType {
		case VarianceType:
			keyExpFac *= params.SecretKeyParams().Variance()
		case WorstCaseType:
			keyExpFac *= params.SecretKeyParams().Bound()
		default:
			panic("unsupported noise estimation type")
		}
	}

	// Compute the decomposition bound.
	var dcmpBound float64
	switch gadParams := params.GadgetParams().(type) {
	case rlwe.RNSGadgetParameters:
		dcmpBound = float64(gadParams.ChunkSize())

		baseMod := params.BaseModulus()
		chunkSize := gadParams.ChunkSize()

		maxBound := new(big.Int)
		chunkBound := new(big.Int)
		modi := new(big.Int)
		for i := 0; i < len(baseMod)/chunkSize; i++ {
			chunkBound.SetInt64(1)
			for j := 0; j < chunkSize; j++ {
				modi.SetUint64(baseMod[i*chunkSize+j].Value())
				chunkBound.Mul(chunkBound, modi)
			}

			if maxBound.Cmp(chunkBound) < 0 {
				maxBound.Set(chunkBound)
			}
		}

		dcmpBound, _ = maxBound.Float64()
	case rlwe.DigitGadgetParameters:
		dcmpBound = math.Exp2(float64(gadParams.LogDigitBase()))
	default:
		panic("unsupported gadget parameters")
	}

	switch estimType {
	case VarianceType:
		dcmpBound = dcmpBound * dcmpBound / 12
	case WorstCaseType:
		dcmpBound /= 2
	default:
		panic("unsupported noise estimation type")
	}

	return &NoiseEstimator{
		params:    params,
		msgMod:    msgMod,
		estimType: estimType,

		dcmp:      rlwe.NewDecomposer(params),
		keyExpFac: keyExpFac,
		dcmpBound: dcmpBound,
	}
}

// KeyExpansionFactor returns the key expansion factor.
func (ne *NoiseEstimator) KeyExpansionFactor() float64 {
	return ne.keyExpFac
}

// RoundNoise returns the rounding noise.
func (ne *NoiseEstimator) RoundNoise() float64 {
	switch ne.estimType {
	case VarianceType:
		return 1 + ne.keyExpFac/12
	case WorstCaseType:
		return 1 + ne.keyExpFac/2
	default:
		panic("unsupported noise estimation type")
	}
}

// TODO: consider the noise from encoding when the plaintext modulus does not divide the ciphertext modulus.

// Encrypt returns the noise of the fresh ciphertext.
func (ne *NoiseEstimator) Encrypt() float64 {
	switch ne.estimType {
	case VarianceType:
		return ne.params.NoiseParams().Variance()
	case WorstCaseType:
		return ne.params.NoiseParams().Bound()
	default:
		panic("unsupported noise estimation type")
	}
}

// RescaleTo rescales the noise of the ciphertext to the target modulus.
func (ne *NoiseEstimator) RescaleTo(noise float64, modLen int) float64 {
	var scale float64
	switch ne.estimType {
	case VarianceType:
		scale = math.Ceil(math.Sqrt(noise / ne.RoundNoise()))
	case WorstCaseType:
		scale = math.Ceil(noise / ne.RoundNoise())
	}

	divMod := float64(ne.params.BaseModulus()[modLen-1].Value())
	if divMod < scale {
		switch ne.estimType {
		case VarianceType:
			return noise/(divMod*divMod) + ne.RoundNoise()
		case WorstCaseType:
			return noise/divMod + ne.RoundNoise()
		default:
			panic("unsupported noise estimation type")
		}
	}

	return noise
}

// ModSwitchTo switches the modulus of the ciphertext to the given length.
func (ne *NoiseEstimator) ModSwitch(noise float64, inLen, outLen int) float64 {
	res := noise
	switch {
	case inLen < outLen:
		for i := inLen; i < outLen; i++ {
			res *= float64(ne.params.BaseModulus()[i].Value())
		}
	case inLen > outLen:
		for i := outLen; i < inLen; i++ {
			res /= float64(ne.params.BaseModulus()[i].Value())
		}
		res += ne.RoundNoise()
	}

	return res
}

// NegTo returns the noise of the negation of a ciphertext.
func (ne *NoiseEstimator) Neg(noise float64) float64 {
	return noise
}

// Add returns the noise of the sum of two ciphertexts.
func (ne *NoiseEstimator) Add(noise0, noise1 float64, len0, len1 int) float64 {
	tarLen := min(len0, len1)

	return ne.ModSwitch(noise0, len0, tarLen) + ne.ModSwitch(noise1, len1, tarLen)
}

// AddPlainTo returns the noise of the sum of a ciphertext and a plaintext.
func (ne *NoiseEstimator) AddPlain(noise float64, pt []uint64) float64 {
	return noise
}

// AddElementTo returns the noise of the sum of a ciphertext and an element.
func (ne *NoiseEstimator) AddElement(noise float64, e *rlwe.Element) float64 {
	return noise
}

// SubTo returns the noise of the difference of two ciphertexts.
func (ne *NoiseEstimator) Sub(noise0, noise1 float64, len0, len1 int) float64 {
	return ne.Add(noise0, noise1, len0, len1)
}

// SubPlainTo returns the noise of the difference of a ciphertext and a plaintext.
func (ne *NoiseEstimator) SubPlain(noise float64, pt []uint64) float64 {
	return noise
}

// SubElementTo returns the noise of the difference of a ciphertext and an element.
func (ne *NoiseEstimator) SubElement(noise float64, e *rlwe.Element) float64 {
	return noise
}

// MulPlainTo returns the noise of the product of a ciphertext and a plaintext.
func (ne *NoiseEstimator) MulPlain(noise float64, pt []uint64) float64 {
	// Compute the tight bound of the noise.
	max := uint64(0)
	halfMsgMod := ne.msgMod.Value() >> 1
	for i := 0; i < len(pt); i++ {
		if pt[i] > halfMsgMod && pt[i]-halfMsgMod > max {
			max = pt[i] - halfMsgMod
		} else if pt[i] <= halfMsgMod && pt[i] > max {
			max = pt[i]
		}
	}

	if len(pt) == 1 {
		return noise * float64(max)
	} else {
		return noise * float64(max) * float64(ne.params.RingParams().ExpFactor())
	}
}

// MulElementTo returns the noise of the product of a ciphertext and an element.
func (ne *NoiseEstimator) MulElement(noise float64, e *rlwe.Element) float64 {
	// When we multiply a ciphertext and an element, we cannot compute the tight bound of the noise
	// without expensive operations, such as basis embedding or inverse NTT.
	// Therefore, we use the half of the message modulus as the noise bound.
	if e.Type() == crt.TypeScalar {
		return noise * float64(ne.msgMod.Value()) / 2
	} else {
		return noise * float64(ne.msgMod.Value()) / 2 * float64(ne.params.RingParams().ExpFactor())
	}
}

// GadgetProd computes the noise of the gadget product of a ciphertext.
func (ne *NoiseEstimator) GadgetProd(modLen int) float64 {
	res := float64(0)

	switch ne.estimType {
	case VarianceType:
		res = ne.params.NoiseParams().Variance()
	case WorstCaseType:
		res = ne.params.NoiseParams().Bound()
	default:
		panic("unsupported noise estimation type")
	}

	res *= float64(ne.params.RingParams().ExpFactor()) * ne.dcmpBound * float64(ne.dcmp.DecomposeLen(modLen))

	if ne.params.HasAuxModulus() {
		auxLen := ne.dcmp.AuxModLen(modLen)
		for i := 0; i < auxLen; i++ {
			res /= float64(ne.params.AuxModulus()[i].Value())
		}
		res += ne.RoundNoise()
	}

	return res
}

// KeySwitchTo returns the noise of the key switched ciphertext.
func (ne *NoiseEstimator) KeySwitch(noise float64, modLen int) float64 {
	return noise + ne.GadgetProd(modLen)
}

// AutTo returns the noise of the automorphism of a ciphertext.
func (ne *NoiseEstimator) Aut(noise float64, modLen int) float64 {
	return noise + ne.GadgetProd(modLen)
}
