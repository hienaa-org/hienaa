package bfv

import (
	"math"
	"math/big"
	"sync"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// NoiseEstimator estimates the noise of the ciphertext.
type NoiseEstimator struct {
	params Parameters

	// TODO: do we need 'heavy' decomposer just to compute the decomposition length?
	dcmp rlwe.Decomposer

	expFac     *big.Float
	roundNoise *big.Float
	dcmpBound  *big.Float

	bigPool *sync.Pool
}

// NewNoiseEstimator creates a new [NoiseEstimator].
func NewNoiseEstimator(params Parameters) *NoiseEstimator {
	rlweParams := params.rlweParams

	prec := 0.0
	for i := range rlweParams.FullModulus() {
		prec += num.Log2(rlweParams.FullModulus()[i].Value())
	}
	prec = math.Ceil(prec)

	// Compute the expansion factor of the noise.
	expFac := big.NewFloat(float64(rlweParams.RingParams().ExpFactor()))
	switch keyType := rlweParams.SecretKeyParams().(type) {
	case *crt.TernarySamplerParameters:
		// TODO: Currently we simply use the expansion factor of a ternary vector with Hamming weight
		// as expFac * Hamming weight / rank. Is this the tight/correct bound?
		if keyType.HammingWeight > 0 {
			expFac.Mul(expFac, big.NewFloat(float64(keyType.HammingWeight)))
			expFac.Quo(expFac, big.NewFloat(float64(rlweParams.Rank())))
		}
	}

	// Compute the rounding noise.
	roundNoise := new(big.Float).SetPrec(uint(prec))
	switch params.estimType {
	case TypeVariance:
		roundNoise.SetInt64(1).Quo(roundNoise, big.NewFloat(12))
		expFac.Mul(expFac, rlweParams.SecretKeyParams().Variance())
	case TypeWorstCase:
		roundNoise.SetInt64(1).Quo(roundNoise, big.NewFloat(2))
		expFac.Mul(expFac, rlweParams.SecretKeyParams().Bound())
	default:
		panic("unsupported noise estimation type")
	}
	expFac.Add(expFac, big.NewFloat(1))
	roundNoise.Mul(roundNoise, expFac)

	// Compute the decomposition bound.
	dcmpBound := big.NewFloat(0)
	switch gadParams := rlweParams.GadgetParams().(type) {
	case rlwe.RNSGadgetParameters:
		dcmpBound.SetInt64(int64(gadParams.ChunkSize()))

		baseMod := rlweParams.BaseModulus()
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

		dcmpBound.Set(big.NewFloat(float64(maxBound.Int64())))
	case rlwe.DigitGadgetParameters:
		dcmpBound.SetInt(new(big.Int).Lsh(big.NewInt(1), uint(gadParams.LogDigitBase())))
	default:
		panic("unsupported gadget parameters")
	}

	switch params.estimType {
	case TypeVariance:
		dcmpBound.Mul(dcmpBound, dcmpBound)
		dcmpBound.Quo(dcmpBound, big.NewFloat(12))
	case TypeWorstCase:
		dcmpBound.Quo(dcmpBound, big.NewFloat(2))
	default:
		panic("unsupported noise estimation type")
	}

	return &NoiseEstimator{
		params: params,

		dcmp: rlwe.NewDecomposer(rlweParams),

		expFac:     expFac.SetInt64(int64(rlweParams.RingParams().ExpFactor())),
		roundNoise: roundNoise,
		dcmpBound:  dcmpBound,

		bigPool: &sync.Pool{
			New: func() any {
				return new(big.Float).SetPrec(uint(prec))
			},
		},
	}
}

// TODO: consider the noise from encoding when the plaintext modulus does not divide the ciphertext modulus.

// Encrypt returns the noise of the fresh ciphertext.
func (ne *NoiseEstimator) EncryptTo(cOut *Ciphertext) {
	switch ne.params.estimType {
	case TypeVariance:
		cOut.noise.Set(ne.params.rlweParams.NoiseParams().Variance())
	case TypeWorstCase:
		cOut.noise.Set(ne.params.rlweParams.NoiseParams().Bound())
	default:
		panic("unsupported noise estimation type")
	}
}

// ModSwitchTo switches the modulus of the ciphertext to the given length.
func (ne *NoiseEstimator) ModSwitchTo(cOut, ct *Ciphertext, l int) {
	ne.modSwitchTo(cOut.noise, ct, l)
}

// modSwitchTo stores the noise of the modulus switched ciphertext in noise.
func (ne *NoiseEstimator) modSwitchTo(noise *big.Float, ct *Ciphertext, l int) {
	modi := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(modi)
	noise.Set(ct.noise)

	switch {
	case ct.ModLen() < l:
		for i := ct.ModLen(); i < l; i++ {
			modi.SetUint64(ne.params.rlweParams.BaseModulus()[i].Value())
			noise.Quo(noise, modi)
		}
		noise.Add(noise, ne.roundNoise)
	case ct.ModLen() > l:
		for i := l; i < ct.ModLen(); i++ {
			modi.SetUint64(ne.params.rlweParams.BaseModulus()[i].Value())
			noise.Mul(noise, modi)
		}
	case ct.ModLen() == l:
		noise.Set(noise)
	}
}

// NegTo returns the noise of the negation of a ciphertext.
func (ne *NoiseEstimator) NegTo(cOut, ct *Ciphertext) {
	cOut.noise.Set(ct.noise)
}

// Add returns the noise of the sum of two ciphertexts.
func (ne *NoiseEstimator) AddTo(cOut, ct0, ct1 *Ciphertext) {
	tarLen := min(ct0.ModLen(), ct1.ModLen())

	noise0 := ne.bigPool.Get().(*big.Float)
	noise1 := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(noise0)
	defer ne.bigPool.Put(noise1)

	ne.modSwitchTo(noise0, ct0, tarLen)
	ne.modSwitchTo(noise1, ct1, tarLen)

	cOut.noise.Add(noise0, noise1)
}

// AddPlainTo returns the noise of the sum of a ciphertext and a plaintext.
func (ne *NoiseEstimator) AddPlainTo(cOut, ct *Ciphertext, pt *Plaintext) {
	cOut.noise.Set(ct.noise)
}

// AddElementTo returns the noise of the sum of a ciphertext and an element.
func (ne *NoiseEstimator) AddElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	cOut.noise.Set(ct.noise)
}

// SubTo returns the noise of the difference of two ciphertexts.
func (ne *NoiseEstimator) SubTo(cOut, ct0, ct1 *Ciphertext) {
	ne.AddTo(cOut, ct0, ct1)
}

// SubPlainTo returns the noise of the difference of a ciphertext and a plaintext.
func (ne *NoiseEstimator) SubPlainTo(cOut, ct *Ciphertext, pt *Plaintext) {
	ne.AddPlainTo(cOut, ct, pt)
}

// SubElementTo returns the noise of the difference of a ciphertext and an element.
func (ne *NoiseEstimator) SubElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	ne.AddElementTo(cOut, ct, e)
}

// Mul returns the noise of the product of two ciphertexts.
func (ne *NoiseEstimator) MulTo(ctOut, ct0, ct1 *Ciphertext) {
	tarLen := min(ct0.ModLen(), ct1.ModLen())

	noise0 := ne.bigPool.Get().(*big.Float)
	noise1 := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(noise0)
	defer ne.bigPool.Put(noise1)

	ne.modSwitchTo(noise0, ct0, tarLen)
	ne.modSwitchTo(noise1, ct1, tarLen)

	// Given phases q/t*m + e + qI and q/t*m' + e' + qI' of lifted ciphertexts,
	// the output noise is given by q/t*mm' + (me'+m'e) + t*(eI' + e'I) + t/q * ee' + e_rnd + keyswitching noise.
	tmp := ne.bigPool.Get().(*big.Float)
	msgMod := ne.bigPool.Get().(*big.Float)
	baseMod := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(tmp)
	defer ne.bigPool.Put(msgMod)
	defer ne.bigPool.Put(baseMod)
	msgMod.SetUint64(ne.params.MessageModulus().Value())
	baseMod.SetInt64(1)
	for i := 0; i < tarLen; i++ {
		tmp.SetUint64(ne.params.rlweParams.BaseModulus()[i].Value())
		baseMod.Mul(baseMod, tmp)
	}

	// key-switching noise.
	// TODO: rounding in tensor form adds extra noise.
	ne.gadgetProdTo(ctOut)

	// t/q * ee'
	tmp.Quo(msgMod, baseMod)
	tmp.Mul(tmp, noise0)
	tmp.Mul(tmp, noise1)
	tmp.Mul(tmp, ne.expFac)
	ctOut.noise.Add(ctOut.noise, tmp)

	// t*(eI' + e'I)
	tmp.Add(noise0, noise1)
	tmp.Mul(tmp, ne.roundNoise) // I has the same bound as the rouding noise.
	tmp.Mul(tmp, ne.expFac)
	tmp.Mul(tmp, msgMod)
	ctOut.noise.Add(ctOut.noise, tmp)

	// me' + m'e
	halfMsgMod := msgMod.Quo(msgMod, tmp.SetInt64(2))
	tmp.Add(noise0, noise1)
	tmp.Mul(tmp, ne.expFac)
	tmp.Mul(tmp, halfMsgMod)
	ctOut.noise.Add(ctOut.noise, tmp)

	// e_rnd
	ctOut.noise.Add(ctOut.noise, ne.roundNoise)
}

// MulPlainTo returns the noise of the product of a ciphertext and a plaintext.
func (ne *NoiseEstimator) MulPlainTo(cOut, ct *Ciphertext, pt *Plaintext) {
	// Compute the tight bound of the noise.
	max := uint64(0)
	halfMsgMod := ne.params.messageModulus.Value() >> 1
	for i := 0; i < pt.Rank(); i++ {
		if pt.Coeffs[0][i] > halfMsgMod && pt.Coeffs[0][i]-halfMsgMod > max {
			max = pt.Coeffs[0][i] - halfMsgMod
		} else if pt.Coeffs[0][i] <= halfMsgMod && pt.Coeffs[0][i] > max {
			max = pt.Coeffs[0][i]
		}
	}
	cOut.noise.SetInt64(int64(max))
	cOut.noise.Mul(cOut.noise, ne.roundNoise)
	if pt.Rank() > 1 {
		cOut.noise.Mul(cOut.noise, ne.expFac)
	}
}

// MulElementTo returns the noise of the product of a ciphertext and an element.
func (ne *NoiseEstimator) MulElementTo(cOut, ct *Ciphertext, e *rlwe.Element) {
	// When we multiply a ciphertext and an element, we cannot compute the tight bound of the noise
	// without expensive operations, such as basis embedding or inverse NTT.
	// Therefore, we use the half of the message modulus as the noise bound.
	halfMsgMod := ne.params.messageModulus.Value() >> 1
	cOut.noise.SetInt64(int64(halfMsgMod))
	cOut.noise.Mul(cOut.noise, ne.roundNoise)
	if e.Rank() > 1 {
		cOut.noise.Mul(cOut.noise, ne.expFac)
	}
}

// gadgetProdTo computes the noise of the gadget product of a ciphertext.
func (ne *NoiseEstimator) gadgetProdTo(cOut *Ciphertext) {
	modLen := cOut.ModLen()

	decLen := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(decLen)
	decLen.SetInt64(int64(ne.dcmp.DecomposeLen(modLen)))

	switch ne.params.estimType {
	case TypeVariance:
		cOut.noise.Set(ne.params.rlweParams.NoiseParams().Variance())
	case TypeWorstCase:
		cOut.noise.Set(ne.params.rlweParams.NoiseParams().Bound())
	default:
		panic("unsupported noise estimation type")
	}

	cOut.noise.Mul(cOut.noise, ne.expFac)
	cOut.noise.Mul(cOut.noise, ne.dcmpBound)
	cOut.noise.Mul(cOut.noise, decLen)

	if ne.params.rlweParams.HasAuxModulus() {
		auxLen := ne.dcmp.AuxModLen(modLen)

		auxMod := ne.bigPool.Get().(*big.Float)
		modi := ne.bigPool.Get().(*big.Float)
		defer ne.bigPool.Put(auxMod)
		defer ne.bigPool.Put(modi)
		auxMod.SetInt64(1)
		for i := 0; i < auxLen; i++ {
			modi.SetUint64(ne.params.rlweParams.AuxModulus()[i].Value())
			auxMod.Mul(auxMod, modi)
		}
		cOut.noise.Quo(cOut.noise, auxMod)
		cOut.noise.Add(cOut.noise, ne.roundNoise)
	}
}

// KeySwitchTo returns the noise of the key switched ciphertext.
func (ne *NoiseEstimator) KeySwitchTo(cOut, ct *Ciphertext) {
	ctNoise := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(ctNoise)
	ctNoise.Set(ct.noise)

	ne.gadgetProdTo(cOut)
	cOut.noise.Add(cOut.noise, ctNoise)
}

// AutTo returns the noise of the automorphism of a ciphertext.
func (ne *NoiseEstimator) AutTo(cOut, ct *Ciphertext) {
	ctNoise := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(ctNoise)
	ctNoise.Set(ct.noise)

	ne.gadgetProdTo(cOut)
	cOut.noise.Add(cOut.noise, ctNoise)
}

// ExtProdTo returns the noise of the external product of a ciphertext.
func (ne *NoiseEstimator) ExtProdTo(cOut *Ciphertext, ct *Ciphertext) {
	ctNoise := ne.bigPool.Get().(*big.Float)
	two := ne.bigPool.Get().(*big.Float)
	defer ne.bigPool.Put(ctNoise)
	defer ne.bigPool.Put(two)
	ctNoise.Set(ct.noise)
	two.SetInt64(2)

	ne.gadgetProdTo(cOut)
	cOut.noise.Mul(cOut.noise, two)
	cOut.noise.Add(cOut.noise, ctNoise)
}
