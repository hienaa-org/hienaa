package bfv

import (
	"cmp"
	"math"
	"slices"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

const (
	VarianceType heint.EstimType = iota
	WorstCaseType
)

// NoiseEstimator estimates the noise of the ciphertext.
type NoiseEstimator struct {
	params    rlwe.Parameters
	ambMod    []*num.Modulus
	msgMod    *num.Modulus
	estimType heint.EstimType

	noise *heint.NoiseEstimator
}

// NewNoiseEstimator creates a new [NoiseEstimator].
func NewNoiseEstimator(params rlwe.Parameters, msgMod *num.Modulus, estimType heint.EstimType) *NoiseEstimator {
	ringParams := params.RingParams()
	extraBits := float64(0)
	for _, mod := range params.BaseModulus() {
		extraBits += num.Log2(mod.Value())
	}
	extraBits = extraBits + num.Log2(ringParams.ExpandFactor())
	extraLen := int(math.Ceil(extraBits / num.MaxModulusBits))
	baseMod := params.BaseModulus()

	var extraMod []*num.Modulus
	var bitlen float64
	for {
		extraMod = make([]*num.Modulus, extraLen)

		gap, err := dft.NTTPrimeGap(ringParams)
		if err != nil {
			panic(err)
		}

		bitlen = extraBits
		start := (uint64(math.Floor(num.MaxModulus/float64(gap))))*gap + 1
		if start > num.MaxModulus {
			for start > num.MaxModulus {
				start -= gap
			}
		}
		prime := num.MustPrevPrime(start, gap)

		cnt := 0
		for cnt < extraLen {
			primemod := num.NewModulus(prime)
			for {
				_, ok := slices.BinarySearchFunc(baseMod, primemod, func(a, b *num.Modulus) int {
					return cmp.Compare(a.Value(), b.Value())
				})

				if !ok {
					break
				} else {
					prime = num.MustPrevPrime(prime, gap)
					primemod = num.NewModulus(prime)
				}
			}
			extraMod[cnt] = primemod
			bitlen -= num.Log2(primemod.Value())
			prime = num.MustPrevPrime(prime, gap)
			cnt++
		}

		if bitlen > 0 {
			extraLen++
		} else {
			break
		}
	}

	ambMod := append(baseMod, extraMod...)

	return &NoiseEstimator{
		params:    params,
		ambMod:    ambMod,
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
	cOut.noise = ne.noise.ModSwitch(ct.noise, ct.ModLen(), l)
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

// getAuxLen returns the length of the auxiliary modulus for the multiplication of two ciphertexts.
func (ne *NoiseEstimator) getAmbLen(ct0, ct1 *Ciphertext) int {
	// Force ct0 to have the larger noise.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	tarLen := ct0.ModLen()
	tarBits := float64(0)
	for i := 0; i < tarLen; i++ {
		tarBits += num.Log2(ne.ambMod[i].Value())
	}

	switch ne.estimType {
	case heint.VarianceType:
		tarBits -= num.Log2(ct0.noise/ne.noise.RoundNoise()) / 2
	case heint.WorstCaseType:
		tarBits -= num.Log2(ct0.noise / ne.noise.RoundNoise())
	}

	auxModLen := 1
	for auxModLen <= len(ne.ambMod)-tarLen {
		tarBits -= num.Log2(ne.ambMod[tarLen+auxModLen-1].Value())
		if tarBits <= 0 {
			break
		}
		auxModLen += 1
	}

	return auxModLen
}

// liftAndTensorTo computes the noise of the tensor product of two ciphertexts.
func (ne *NoiseEstimator) liftAndTensor(ct0, ct1 *Ciphertext, ambLen int) float64 {
	// Force ct0 to have the smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	tarLen := ct0.ModLen()
	noise0 := ct0.noise
	noise1 := ct1.noise

	// Compute the auxiliary and base moduli.
	logBaseMod := float64(0)
	for i := 0; i < tarLen; i++ {
		logBaseMod += num.Log2(float64(ne.ambMod[i].Value()))
	}
	logAmbMod := float64(0)
	for i := tarLen; i < tarLen+ambLen; i++ {
		logAmbMod += num.Log2(float64(ne.ambMod[i].Value()))
	}

	// Modulus-switch ct0.
	rdNoise := ne.noise.RoundNoise()
	scale := math.Exp2(logAmbMod - logBaseMod)
	switch ne.estimType {
	case heint.VarianceType:
		noise1 = noise1 * scale * scale
	case heint.WorstCaseType:
		noise1 = noise1 * scale
	}
	noise1 += rdNoise

	// After tensoring we have PQ/t*mm' + Q(m'+tI') * e + P(m+tI) * e' + tee'
	// For floating point arithmetic, we scale this by P.
	// Essentially the noise term will be Q/P(m'+tI') * e + (m+tI) * e' + t/P * ee'.
	msgMod := float64(ne.msgMod.Value())
	expFac := float64(ne.params.RingParams().ExpandFactor())
	scale = math.Exp2(logBaseMod - logAmbMod)
	res := float64(0)
	switch ne.estimType {
	case heint.VarianceType:
		res += scale * scale * (msgMod + msgMod*msgMod*rdNoise) * noise1 * expFac
		res += (msgMod + msgMod*msgMod*rdNoise) * noise0 * expFac
		res += float64(1) / float64(12)
	case heint.WorstCaseType:
		res += scale * (msgMod + msgMod*rdNoise) * noise1 * expFac
		res += (msgMod + msgMod*rdNoise) * noise0 * expFac
		res += 0.5
	}

	return res
}

// TODO: Output is incorrect. Revise.
// Mul returns the noise of the product of two ciphertexts.
func (ne *NoiseEstimator) MulTo(ctOut, ct0, ct1 *Ciphertext) {
	// Force ct0 to have the smaller modulus.
	if ct0.ModLen() > ct1.ModLen() || (ct0.ModLen() == ct1.ModLen() && ct0.noise < ct1.noise) {
		ct0, ct1 = ct1, ct0
	}

	tarLen := ct0.ModLen()
	ambLen := ne.getAmbLen(ct0, ct1)
	res := ne.liftAndTensor(ct0, ct1, ambLen)
	res += ne.noise.RoundNoise()
	res += ne.noise.GadgetProd(tarLen)

	ctOut.noise = res
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

// MulPlainMatrixTo returns the noise of the product of a ciphertext and a plaintext matrix.
func (ne *NoiseEstimator) MulPlainMatrixTo(cOut, ct *Ciphertext, mat *rlwe.PlainMatrix) {
	cOut.noise = ne.noise.MulPlainMatrix(ct.noise, mat, ct.ModLen())
}
