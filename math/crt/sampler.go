package crt

import (
	"crypto/rand"
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

const (
	GaussianTailCut = 12
)

// SamplerParameters is the parameters for [Sampler].
type SamplerParameters interface {
	// Sampler returns the [Sampler].
	Sampler() Sampler
	// Bound returns the inf-norm bound of the sampler.
	Bound() float64
	// Variance returns the inf-norm bound of the variance of the sampler.
	Variance() float64
}

// UniformSamplerParameters is the parameters for [UniformSampler].
type UniformSamplerParameters struct {
	// BoundMin is the minimum bound that is sampled.
	// If BoundMax is set but BoundMin is nil, then it is set to -q/2.
	// If both BoundMin and BoundMax are nil, then it is ignored.
	BoundMin *big.Int
	// BoundMax is the maximum bound (exclusive) that is sampled.
	// If BoundMin is set but BoundMax is nil, then it is set to q/2.
	// If both BoundMin and BoundMax are nil, then it is ignored.
	BoundMax *big.Int
}

// Sampler returns the [Sampler].
func (p UniformSamplerParameters) Sampler() Sampler {
	return &UniformSampler{
		baseSampler: csprng.NewUniformSampler(),

		boundMin: p.BoundMin,
		boundMax: p.BoundMax,
	}
}

// Bound returns the inf-norm bound of the sampler.
func (p UniformSamplerParameters) Bound() float64 {
	if p.BoundMin == nil || p.BoundMax == nil {
		panic("bound min and max must be set")
	}

	bMax, _ := p.BoundMax.Float64()
	bMin, _ := p.BoundMin.Float64()
	return max(math.Abs(bMax-1), math.Abs(bMin))
}

// Variance returns the inf-norm bound of the variance of the sampler.
func (p UniformSamplerParameters) Variance() float64 {
	if p.BoundMin == nil || p.BoundMax == nil {
		panic("bound min and max must be set")
	}

	bMax, _ := p.BoundMax.Float64()
	bMin, _ := p.BoundMin.Float64()

	return ((bMax-bMin)*(bMax-bMin) - 1) / 12
}

// TernarySamplerParameters is the parameters for [TernarySampler].
type TernarySamplerParameters struct {
	// Positive is the probability of sampling a 1.
	// Panics if not in [0, 1] or Positive+Negative > 1.
	Positive float64
	// Negative is the probability of sampling a -1.
	// Panics if not in [0, 1] or Positive+Negative > 1.
	Negative float64
	// HammingWeight is the hamming weight of the polynomial.
	// Panics if not in [0, Rank].
	// If set to 0, then it is ignored.
	HammingWeight int
}

// Sampler returns the [Sampler].
func (p TernarySamplerParameters) Sampler() Sampler {
	if p.Positive < 0 || p.Positive > 1 {
		panic("positive must be in [0, 1]")
	} else if p.Negative < 0 || p.Negative > 1 {
		panic("negative must be in [0, 1]")
	} else if p.Positive+p.Negative > 1 {
		panic("positive + negative must be in [0, 1]")
	} else if p.HammingWeight < 0 {
		panic("hamming weight must be in [0, rank]")
	} else if p.HammingWeight > 0 && p.Positive+p.Negative == 0 {
		panic("positive+negative must be nonzero if hamming weight is set")
	}

	return &TernarySampler{
		baseSampler: csprng.NewUniformSampler(),

		posFloat: p.Positive,
		negFloat: p.Negative,

		pos: uint64(p.Positive * math.Exp2(63)),
		neg: uint64(p.Negative * math.Exp2(63)),
		hw:  p.HammingWeight,
	}
}

// Bound returns the inf-norm bound of the sampler.
func (p TernarySamplerParameters) Bound() float64 {
	return 1
}

// Variance returns the inf-norm bound of the variance of the sampler.
func (p TernarySamplerParameters) Variance() float64 {
	if p.HammingWeight > 0 {
		panic("hamming weight must be 0")
	}

	return p.Positive + p.Negative
}

// RoundedGaussianSamplerParameters is the parameters for [RoundedGaussianSampler].
type RoundedGaussianSamplerParameters struct {
	// StdDev is the standard deviation of the Gaussian distribution.
	StdDev float64
}

// Sampler returns the [Sampler].
func (p RoundedGaussianSamplerParameters) Sampler() Sampler {
	return &RoundedGaussianSampler{
		baseSampler: csprng.NewRoundedGaussianSampler(),
		stdDev:      p.StdDev,
	}
}

// Bound returns the inf-norm bound of the sampler.
func (p RoundedGaussianSamplerParameters) Bound() float64 {
	return p.StdDev * GaussianTailCut
}

// Variance returns the inf-norm bound of the variance of the sampler.
func (p RoundedGaussianSamplerParameters) Variance() float64 {
	return p.StdDev*p.StdDev + 0.5
}

// Sampler is an interface for sampling *[Poly].
type Sampler interface {
	// SamplerParams returns the sampler parameters.
	SamplerParams() SamplerParameters
	// Sample samples a *[Poly] with given rank and modulus in standard domain.
	Sample(rank int, mod []*num.Modulus) *Poly
	// SampleTo samples a *[Poly] to pOut in standard domain.
	SampleTo(pOut *Poly, mod []*num.Modulus)
}

// UniformSampler is a sampler for uniform distribution.
type UniformSampler struct {
	baseSampler *csprng.UniformSampler
	boundMin    *big.Int
	boundMax    *big.Int
}

// SamplerParams returns the sampler parameters.
func (s *UniformSampler) SamplerParams() SamplerParameters {
	return UniformSamplerParameters{
		BoundMin: s.boundMin,
		BoundMax: s.boundMax,
	}
}

// Sample samples a *[Poly] with given rank and modulus in standard domain.
func (s *UniformSampler) Sample(rank int, mod []*num.Modulus) *Poly {
	pOut := NewPoly(rank, len(mod))
	s.SampleTo(pOut, mod)
	return pOut
}

// SampleTo samples a *[Poly] to pOut in standard domain.
func (s *UniformSampler) SampleTo(pOut *Poly, mod []*num.Modulus) {
	checkPolyShape(pOut.Rank(), len(mod), pOut)

	if s.boundMin != nil || s.boundMax != nil {
		s.sampleToBounded(pOut, mod)
		pOut.Form = dft.FormCoeff
		return
	}

	for l := range mod {
		for i := range pOut.Coeffs[l] {
			pOut.Coeffs[l][i] = s.baseSampler.SampleN(0, mod[l].Value())
		}
	}

	pOut.Form = dft.FormCoeff
}

func (s *UniformSampler) sampleToBounded(pOut *Poly, mod []*num.Modulus) {
	modBig := make([]*big.Int, len(mod))
	modProd := big.NewInt(1)
	for i := range mod {
		modBig[i] = new(big.Int).SetUint64(mod[i].Value())
		modProd.Mul(modProd, modBig[i])
	}

	boundMin := s.boundMin
	if boundMin == nil {
		boundMin = new(big.Int).Rsh(modProd, 1)
		boundMin.Neg(boundMin)
	}

	boundMax := s.boundMax
	if boundMax == nil {
		boundMax = new(big.Int).Rsh(modProd, 1)
	}

	boundDiff := new(big.Int).Sub(boundMax, boundMin)
	for j := 0; j < pOut.Rank(); j++ {
		cInt, _ := rand.Int(s.baseSampler, boundDiff)
		cInt.Add(cInt, boundMin)

		for i := range mod {
			pOut.Coeffs[i][j] = new(big.Int).Mod(cInt, modBig[i]).Uint64()
		}
	}
}

// TernarySampler is a sampler for ternary distribution.
type TernarySampler struct {
	baseSampler *csprng.UniformSampler

	posFloat float64
	negFloat float64

	pos uint64
	neg uint64
	hw  int
}

// SamplerParams returns the sampler parameters.
func (s *TernarySampler) SamplerParams() SamplerParameters {
	return TernarySamplerParameters{
		Positive:      s.posFloat,
		Negative:      s.negFloat,
		HammingWeight: s.hw,
	}
}

// Sample samples a *[Poly] with given rank and modulus in standard domain.
func (s *TernarySampler) Sample(rank int, mod []*num.Modulus) *Poly {
	pOut := NewPoly(rank, len(mod))
	s.SampleTo(pOut, mod)
	return pOut
}

// SampleTo samples a *[Poly] to pOut in standard domain.
func (s *TernarySampler) SampleTo(pOut *Poly, mod []*num.Modulus) {
	checkPolyShape(pOut.Rank(), len(mod), pOut)
	if s.hw > pOut.Rank() {
		panic("hamming weight must be less than or equal to rank")
	}

	var c int64

	if s.hw == 0 {
		for j := 0; j < pOut.Rank(); j++ {
			r := s.baseSampler.Sample[uint64]() >> 1
			switch {
			case r < s.pos:
				c = 1
			case r < s.pos+s.neg:
				c = -1
			default:
				c = 0
			}

			switch c {
			case 1:
				for i := range mod {
					pOut.Coeffs[i][j] = 1
				}
			case -1:
				for i := range mod {
					pOut.Coeffs[i][j] = mod[i].Value() - 1
				}
			case 0:
				for i := range mod {
					pOut.Coeffs[i][j] = 0
				}
			}
		}

		pOut.Form = dft.FormCoeff

		return
	}

	for j := 0; j < s.hw; j++ {
		r := s.baseSampler.SampleN(0, s.pos+s.neg)
		switch {
		case r < s.pos:
			c = 1
		default:
			c = -1
		}

		switch c {
		case 1:
			for i := range mod {
				pOut.Coeffs[i][j] = 1
			}
		case -1:
			for i := range mod {
				pOut.Coeffs[i][j] = mod[i].Value() - 1
			}
		}
	}

	for i := range mod {
		clear(pOut.Coeffs[i][s.hw:])
	}

	for j := 1; j < pOut.Rank(); j++ {
		k := s.baseSampler.SampleN(0, uint64(j+1))
		for i := range mod {
			pOut.Coeffs[i][j], pOut.Coeffs[i][k] = pOut.Coeffs[i][k], pOut.Coeffs[i][j]
		}
	}

	pOut.Form = dft.FormCoeff
}

// RoundedGaussianSampler is a sampler for rounded Gaussian distribution.
type RoundedGaussianSampler struct {
	baseSampler *csprng.RoundedGaussianSampler
	stdDev      float64
}

// SamplerParams returns the sampler parameters.
func (s *RoundedGaussianSampler) SamplerParams() SamplerParameters {
	return RoundedGaussianSamplerParameters{
		StdDev: s.stdDev,
	}
}

// Sample samples a *[Poly] with given rank and modulus in standard domain.
func (s *RoundedGaussianSampler) Sample(rank int, mod []*num.Modulus) *Poly {
	pOut := NewPoly(rank, len(mod))
	s.SampleTo(pOut, mod)
	return pOut
}

// SampleTo samples a *[Poly] to pOut in standard domain.
func (s *RoundedGaussianSampler) SampleTo(pOut *Poly, mod []*num.Modulus) {
	checkPolyShape(pOut.Rank(), len(mod), pOut)

	if math.Log2(s.stdDev) <= num.MaxModulusBits {
		s.sampleFloat64To(pOut, mod, s.stdDev)
	} else {
		s.sampleBigFloatTo(pOut, mod, s.stdDev)
	}
}

func (s *RoundedGaussianSampler) sampleFloat64To(pOut *Poly, mod []*num.Modulus, stdDev float64) {
	for j := 0; j < pOut.Rank(); j++ {
		c := s.baseSampler.Sample(stdDev)
		for i := range mod {
			pOut.Coeffs[i][j] = num.Reduce(c, mod[i])
		}
	}

	pOut.Form = dft.FormCoeff
}

func (s *RoundedGaussianSampler) sampleBigFloatTo(pOut *Poly, mod []*num.Modulus, stdDev float64) {
	modBig := make([]*big.Int, len(mod))
	for i := range mod {
		modBig[i] = new(big.Int).SetUint64(mod[i].Value())
	}

	for j := 0; j < pOut.Rank(); j++ {
		cInt := s.baseSampler.SampleBig(stdDev)
		for i := range mod {
			pOut.Coeffs[i][j] = new(big.Int).Mod(cInt, modBig[i]).Uint64()
		}
	}

	pOut.Form = dft.FormCoeff
}
