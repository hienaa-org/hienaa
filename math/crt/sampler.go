package crt

import (
	"crypto/rand"
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/num"
)

// SamplerParameters is the parameters for [Sampler].
type SamplerParameters interface {
	// Sampler returns the [Sampler].
	Sampler() Sampler
}

// UniformSamplerParameters is the parameters for [UniformSampler].
type UniformSamplerParameters struct {
	// BoundMin is the minimum bound that is sampled.
	// If BoundMax is set but BoundMin is nil, then it is set to -q/2.
	// If both BoundMin and BoundMax are nil, then it is ignored.
	BoundMin *big.Int
	// BoundMax is the maximum bound that is sampled.
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

// RoundedGaussianSamplerParameters is the parameters for [RoundedGaussianSampler].
type RoundedGaussianSamplerParameters[T float64 | *big.Float] struct {
	// Center is the center of the Gaussian distribution.
	// If nil, it is set to 0.
	Center T
	// StdDev is the standard deviation of the Gaussian distribution.
	StdDev T
}

// Sampler returns the [Sampler].
func (p RoundedGaussianSamplerParameters[T]) Sampler() Sampler {
	switch any(p.StdDev).(type) {
	case float64:
		stdDev := any(p.StdDev).(float64)
		if stdDev <= 0 || stdDev >= math.Exp2(64) {
			panic("standard deviation must be in (0, 2^64)")
		}

		return &RoundedGaussianSampler[float64]{
			baseSampler: csprng.NewRoundedGaussianSampler(),

			center: any(p.Center).(float64),
			stdDev: stdDev,
		}

	case *big.Float:
		stdDev := any(p.StdDev).(*big.Float)
		if stdDev.Cmp(big.NewFloat(0)) != 1 {
			panic("standard deviation must be positive")
		}

		center := any(p.Center).(*big.Float)
		if center == nil {
			center = big.NewFloat(0)
		}

		return &RoundedGaussianSampler[*big.Float]{
			baseSampler: csprng.NewRoundedGaussianSampler(),

			center: center,
			stdDev: stdDev,
		}
	}

	panic("unsupported parameters")
}

// Sampler is an interface for sampling polynomials.
type Sampler interface {
	// SamplerParams returns the sampler parameters.
	SamplerParams() SamplerParameters
	// Sample samples an [Element] with given rank and modulus in standard domain.
	Sample(rank int, mod []*num.Modulus) *Element
	// SampleTo samples an [Element] to eOut in standard domain.
	SampleTo(eOut *Element, mod []*num.Modulus)
	// SafeCopy returns a thread-safe copy.
	// This always returns a freshly seeded sampler.
	SafeCopy() Sampler
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

// Sample samples an [Element] with given rank and modulus in standard domain.
func (s *UniformSampler) Sample(rank int, mod []*num.Modulus) *Element {
	eOut := NewPoly(rank, len(mod))
	s.SampleTo(eOut, mod)
	return eOut
}

// SampleTo samples an [Element] to eOut in standard domain.
func (s *UniformSampler) SampleTo(eOut *Element, mod []*num.Modulus) {
	checkShape(eOut.Rank(), len(mod), eOut)

	if s.boundMin != nil || s.boundMax != nil {
		s.sampleToBounded(eOut, mod)
		eOut.IsNTT = false
		return
	}

	for i := range mod {
		for j := range eOut.Coeffs[i] {
			eOut.Coeffs[i][j] = s.baseSampler.SampleN(mod[i].Value())
		}
	}

	eOut.IsNTT = false
}

func (s *UniformSampler) sampleToBounded(eOut *Element, mod []*num.Modulus) {
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
	for j := 0; j < eOut.Rank(); j++ {
		cInt, _ := rand.Int(s.baseSampler, boundDiff)
		cInt.Add(cInt, boundMin)

		for i := range mod {
			eOut.Coeffs[i][j] = new(big.Int).Mod(cInt, modBig[i]).Uint64()
		}
	}
}

// SafeCopy returns a thread-safe copy.
// This always returns a freshly seeded sampler.
func (s *UniformSampler) SafeCopy() Sampler {
	return &UniformSampler{
		baseSampler: csprng.NewUniformSampler(),
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

// Sample samples an [Element] with given rank and modulus in standard domain.
func (s *TernarySampler) Sample(rank int, mod []*num.Modulus) *Element {
	eOut := NewPoly(rank, len(mod))
	s.SampleTo(eOut, mod)
	return eOut
}

// SampleTo samples an [Element] to eOut in standard domain.
func (s *TernarySampler) SampleTo(eOut *Element, mod []*num.Modulus) {
	checkShape(eOut.Rank(), len(mod), eOut)
	if s.hw > eOut.Rank() {
		panic("hamming weight must be less than or equal to rank")
	}

	var c int64

	if s.hw == 0 {
		for j := 0; j < eOut.Rank(); j++ {
			r := s.baseSampler.Sample() >> 1
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
					eOut.Coeffs[i][j] = 1
				}
			case -1:
				for i := range mod {
					eOut.Coeffs[i][j] = mod[i].Value() - 1
				}
			case 0:
				for i := range mod {
					eOut.Coeffs[i][j] = 0
				}
			}
		}

		return
	}

	for j := 0; j < s.hw; j++ {
		r := s.baseSampler.SampleN(s.pos + s.neg)
		switch {
		case r < s.pos:
			c = 1
		default:
			c = -1
		}

		switch c {
		case 1:
			for i := range mod {
				eOut.Coeffs[i][j] = 1
			}
		case -1:
			for i := range mod {
				eOut.Coeffs[i][j] = mod[i].Value() - 1
			}
		}
	}

	for i := range mod {
		clear(eOut.Coeffs[i][s.hw:])
	}

	for j := 1; j < eOut.Rank(); j++ {
		k := s.baseSampler.SampleN(uint64(j + 1))
		for i := range mod {
			eOut.Coeffs[i][j], eOut.Coeffs[i][k] = eOut.Coeffs[i][k], eOut.Coeffs[i][j]
		}
	}

	eOut.IsNTT = false
}

// SafeCopy returns a thread-safe copy.
// This always returns a freshly seeded sampler.
func (s *TernarySampler) SafeCopy() Sampler {
	return &TernarySampler{
		baseSampler: csprng.NewUniformSampler(),

		posFloat: s.posFloat,
		negFloat: s.negFloat,

		pos: s.pos,
		neg: s.neg,
		hw:  s.hw,
	}
}

// RoundedGaussianSampler is a sampler for rounded Gaussian distribution.
type RoundedGaussianSampler[T float64 | *big.Float] struct {
	baseSampler *csprng.RoundedGaussianSampler

	center T
	stdDev T
}

// SamplerParams returns the sampler parameters.
func (s *RoundedGaussianSampler[T]) SamplerParams() SamplerParameters {
	return RoundedGaussianSamplerParameters[T]{
		Center: s.center,
		StdDev: s.stdDev,
	}
}

// Sample samples an [Element] with given rank and modulus in standard domain.
func (s *RoundedGaussianSampler[T]) Sample(rank int, mod []*num.Modulus) *Element {
	eOut := NewPoly(rank, len(mod))
	s.SampleTo(eOut, mod)
	return eOut
}

// SampleTo samples an [Element] to eOut in standard domain.
func (s *RoundedGaussianSampler[T]) SampleTo(eOut *Element, mod []*num.Modulus) {
	checkShape(eOut.Rank(), len(mod), eOut)

	switch any(s.stdDev).(type) {
	case float64:
		s.sampleFloat64To(eOut, mod, any(s.center).(float64), any(s.stdDev).(float64))
	case *big.Float:
		s.sampleBigFloatTo(eOut, mod, any(s.center).(*big.Float), any(s.stdDev).(*big.Float))
	}
}

func (s *RoundedGaussianSampler[T]) sampleFloat64To(eOut *Element, mod []*num.Modulus, center, stdDev float64) {
	for j := 0; j < eOut.Rank(); j++ {
		c := s.baseSampler.Sample(center, stdDev)
		for i := range mod {
			eOut.Coeffs[i][j] = num.Reduce(c, mod[i])
		}
	}

	eOut.IsNTT = false
}

func (s *RoundedGaussianSampler[T]) sampleBigFloatTo(eOut *Element, mod []*num.Modulus, center, stdDev *big.Float) {
	modBig := make([]*big.Int, len(mod))
	for i := range mod {
		modBig[i] = new(big.Int).SetUint64(mod[i].Value())
	}

	for j := 0; j < eOut.Rank(); j++ {
		cInt := s.baseSampler.SampleBig(center, stdDev)
		for i := range mod {
			eOut.Coeffs[i][j] = new(big.Int).Mod(cInt, modBig[i]).Uint64()
		}
	}

	eOut.IsNTT = false
}

// SafeCopy returns a thread-safe copy.
// This always returns a freshly seeded sampler.
func (s *RoundedGaussianSampler[T]) SafeCopy() Sampler {
	return &RoundedGaussianSampler[T]{
		baseSampler: csprng.NewRoundedGaussianSampler(),

		center: s.center,
		stdDev: s.stdDev,
	}
}
