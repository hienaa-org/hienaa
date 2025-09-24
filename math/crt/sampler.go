package crt

import (
	"crypto/rand"
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// SamplerParameters is the parameters for [Sampler].
type SamplerParameters interface {
	// Sampler returns the [Sampler].
	Sampler() Sampler
}

// UniformSamplerParameters is the parameters for [UniformSampler].
type UniformSamplerParameters struct {
	Params  dft.RingParameters
	Modulus []*num.Modulus

	// BoundMin is the minimum bound that is sampled.
	// If BoundMax is set but BoundMin is nil, then it is set to -q/2.
	// If both BoundMin and BoundMax are nil, then it is ignored.
	BoundMin *big.Int
	// BoundMax is the maximum bound that is sampled.
	// If BoundMin is set but BoundMax is nil, then it is set to q/2.
	// If both BoundMin and BoundMax are nil, then it is ignored.
	BoundMax *big.Int
}

func (p UniformSamplerParameters) Sampler() Sampler {
	var modBig []*big.Int
	var modProd *big.Int
	if p.BoundMin != nil || p.BoundMax != nil {
		modBig = make([]*big.Int, len(p.Modulus))
		modProd = big.NewInt(1)
		for i := range p.Modulus {
			modBig[i] = new(big.Int).SetUint64(p.Modulus[i].Value())
			modProd.Mul(modProd, modBig[i])
		}
	}

	return &UniformSampler{
		params:  p.Params,
		mod:     p.Modulus,
		modBig:  modBig,
		modProd: modProd,

		baseSampler: csprng.NewUniformSampler(),

		boundMin: p.BoundMin,
		boundMax: p.BoundMax,
	}
}

// TernarySamplerParameters is the parameters for [TernarySampler].
type TernarySamplerParameters struct {
	Params  dft.RingParameters
	Modulus []*num.Modulus

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

func (p TernarySamplerParameters) Sampler() Sampler {
	switch {
	case p.Positive < 0 || p.Positive > 1:
		panic("TernarySamplerParameters: Positive must be in [0, 1]")
	case p.Negative < 0 || p.Negative > 1:
		panic("TernarySamplerParameters: Negative must be in [0, 1]")
	case p.Positive+p.Negative > 1:
		panic("TernarySamplerParameters: Positive + Negative must be in [0, 1]")
	case p.HammingWeight < 0 || p.HammingWeight > p.Params.Rank():
		panic("TernarySamplerParameters: HammingWeight must be in [0, Rank]")
	}

	return &TernarySampler{
		params: p.Params,
		mod:    p.Modulus,

		baseSampler: csprng.NewUniformSampler(),

		posFloat: p.Positive,
		negFloat: p.Negative,

		pos: uint64(p.Positive * math.Exp2(63)),
		neg: uint64(p.Negative * math.Exp2(63)),
		hw:  p.HammingWeight,
	}
}

// RoundedGaussianSamplerParameters is the parameters for [RoundedGaussianSampler].
// RoundedGaussianSampler only supports small standard deviations.
// For large devations, use [BigRoundedGaussianSamplerParameters].
type RoundedGaussianSamplerParameters[T float64 | *big.Float] struct {
	Params  dft.RingParameters
	Modulus []*num.Modulus

	// Center is the center of the Gaussian distribution.
	// If nil, it is set to 0.
	Center T
	// StdDev is the standard deviation of the Gaussian distribution.
	StdDev T
}

func (p RoundedGaussianSamplerParameters[T]) Sampler() Sampler {
	switch any(p.StdDev).(type) {
	case float64:
		stdDev := any(p.StdDev).(float64)
		if stdDev <= 0 || stdDev >= math.Exp2(64) {
			panic("RoundedGaussianSamplerParameters: StdDev must be in (0, 2^64)")
		}

		return &RoundedGaussianSampler[float64]{
			params: p.Params,
			mod:    p.Modulus,

			baseSampler: csprng.NewRoundedGaussianSampler(),

			center: any(p.Center).(float64),
			stdDev: stdDev,
		}

	case *big.Float:
		stdDev := any(p.StdDev).(*big.Float)
		if stdDev.Cmp(big.NewFloat(0)) != 1 {
			panic("RoundedGaussianSamplerParameters: StdDev must be positive")
		}

		center := any(p.Center).(*big.Float)
		if center == nil {
			center = big.NewFloat(0)
		}

		modBig := make([]*big.Int, len(p.Modulus))
		for i := range p.Modulus {
			modBig[i] = new(big.Int).SetUint64(p.Modulus[i].Value())
		}

		return &RoundedGaussianSampler[*big.Float]{
			params: p.Params,
			mod:    p.Modulus,
			modBig: modBig,

			baseSampler: csprng.NewRoundedGaussianSampler(),

			center: center,
			stdDev: stdDev,
		}
	}

	panic("RoundedGaussianSamplerParameters: unsupported parameters")
}

// Sampler is an interface for sampling polynomials.
type Sampler interface {
	// Params returns the ring parameters.
	Params() dft.RingParameters
	// Modulus returns the modulus.
	Modulus() []*num.Modulus
	// SamplerParams returns the sampler parameters.
	SamplerParams() SamplerParameters
	// Sample samples a polynomial.
	// Output is always Standard form.
	Sample() *Poly
	// SampleTo samples a polynomial to pOut.
	// pOut is set to Standard form.
	SampleTo(pOut *Poly)
	// SafeCopy returns a thread-safe copy.
	// This always returns a freshly seeded sampler.
	SafeCopy() Sampler
}

// UniformSampler is a sampler for uniform distribution.
type UniformSampler struct {
	params  dft.RingParameters
	mod     []*num.Modulus
	modBig  []*big.Int
	modProd *big.Int

	baseSampler *csprng.UniformSampler

	boundMin *big.Int
	boundMax *big.Int
}

func (s *UniformSampler) Params() dft.RingParameters {
	return s.params
}

func (s *UniformSampler) Modulus() []*num.Modulus {
	return s.mod
}

func (s *UniformSampler) SamplerParams() SamplerParameters {
	return UniformSamplerParameters{
		Params:  s.params,
		Modulus: s.mod,

		BoundMin: s.boundMin,
		BoundMax: s.boundMax,
	}
}

func (s *UniformSampler) Sample() *Poly {
	p := NewPoly(s.params.Rank(), len(s.mod))
	s.SampleTo(p)
	return p
}

func (s *UniformSampler) SampleTo(pOut *Poly) {
	if !isConsistent(s.params.Rank(), len(s.mod), pOut) {
		panic("SampleTo: inputs not consistent")
	}

	if s.boundMin != nil || s.boundMax != nil {
		s.sampleToBounded(pOut)
		pOut.isNTT = false
		return
	}

	for i := range s.mod {
		for j := range pOut.Coeffs[i] {
			pOut.Coeffs[i][j] = s.baseSampler.SampleN(s.mod[i].Value())
		}
	}

	pOut.isNTT = false
}

func (s *UniformSampler) sampleToBounded(pOut *Poly) {
	boundMin := s.boundMin
	if boundMin == nil {
		boundMin = new(big.Int).Rsh(s.modProd, 1)
		boundMin.Neg(boundMin)
	}

	boundMax := s.boundMax
	if boundMax == nil {
		boundMax = new(big.Int).Rsh(s.modProd, 1)
	}

	boundDiff := new(big.Int).Sub(boundMax, boundMin)
	for j := 0; j < s.params.Rank(); j++ {
		cInt, _ := rand.Int(s.baseSampler, boundDiff)
		cInt.Add(cInt, boundMin)

		for i := range s.mod {
			pOut.Coeffs[i][j] = new(big.Int).Mod(cInt, s.modBig[i]).Uint64()
		}
	}
}

func (s *UniformSampler) SafeCopy() Sampler {
	return &UniformSampler{
		params:      s.params,
		mod:         s.mod,
		baseSampler: csprng.NewUniformSampler(),
	}
}

type TernarySampler struct {
	params dft.RingParameters
	mod    []*num.Modulus

	baseSampler *csprng.UniformSampler

	posFloat float64
	negFloat float64

	pos uint64
	neg uint64
	hw  int
}

func (s *TernarySampler) Params() dft.RingParameters {
	return s.params
}

func (s *TernarySampler) Modulus() []*num.Modulus {
	return s.mod
}

func (s *TernarySampler) SamplerParams() SamplerParameters {
	return TernarySamplerParameters{
		Params:  s.params,
		Modulus: s.mod,

		Positive:      s.posFloat,
		Negative:      s.negFloat,
		HammingWeight: s.hw,
	}
}

func (s *TernarySampler) Sample() *Poly {
	p := NewPoly(s.params.Rank(), len(s.mod))
	s.SampleTo(p)
	return p
}

func (s *TernarySampler) SampleTo(pOut *Poly) {
	if !isConsistent(s.params.Rank(), len(s.mod), pOut) {
		panic("SampleTo: inputs not consistent")
	}

	var c int64

	if s.hw == 0 {
		for j := 0; j < s.params.Rank(); j++ {
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
				for i := range s.mod {
					pOut.Coeffs[i][j] = 1
				}
			case -1:
				for i := range s.mod {
					pOut.Coeffs[i][j] = s.mod[i].Value() - 1
				}
			case 0:
				for i := range s.mod {
					pOut.Coeffs[i][j] = 0
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
			for i := range s.mod {
				pOut.Coeffs[i][j] = 1
			}
		case -1:
			for i := range s.mod {
				pOut.Coeffs[i][j] = s.mod[i].Value() - 1
			}
		}
	}

	for i := range s.mod {
		clear(pOut.Coeffs[i][s.hw:])
	}

	for j := 1; j < s.params.Rank(); j++ {
		k := s.baseSampler.SampleN(uint64(j + 1))
		for i := range s.mod {
			pOut.Coeffs[i][j], pOut.Coeffs[i][k] = pOut.Coeffs[i][k], pOut.Coeffs[i][j]
		}
	}

	pOut.isNTT = false
}

func (s *TernarySampler) SafeCopy() Sampler {
	return &TernarySampler{
		params: s.params,
		mod:    s.mod,

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
	params dft.RingParameters
	mod    []*num.Modulus
	modBig []*big.Int

	baseSampler *csprng.RoundedGaussianSampler

	center T
	stdDev T
}

func (s *RoundedGaussianSampler[T]) Params() dft.RingParameters {
	return s.params
}

func (s *RoundedGaussianSampler[T]) Modulus() []*num.Modulus {
	return s.mod
}

func (s *RoundedGaussianSampler[T]) SamplerParams() SamplerParameters {
	return RoundedGaussianSamplerParameters[T]{
		Params:  s.params,
		Modulus: s.mod,

		Center: s.center,
		StdDev: s.stdDev,
	}
}

func (s *RoundedGaussianSampler[T]) Sample() *Poly {
	p := NewPoly(s.params.Rank(), len(s.mod))
	s.SampleTo(p)
	return p
}

func (s *RoundedGaussianSampler[T]) SampleTo(pOut *Poly) {
	if !isConsistent(s.params.Rank(), len(s.mod), pOut) {
		panic("SampleTo: inputs not consistent")
	}

	switch any(s.stdDev).(type) {
	case float64:
		s.sampleFloat64To(pOut, any(s.center).(float64), any(s.stdDev).(float64))
	case *big.Float:
		s.sampleBigFloatTo(pOut, any(s.center).(*big.Float), any(s.stdDev).(*big.Float))
	}
}

func (s *RoundedGaussianSampler[T]) sampleFloat64To(pOut *Poly, center, stdDev float64) {
	for j := 0; j < s.params.Rank(); j++ {
		c := s.baseSampler.Sample(center, stdDev)
		for i := range s.mod {
			pOut.Coeffs[i][j] = num.Reduce(c, s.mod[i])
		}
	}

	pOut.isNTT = false
}

func (s *RoundedGaussianSampler[T]) sampleBigFloatTo(pOut *Poly, center, stdDev *big.Float) {
	for j := 0; j < s.params.Rank(); j++ {
		cInt := s.baseSampler.SampleBig(center, stdDev)
		for i := range s.mod {
			pOut.Coeffs[i][j] = new(big.Int).Mod(cInt, s.modBig[i]).Uint64()
		}
	}

	pOut.isNTT = false
}

func (s *RoundedGaussianSampler[T]) SafeCopy() Sampler {
	return &RoundedGaussianSampler[T]{
		params: s.params,
		mod:    s.mod,

		stdDev:      s.stdDev,
		baseSampler: csprng.NewRoundedGaussianSampler(),
	}
}
