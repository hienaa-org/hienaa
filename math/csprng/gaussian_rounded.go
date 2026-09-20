package csprng

import (
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

const (
	blockSize = 128
	floatPrec = 52
	roundPrec = 49

	// From "The Ziggurat Method for Generating Random Variables"
	// by Marsaglia and Tsang (2000)
	rn = 3.442619855899
)

var (
	kn [blockSize]uint64
	wn [blockSize]float64
	fn [blockSize]float64
)

func init() {
	v := rn*normal(rn) + normalIntegral(rn)

	var xn [blockSize]float64
	xn[blockSize-1] = rn
	for i := blockSize - 2; i >= 1; i-- {
		xn[i] = normalInv(v/xn[i+1] + normal(xn[i+1]))
	}

	const scale = 1 << floatPrec
	for i := 1; i < blockSize; i++ {
		kn[i] = uint64((xn[i-1] / xn[i]) * scale)
		wn[i] = xn[i] / scale
		fn[i] = normal(xn[i])
	}
	kn[0] = uint64((rn * normal(rn) / v) * scale)
	wn[0] = (v / normal(rn)) / scale
	fn[0] = 1
}

func normal(x float64) float64 {
	return math.Exp(-0.5 * x * x)
}

func normalIntegral(x float64) float64 {
	return math.Sqrt(math.Pi/2) * math.Erfc(x/math.Sqrt(2))
}

func normalInv(x float64) float64 {
	return math.Sqrt(-2 * math.Log(x))
}

// RoundedGaussianSampler samples from Rounded Gaussian Distribution.
// The underlying normal draw uses float64, with at most 53 bits of precision.
type RoundedGaussianSampler struct {
	baseSampler *UniformSampler
}

// NewRoundedGaussianSampler creates a new [RoundedGaussianSampler].
//
// Panics when read from crypto/rand or AES initialization fails.
func NewRoundedGaussianSampler() *RoundedGaussianSampler {
	return &RoundedGaussianSampler{
		baseSampler: NewUniformSampler(),
	}
}

// NewRoundedGaussianSamplerWithSeed creates a new [RoundedGaussianSampler], with user supplied seed.
//
// Panics when AES initialization fails.
func NewRoundedGaussianSamplerWithSeed(seed []byte) *RoundedGaussianSampler {
	return &RoundedGaussianSampler{
		baseSampler: NewUniformSamplerWithSeed(seed),
	}
}

// normFloat samples float64 from normal distribution.
func (s *RoundedGaussianSampler) normFloat() float64 {
	for {
		r := s.baseSampler.Sample[uint64]()

		b := r >> 63
		i := r % (1 << 7)
		j := (r >> 7) % (1 << floatPrec)

		x := float64(int64((j^-b)+b)) * wn[i]

		if j < kn[i] {
			return x
		}

		if i == 0 {
			var u, v float64
			for {
				u = -math.Log(s.baseSampler.SampleFloat()) * (1.0 / rn)
				v = -math.Log(s.baseSampler.SampleFloat())
				if v+v >= u*u {
					break
				}
			}
			u += rn
			if math.IsInf(u, 0) {
				continue
			}

			if b == 1 {
				return -u
			}
			return u
		}

		f0, f1 := fn[i-1], fn[i]
		if s.baseSampler.SampleFloat()*(f0-f1) < math.Exp(-0.5*x*x)-f1 {
			return x
		}
	}
}

// Sample returns a int64 value sampled from rounded gaussian distribution
// with given standard deviation.
// For large standard deviation, it heuristically samples the lower bits as uniform.
// If the output is not representable in int64, it overflows.
//
// Panics when stdDev is negative.
func (s *RoundedGaussianSampler) Sample(stdDev float64) int64 {
	if stdDev < 0 {
		panic("stdDev must be nonnegative")
	}

	if stdDev <= math.Exp2(roundPrec) {
		return int64(math.Round(s.normFloat() * stdDev))
	}

	stdDevMant, stdDevExp := math.Frexp(stdDev)
	logBound := stdDevExp - roundPrec
	if stdDevMant == 0.5 {
		logBound--
	}
	if logBound >= 63 {
		panic("stdDev too large")
	}

	bound := int64(1 << logBound)

	hi := int64(math.Round(s.normFloat()*math.Ldexp(stdDev, -logBound) + math.Ldexp(0.5, -logBound)))
	lo := s.baseSampler.SampleN(-bound/2, bound/2)

	return hi<<logBound + lo
}

// SampleBig returns a [big.Int] value sampled from an approximate rounded gaussian distribution
// with given standard deviation.
// For large standard deviation, it heuristically samples the lower bits as uniform.
//
// Panics when stdDev is negative.
func (s *RoundedGaussianSampler) SampleBig(stdDev float64) *big.Int {
	if stdDev < 0 || math.IsNaN(stdDev) {
		panic("stdDev must be nonnegative")
	}

	if stdDev <= math.Exp2(roundPrec) {
		return big.NewInt(s.Sample(stdDev))
	}

	stdDevMant, stdDevExp := math.Frexp(stdDev)
	logBound := stdDevExp - roundPrec
	if stdDevMant == 0.5 {
		logBound--
	}

	hi := big.NewInt(int64(math.Round(s.normFloat()*math.Ldexp(stdDev, -logBound) + math.Ldexp(0.5, -logBound))))

	loBytes := make([]byte, num.DivCeil(logBound, 8))
	s.baseSampler.Read(loBytes)
	loBytes[0] &= byte(0xff >> (len(loBytes)*8 - logBound))
	lo := new(big.Int).SetBytes(loBytes)

	boundHalf := big.NewInt(1)
	boundHalf.Lsh(boundHalf, uint(logBound-1))
	lo.Sub(lo, boundHalf)

	return hi.Add(hi.Lsh(hi, uint(logBound)), lo)
}
