package csprng_test

import (
	"math"
	"math/big"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

func entropy(v []float64) float64 {
	hist := make(map[float64]uint64)
	for i := range v {
		hist[v[i]] += 1
	}

	var e float64
	for k := range hist {
		p := float64(hist[k]) / float64(len(v))
		e += -p * math.Log2(p)
	}
	return e
}

func meanStdDev(v []float64) (mean, stdDev float64) {
	N := float64(len(v))
	var sum float64
	for i := range v {
		sum += v[i]
	}
	mean = sum / N

	var vari float64
	for i := range v {
		vari += (v[i] - mean) * (v[i] - mean)
	}
	stdDev = math.Sqrt(vari / N)

	return
}

func TestUniform(t *testing.T) {
	s := csprng.NewUniformSamplerWithSeed(nil)

	v := make([]uint64, 1<<10)
	for i := range v {
		v[i] = s.Sample()
	}

	assert.InDelta(t, num.Log2(len(v)), entropy(vec.Cast[float64](v)), 0.1)
}

func TestRoundedGaussian(t *testing.T) {
	t.Run("float64", func(t *testing.T) {
		s := csprng.NewRoundedGaussianSamplerWithSeed(nil)
		centerRef := 8.0
		stdDevRef := 16.0

		v := make([]int64, 1<<10)
		for i := range v {
			v[i] = s.Sample(centerRef, stdDevRef)
		}

		mean, stdDev := meanStdDev(vec.Cast[float64](v))
		delta := stdDev / math.Sqrt(float64(len(v)))

		assert.InDelta(t, centerRef, mean, delta)
		assert.InDelta(t, stdDevRef, stdDev, delta)
	})

	t.Run("big.Float", func(t *testing.T) {
		s := csprng.NewRoundedGaussianSamplerWithSeed(nil)
		centerRef := 8.0
		stdDevRef := 16.0
		centerRefBig := big.NewFloat(centerRef)
		stdDevRefBig := big.NewFloat(stdDevRef)

		v := make([]int64, 1<<10)
		for i := range v {
			v[i] = s.SampleBig(centerRefBig, stdDevRefBig).Int64()
		}

		mean, stdDev := meanStdDev(vec.Cast[float64](v))
		delta := stdDev / math.Sqrt(float64(len(v)))

		assert.InDelta(t, centerRef, mean, delta)
		assert.InDelta(t, stdDevRef, stdDev, delta)
	})
}
