package crt_test

import (
	"math"
	"math/big"
	"slices"
	"testing"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

func entropy(v []*big.Int) float64 {
	hist := make(map[*big.Int]uint64)
	for i := range v {
		for k := range hist {
			if k.Cmp(v[i]) == 0 {
				hist[k] += 1
				goto next
			}
		}
		hist[new(big.Int).Set(v[i])] = 1
	next:
	}

	var e float64
	for k := range hist {
		p := float64(hist[k]) / float64(len(v))
		e += -p * math.Log2(p)
	}
	return e
}

// func meanStdDev(v []*big.Int) (mean, stdDev float64) {
// 	N := big.NewFloat(0).SetInt64(int64(len(v)))
// 	sumBig := big.NewFloat(0)
// 	for i := range v {
// 		sumBig.Add(sumBig, big.NewFloat(0).SetInt(v[i]))
// 	}
// 	meanBig := sumBig.Quo(sumBig, N)

// 	variBig := big.NewFloat(0)
// 	for i := range v {
// 		diff := big.NewFloat(0).SetInt(v[i])
// 		diff.Sub(diff, meanBig)
// 		diff.Mul(diff, diff)
// 		variBig.Add(variBig, diff)
// 	}
// 	variBig.Quo(variBig, N)
// 	stdDevBig := big.NewFloat(0).Sqrt(variBig)

// 	mean, _ = meanBig.Float64()
// 	stdDev, _ = stdDevBig.Float64()

// 	return
// }

func TestSampler(t *testing.T) {
	rP := dft.NewCyclicParameters(1 << 10)
	q := dft.MustFindPrevNTTPrimes(rP, num.MaxModulusBits, 2)
	pev := crt.NewPolyEvaluator(rP, q)

	t.Run("type=Uniform", func(t *testing.T) {
		s := crt.UniformSamplerParameters{
			Params:  rP,
			Modulus: q,
		}.Sampler()

		pOut := s.Sample()
		vOut := pev.AsBig(pOut)

		assert.InDelta(t, num.Log2(len(vOut)), entropy(vOut), 0.1)
	})

	t.Run("type=BoundedUniform", func(t *testing.T) {
		boundMin := big.NewInt(-10)
		boundMax := big.NewInt(10)
		s := crt.UniformSamplerParameters{
			Params:  rP,
			Modulus: q,

			BoundMin: boundMin,
			BoundMax: boundMax,
		}.Sampler()

		pOut := s.Sample()
		vOut := pev.AsBig(pOut)

		assert.True(t, slices.MinFunc(vOut, (*big.Int).Cmp).Cmp(boundMin) >= 0)
		assert.True(t, slices.MaxFunc(vOut, (*big.Int).Cmp).Cmp(boundMax) < 0)
	})

	t.Run("type=Ternary", func(t *testing.T) {
		s := crt.TernarySamplerParameters{
			Params:  rP,
			Modulus: q,

			Positive: 1.0 / 3.0,
			Negative: 1.0 / 3.0,
		}.Sampler()

		pOut := s.Sample()
		vOut := pev.AsBig(pOut)

		assert.InDelta(t, entropy(vOut), math.Log2(3), 0.1)
		assert.True(t, slices.MaxFunc(vOut, (*big.Int).Cmp).Cmp(big.NewInt(1)) <= 0)
		assert.True(t, slices.MinFunc(vOut, (*big.Int).Cmp).Cmp(big.NewInt(-1)) >= 0)
	})

	t.Run("type=TernaryFixedHammingWeight", func(t *testing.T) {
		hwRef := rP.Rank() / 4
		s := crt.TernarySamplerParameters{
			Params:  rP,
			Modulus: q,

			Positive:      1.0 / 3.0,
			Negative:      1.0 / 3.0,
			HammingWeight: hwRef,
		}.Sampler()

		pOut := s.Sample()
		vOut := pev.AsBig(pOut)

		hw := 0
		for i := range vOut {
			if vOut[i].Cmp(big.NewInt(0)) != 0 {
				hw++
			}
		}

		assert.Equal(t, hwRef, hw)
		assert.True(t, slices.MaxFunc(vOut, (*big.Int).Cmp).Cmp(big.NewInt(1)) <= 0)
		assert.True(t, slices.MinFunc(vOut, (*big.Int).Cmp).Cmp(big.NewInt(-1)) >= 0)
	})
}
