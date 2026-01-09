package rlwe_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{12, 13, 14, 15}
)

func Recompose(lweP rlwe.Parameters, d rlwe.Decomposer, pOut *rlwe.Tensor) *crt.Poly {
	rP := lweP.RingParams()
	q := lweP.Modulus()
	auxQ := lweP.AuxModulus()
	modLen := pOut.ModLen() - len(auxQ)

	sum := crt.NewPolyCustom(rP.Rank(), modLen+len(auxQ), false)
	res := crt.NewPolyCustom(rP.Rank(), modLen, false)

	nowQ := q[:modLen]
	gadgetMod := append(auxQ, nowQ...)
	eval := crt.NewPolyEvaluator(rP, gadgetMod)
	for i := range pOut.Value {
		gi := d.GadgetVector()[i][:modLen+len(auxQ)]
		eval.ScalarMulAddTo(sum, pOut.Value[i], gi)
	}

	sc := crt.NewScaler(nowQ, gadgetMod)
	sc.ScaleVecTo(res.Coeffs, sum.Coeffs)

	return res
}

func TestDecompose(t *testing.T) {
	t.Run("type=RNS", func(t *testing.T) {
		N := 1 << 10
		rP := dft.NewCyclotomicParameters(N << 1)

		q, auxQ := rlwe.FindNTTPrimesFromBits(rP, 400, 100)
		gP := rlwe.NewRNSGadgetParameters(2)
		sP := crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		}

		lweP := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    q,
			AuxModulus: auxQ,

			GadgetParams:    gP,
			SecretKeyParams: sP,
			NoiseParams:     sP,
		}.Compile()

		modLen := rSrc.SampleN(uint64(len(q)-1)) + 1
		s := crt.UniformSamplerParameters{}.Sampler()

		d := rlwe.NewDecomposer(lweP)
		p := s.Sample(rP.Rank(), q[:modLen])
		pDec := d.Decompose(p)
		pRef := Recompose(lweP, d, pDec)

		assert.Equal(t, pRef, p)
	})

	t.Run("type=Digit", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 3)
		q, auxQ := rlwe.FindNTTPrimesFromBits(rP, 400, 50)
		gP := rlwe.NewDigitGadgetParameters(10)
		sP := crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		}

		lweP := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    q,
			AuxModulus: auxQ,

			GadgetParams:    gP,
			SecretKeyParams: sP,
			NoiseParams:     sP,
		}.Compile()

		modLen := rSrc.SampleN(uint64(len(q)-1)) + 1
		s := crt.UniformSamplerParameters{}.Sampler()

		d := rlwe.NewDecomposer(lweP)
		p := s.Sample(rP.Rank(), q[:modLen])
		pDec := d.Decompose(p)
		pRef := Recompose(lweP, d, pDec)

		assert.Equal(t, pRef, p)
	})
}

func BenchmarkDecompose(b *testing.B) {
	b.Run("type=RNS", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(N << 1)

			q, auxQ := rlwe.FindNTTPrimesFromBits(rP, 400, 100)
			gP := rlwe.NewRNSGadgetParameters(2)
			sP := crt.TernarySamplerParameters{
				Positive: 0.25,
				Negative: 0.25,
			}

			lweP := rlwe.ParametersLiteral{
				RingParams: rP,
				Modulus:    q,
				AuxModulus: auxQ,

				GadgetParams:    gP,
				SecretKeyParams: sP,
				NoiseParams:     sP,
			}.Compile()

			s := crt.UniformSamplerParameters{}.Sampler()

			d := rlwe.NewDecomposer(lweP)
			p := s.Sample(rP.Rank(), q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				pDec := d.Decompose(p)
				_ = pDec
			})
		}
	})

	b.Run("type=Digit", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(N << 1)

			q, auxQ := rlwe.FindNTTPrimesFromBits(rP, 400, 50)
			gP := rlwe.NewDigitGadgetParameters(30)
			sP := crt.TernarySamplerParameters{
				Positive: 0.25,
				Negative: 0.25,
			}

			lweP := rlwe.ParametersLiteral{
				RingParams: rP,
				Modulus:    q,
				AuxModulus: auxQ,

				GadgetParams:    gP,
				SecretKeyParams: sP,
				NoiseParams:     sP,
			}.Compile()

			s := crt.UniformSamplerParameters{}.Sampler()

			d := rlwe.NewDecomposer(lweP)
			p := s.Sample(rP.Rank(), q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				pDec := d.Decompose(p)
				_ = pDec
			})
		}
	})
}
