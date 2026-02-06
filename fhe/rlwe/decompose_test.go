package rlwe_test

import (
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = csprng.NewUniformSamplerWithSeed(nil)
)

func Recompose(d rlwe.Decomposer, dcmp *rlwe.Tensor) *crt.Element {
	auxMod := d.Params().AuxModulus()
	modLen := dcmp.ModLen() - len(auxMod)

	pSum := crt.NewPolyCustom(d.Params().RingParams().Rank(), modLen+len(auxMod), false)
	pOut := crt.NewPolyCustom(d.Params().RingParams().Rank(), modLen, false)

	currMod := d.Params().Modulus()[:modLen]
	fullMod := d.Params().FullModulus()[:len(auxMod)+len(currMod)]
	eval := crt.NewOperator(d.Params().RingParams(), fullMod)
	for i := range dcmp.Value {
		g := d.GadgetVector()[i].WithModIdx(vec.Range(0, modLen+len(auxMod))...)
		eval.MulAddTo(pSum, dcmp.Value[i], g)
	}

	scaler := crt.NewScaler(currMod, fullMod)
	scaler.ScaleVecTo(pOut.Coeffs, pSum.Coeffs)

	return pOut
}

func TestDecompose(t *testing.T) {
	N := 1 << 10
	rP := dft.NewCyclotomicParameters(N << 1)
	q, qAux := rlwe.FindNTTPrimesFromBits(rP, 400, 100)

	paramsLiteral := rlwe.ParametersLiteral{
		RingParams: rP,
		Modulus:    q,
		AuxModulus: qAux,

		GadgetParams: rlwe.RNSGadgetParametersLiteral{
			ChunkSize: 2,
		},
		SecretKeyParams: crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		},
		NoiseParams: crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		},
	}

	t.Run("type=RNS", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.RNSGadgetParametersLiteral{
			ChunkSize: 2,
		}
		params := paramsLiteral.Compile()

		modLen := rSrc.SampleN(uint64(len(q)-1)) + 1
		us := crt.UniformSamplerParameters{}.Sampler()

		dcmp := rlwe.NewDecomposer(params)
		p := us.Sample(rP.Rank(), q[:modLen])
		pDec := dcmp.Decompose(p)
		pRef := Recompose(dcmp, pDec)

		assert.Equal(t, pRef, p)
	})

	t.Run("type=Digit", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.DigitGadgetParametersLiteral{
			LogDigitBase: 10,
		}
		params := paramsLiteral.Compile()

		modLen := rSrc.SampleN(uint64(len(q)-1)) + 1
		us := crt.UniformSamplerParameters{}.Sampler()

		dcmp := rlwe.NewDecomposer(params)
		p := us.Sample(rP.Rank(), q[:modLen])
		pDec := dcmp.Decompose(p)
		pRef := Recompose(dcmp, pDec)

		assert.Equal(t, pRef, p)
	})
}
