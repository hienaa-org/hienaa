package rlwe_test

import (
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

func recompose(params rlwe.Parameters, pDcmp *rlwe.Vector, gadVec []*rlwe.Element) *rlwe.Element {
	fullOp := params.Op().WithModSlice(0, pDcmp.FullModLen())
	pSum := crt.NewPoly(params.Rank(), pDcmp.FullModLen(), dft.FormNTT)
	for i := range pDcmp.Value {
		fullOp.MulAddTo(pSum, pDcmp.Value[i].Value, gadVec[i].Value.WithModSlice(0, pDcmp.FullModLen()))
	}

	scaler := crt.NewScaler(fullOp.WithModSlice(pDcmp.AuxModLen(), pDcmp.FullModLen()), fullOp)
	return rlwe.NewElementFrom(scaler.Scale(pSum, dft.FormNTT), 0)
}

func TestDecompose(t *testing.T) {
	paramsLiteral := rlwe.ParametersLiteral{
		RingParams: dft.NewCyclotomicParameters(rSrc.SampleN(1<<2, 1<<10)),
		LogBaseMod: float64(rSrc.SampleN(1, 10) * num.MaxModulusBits),
		LogAuxMod:  float64(rSrc.SampleN(1, 10) * num.MaxModulusBits),

		SecretKeyParams: crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		},
		NoiseParams: crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		},
	}

	us := crt.UniformSamplerParameters{}.Sampler()

	t.Run("type=RNS", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.RNSGadgetParametersLiteral{
			ChunkSize: rSrc.SampleN(1, int(paramsLiteral.LogBaseMod/num.MaxModulusBits)+1),
		}
		params := paramsLiteral.Compile()

		dcmp := rlwe.NewDecomposer(params)

		baseLen := rSrc.SampleN(1, params.BaseModLen()+1)
		baseOp := params.Op().WithModSlice(params.AuxModLen(), params.AuxModLen()+baseLen)

		p := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
		us.SampleTo(p.Value, params.BaseMod()[:baseLen])
		pDcmp := dcmp.Decompose(p, dft.FormNTT)
		pRcmp := recompose(params, pDcmp, params.GadgetVec())

		baseOp.FwdNTTTo(p.Value, p.Value)
		assert.Equal(t, p, pRcmp)
	})

	t.Run("type=Digit", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.DigitGadgetParametersLiteral{
			LogDigitBase: rSrc.SampleN(1, int(num.MaxModulusBits)+1),
		}
		params := paramsLiteral.Compile()

		dcmp := rlwe.NewDecomposer(params)

		baseLen := rSrc.SampleN(1, params.BaseModLen()+1)
		baseOp := params.Op().WithModSlice(params.AuxModLen(), params.AuxModLen()+baseLen)

		p := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
		us.SampleTo(p.Value, params.BaseMod()[:baseLen])
		pDcmp := dcmp.Decompose(p, dft.FormNTT)
		pRcmp := recompose(params, pDcmp, params.GadgetVec())

		baseOp.FwdNTTTo(p.Value, p.Value)
		assert.Equal(t, p, pRcmp)
	})
}
