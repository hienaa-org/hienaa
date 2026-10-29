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
	M := rSrc.SampleN(4, 1<<10)
	ringParams := dft.NewCyclotomicParameters(M)
	q, qAux := rlwe.FindNTTPrimes(ringParams, float64(rSrc.SampleN(100, 800)), float64(rSrc.SampleN(100, 800)))

	paramsLiteral := rlwe.ParametersLiteral{
		RingParams:  ringParams,
		BaseModulus: q,
		AuxModulus:  qAux,

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

	us := crt.UniformSamplerParameters{}.Sampler()

	t.Run("type=RNS", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.RNSGadgetParametersLiteral{
			ChunkSize: rSrc.SampleN(1, len(q)),
		}
		params := paramsLiteral.Compile()

		dcmp := rlwe.NewDecomposer(params)

		baseLen := rSrc.SampleN(1, len(q))
		baseOp := params.Op().WithModSlice(len(params.AuxMod()), len(params.AuxMod())+baseLen)

		p := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
		us.SampleTo(p.Value, q[:baseLen])
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

		baseLen := rSrc.SampleN(1, len(q))
		baseOp := params.Op().WithModSlice(len(params.AuxMod()), len(params.AuxMod())+baseLen)

		p := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
		us.SampleTo(p.Value, q[:baseLen])
		pDcmp := dcmp.Decompose(p, dft.FormNTT)
		pRcmp := recompose(params, pDcmp, params.GadgetVec())

		baseOp.FwdNTTTo(p.Value, p.Value)
		assert.Equal(t, p, pRcmp)
	})
}
