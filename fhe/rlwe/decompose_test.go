package rlwe_test

import (
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = csprng.NewUniformSamplerWithSeed(nil)
)

func Recompose(d rlwe.Decomposer, dcmp *rlwe.Vector) *crt.Element {
	auxMod := d.Params().AuxModulus()
	baseLen := dcmp.BaseModLen()
	auxLen := d.AuxModLen(baseLen)

	pSum := crt.NewPolyCustom(d.Params().RingParams().Rank(), baseLen+auxLen, false)
	pOut := crt.NewPolyCustom(d.Params().RingParams().Rank(), baseLen, false)

	currMod := d.Params().BaseModulus()[:baseLen]
	fullMod := d.Params().FullModulus()[:len(auxMod)+len(currMod)]
	op := crt.NewOperator(d.Params().RingParams(), fullMod)
	for i := range dcmp.Value {
		g := d.GadgetVector()[i].WithModLen(baseLen, auxLen)
		op.MulAddTo(pSum, dcmp.Value[i].Value, g.Value)
	}

	scaler := crt.NewVecScaler(currMod, fullMod)
	scaler.ScaleTo(pOut.Coeffs, pSum.Coeffs)

	return pOut
}

func TestDecompose(t *testing.T) {
	N := 1 << 10
	rP := dft.NewCyclotomicParameters(N << 1)
	q, qAux := rlwe.FindNTTPrimes(rP, 400, 100)

	paramsLiteral := rlwe.ParametersLiteral{
		RingParams:  rP,
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

	t.Run("type=RNS", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.RNSGadgetParametersLiteral{
			ChunkSize: 2,
		}
		params := paramsLiteral.Compile()

		dcmp := rlwe.NewDecomposer(params)

		baseLen := int(rSrc.SampleN(uint64(len(q)-1))) + 1
		us := crt.UniformSamplerParameters{}.Sampler()

		p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		us.SampleTo(p.Value, q[:baseLen])
		pDec := dcmp.Decompose(p, false)
		pRef := Recompose(dcmp, pDec)

		assert.Equal(t, pRef, p.Value)
	})

	t.Run("type=Digit", func(t *testing.T) {
		paramsLiteral := paramsLiteral
		paramsLiteral.GadgetParams = rlwe.DigitGadgetParametersLiteral{
			LogDigitBase: 10,
		}
		params := paramsLiteral.Compile()

		dcmp := rlwe.NewDecomposer(params)

		baseLen := int(rSrc.SampleN(uint64(len(q)-1))) + 1
		us := crt.UniformSamplerParameters{}.Sampler()

		p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		us.SampleTo(p.Value, q[:baseLen])
		pDec := dcmp.Decompose(p, false)
		pRef := Recompose(dcmp, pDec)

		assert.Equal(t, pRef, p.Value)
	})
}
