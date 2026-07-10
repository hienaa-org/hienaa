package rlwe_test

import (
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	skParams = crt.TernarySamplerParameters{
		Positive: 1.0 / 3.0,
		Negative: 1.0 / 3.0,
	}
	noiseParams = crt.RoundedGaussianSamplerParameters[float64]{
		Center: 0,
		StdDev: 3.2,
	}
	baseBits = float64(200)
	auxBits  = float64(100)
)

func testOperator(t *testing.T, rP dft.RingParameters) {
	baseMod, auxMod := rlwe.FindNTTPrimes(rP, baseBits, auxBits)
	params := rlwe.ParametersLiteral{
		RingParams:  rP,
		BaseModulus: baseMod,
		AuxModulus:  auxMod,

		GadgetParams: rlwe.RNSGadgetParametersLiteral{
			ChunkSize: 1,
		},
		SecretKeyParams: skParams,
		NoiseParams:     noiseParams,
	}.Compile()

	enc := rlwe.NewEncryptor(params)
	op := rlwe.NewOperator(params)
	pOp := op.PlainOperator()

	t.Run("Add", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt0, baseMod[:baseLen])
		randomElementTo(pt1, baseMod[:baseLen])

		ptRef := pOp.Add(pt0, pt1)

		ct0 := enc.Encrypt(pt0, true)
		ct1 := enc.Encrypt(pt1, true)

		ctOut := op.Add(ct0, ct1)
		ptOut := enc.Phase(ctOut)

		diff := pOp.Sub(ptOut, ptRef)

		noiseBound := 2 * params.NoiseParams().Bound()
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("AddElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			s := rlwe.NewElement(1, baseLen, 0, false)
			randomElementTo(pt, baseMod[:baseLen])
			randomElementTo(s, baseMod[:baseLen])

			ptRef := pOp.Add(pt, s)

			ct := enc.Encrypt(pt, true)

			ctOut := op.AddElement(ct, s)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt0, baseMod[:baseLen])
			randomElementTo(pt1, baseMod[:baseLen])

			ptRef := pOp.Add(pt0, pt1)

			ct := enc.Encrypt(pt0, false)
			ctOut := op.AddElement(ct, pt1)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("SubElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			s := rlwe.NewElement(1, baseLen, 0, false)
			randomElementTo(pt, baseMod[:baseLen])
			randomElementTo(s, baseMod[:baseLen])

			ptRef := pOp.Sub(pt, s)

			ct := enc.Encrypt(pt, true)

			ctOut := op.SubElement(ct, s)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt0, baseMod[:baseLen])
			randomElementTo(pt1, baseMod[:baseLen])

			ptRef := pOp.Sub(pt0, pt1)

			ct := enc.Encrypt(pt0, false)
			ctOut := op.SubElement(ct, pt1)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("Neg", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt, baseMod[:baseLen])
		ptRef := pOp.Neg(pt)

		ct := enc.Encrypt(pt, true)

		ctOut := op.Neg(ct)
		ptOut := enc.Phase(ctOut)
		diff := pOp.Sub(ptOut, ptRef)

		noiseBound := params.NoiseParams().Bound()
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("MulElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt, baseMod[:baseLen])
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			ptRef := pOp.Mul(pt, s)

			ct := enc.Encrypt(pt, true)

			ctOut := op.MulElement(ct, s)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := 2 * params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("PolyMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			ts := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt0, baseMod[:baseLen])
			pOp.FwdNTTTo(pt0, pt0)

			pt1 := rlwe.NewElementFrom(ts.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt1, pt1)

			ptRef := pOp.Mul(pt0, pt1)
			pOp.InvNTTTo(ptRef, ptRef)

			ct := enc.Encrypt(pt0, true)

			ctOut := op.MulElement(ct, pt1)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound() * float64(rP.ExpandFactor())
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("MulAddElement", func(t *testing.T) {
		t.Run("ScalarMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt0, baseMod[:baseLen])
			randomElementTo(pt1, baseMod[:baseLen])

			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			ptRef := pt0.Copy()
			pOp.MulAddTo(ptRef, pt1, s)

			ct0 := enc.Encrypt(pt0, true)
			ct1 := enc.Encrypt(pt1, true)

			op.MulAddElementTo(ct0, ct1, s)
			ptOut := enc.Phase(ct0)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := 3 * params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("PolyMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			ts := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt0, baseMod[:baseLen])
			randomElementTo(pt1, baseMod[:baseLen])
			pOp.FwdNTTTo(pt0, pt0)
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(ts.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt2, pt2)

			ptRef := pt0.Copy()
			pOp.MulAddTo(ptRef, pt2, pt1)
			pOp.InvNTTTo(ptRef, ptRef)

			ct0 := enc.Encrypt(pt0, true)
			ct1 := enc.Encrypt(pt1, true)

			op.MulAddElementTo(ct0, ct1, pt2)
			ptOut := enc.Phase(ct0)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound() * float64(rP.ExpandFactor()+1)
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("DivByAux", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
		auxLen := len(auxMod)

		pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
		randomElementTo(pt, params.FullModulus()[:baseLen+auxLen])

		ptRef := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		sc := crt.NewVecScaler(baseMod[:baseLen], params.FullModulus()[:baseLen+auxLen])
		sc.ScaleTo(ptRef.Value.Coeffs, pt.Value.Coeffs)

		ct := enc.Encrypt(pt, true)
		ctOut := op.DivByAuxModulus(ct, true)
		ptOut := enc.Phase(ctOut)
		diff := pOp.Sub(ptOut, ptRef)

		// Assume that the auxiliary modulus is large enough.
		noiseBound := float64(rP.ExpandFactor()) + 1
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("Scale", func(t *testing.T) {
		t.Run("l < oldLen", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-3)) + 3)
			newLen := int(rSrc.SampleN(uint64(oldLen-2)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			randomElementTo(pt, baseMod[:oldLen])

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewVecScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value.Coeffs, pt.Value.Coeffs)

			ct := enc.Encrypt(pt, true)

			ctOut := op.Scale(ct, newLen, false)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := float64(rP.ExpandFactor()) + 1
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("l > oldLen", func(t *testing.T) {
			newLen := int(rSrc.SampleN(uint64(len(baseMod)-3)) + 3)
			oldLen := int(rSrc.SampleN(uint64(newLen-2)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			randomElementTo(pt, baseMod[:oldLen])

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewVecScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value.Coeffs, pt.Value.Coeffs)

			ct := enc.Encrypt(pt, true)

			ctOut := op.Scale(ct, newLen, false)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound()
			for i := oldLen; i < newLen; i++ {
				noiseBound *= float64(baseMod[i].Value())
			}
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("l == oldLen", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-2)) + 2)
			newLen := oldLen

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			randomElementTo(pt, baseMod[:oldLen])

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewVecScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value.Coeffs, pt.Value.Coeffs)

			ct := enc.Encrypt(pt, true)

			ctOut := op.Scale(ct, newLen, false)
			ptOut := enc.Phase(ctOut)
			diff := pOp.Sub(ptOut, ptRef)

			noiseBound := params.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("GadgetProduct", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
		auxLen := len(auxMod)

		ts := crt.TernarySamplerParameters{
			Positive:      float64(1) / float64(3),
			Negative:      float64(1) / float64(3),
			HammingWeight: 0,
		}.Sampler()

		pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt0, baseMod[:baseLen])
		pOp.FwdNTTTo(pt0, pt0)

		pt1 := rlwe.NewElementFrom(ts.Sample(rP.Rank(), params.FullModulus()), auxLen)
		pOp.FwdNTTTo(pt1, pt1)

		pt1Mod := pt1.WithModLen(baseLen, 0)
		ptRef := pOp.Mul(pt0, pt1Mod)
		pOp.InvNTTTo(ptRef, ptRef)

		ctGad := enc.GadgetEncrypt(pt1, true)

		ctOut := op.GadgetProd(pt0, ctGad, true)
		ptOut := enc.Phase(ctOut)
		diff := pOp.Sub(ptOut, ptRef)

		noiseBound := float64(rP.ExpandFactor() + 1)
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("Relinearisation", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		ct0 := enc.Encrypt(pt, true)
		ct1 := enc.Encrypt(pt, true)
		rlk := enc.NewRelinKey()

		ctVec := rlwe.NewVectorCustom(rP.Rank(), baseLen, 0, 3, false)

		pOp.MulTo(ctVec.Value[0], ct0.Body, ct1.Body)
		pOp.MulTo(ctVec.Value[1], ct0.Mask, ct1.Body)
		pOp.MulAddTo(ctVec.Value[1], ct0.Body, ct1.Mask)
		pOp.MulTo(ctVec.Value[2], ct0.Mask, ct1.Mask)

		ctOut := op.Relin(ctVec, rlk, false)
		ptOut := enc.Phase(ctOut)

		noiseBound := params.NoiseParams().Bound()*params.NoiseParams().Bound()*float64(rP.ExpandFactor()) + float64(rP.ExpandFactor()) + 1
		assert.True(t, checkBound(pOp.AsBig(ptOut), noiseBound))
	})

	t.Run("KeySwitch", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		encIn := rlwe.NewEncryptor(params)
		ksk := enc.NewKeySwitchKey(encIn.SecretKey())

		pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		ct := encIn.Encrypt(pt, true)

		ctOut := op.KeySwitch(ct, ksk, true)
		ptOut := enc.Phase(ctOut)

		noiseBound := float64(rP.ExpandFactor() + 1)
		assert.True(t, checkBound(pOp.AsBig(ptOut), noiseBound))
	})

	t.Run("Automorphism", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		var idx int
		for {
			idx = int(rSrc.SampleN(uint64(rP.CycloIndex() - 1)))
			if pOp.CanAut(idx) {
				break
			}
		}
		atk := enc.NewAutomorphismKey(idx)

		pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt, baseMod[:baseLen])
		ptRef := pOp.Aut(pt, atk.Idx)

		ct := enc.Encrypt(pt, true)

		ctOut := op.Aut(ct, atk, true)
		ptOut := enc.Phase(ctOut)
		diff := pOp.Sub(ptOut, ptRef)

		noiseBound := float64(rP.ExpandFactor() + 1)
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("ExternalProduct", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
		auxLen := len(auxMod)

		ts := crt.TernarySamplerParameters{
			Positive:      float64(1) / float64(3),
			Negative:      float64(1) / float64(3),
			HammingWeight: 0,
		}.Sampler()

		pt0 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt0, baseMod[:baseLen])
		pOp.FwdNTTTo(pt0, pt0)

		pt1 := rlwe.NewElementFrom(ts.Sample(rP.Rank(), params.FullModulus()), auxLen)
		pOp.FwdNTTTo(pt1, pt1)

		pt1Mod := pt1.WithModLen(baseLen, 0)
		ptRef := pOp.Mul(pt0, pt1Mod)
		pOp.InvNTTTo(ptRef, ptRef)

		ct := enc.Encrypt(pt0, true)
		ctRGSW := enc.RGSWEncrypt(pt1, true)

		ctOut := op.ExtProd(ct, ctRGSW, true)
		ptOut := enc.Phase(ctOut)
		diff := pOp.Sub(ptOut, ptRef)

		noiseBound := params.NoiseParams().Bound() * float64(rP.ExpandFactor()+1)
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})
}

func TestOperator(t *testing.T) {
	t.Run("type=CyclotomicPow2", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 11)
		testOperator(t, rP)
	})

	t.Run("type=CyclotomicAny", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))
		testOperator(t, rP)
	})

	t.Run("type=AutFixedPow2", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<13, 1<<11)
		testOperator(t, rP)
	})

	t.Run("type=AutFixedPrime", func(t *testing.T) {
		M := num.MustNextPrime(int(rSrc.SampleN(1<<12)), 1)
		primes, exps := num.Factor(M - 1)
		fold := 1
		for i := range primes {
			e := int(rSrc.SampleN(uint64(exps[i])))
			for j := 0; j < e; j++ {
				fold *= primes[i]
			}
		}
		N := (M - 1) / fold
		rP := dft.NewAutFixedParameters(M, N)

		testOperator(t, rP)
	})
}
