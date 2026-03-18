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
	p := rlwe.ParametersLiteral{
		RingParams:  rP,
		BaseModulus: baseMod,
		AuxModulus:  auxMod,

		GadgetParams: rlwe.RNSGadgetParametersLiteral{
			ChunkSize: 1,
		},
		SecretKeyParams: skParams,
		NoiseParams:     noiseParams,
	}.Compile()

	e := rlwe.NewEncryptor(p)
	o := rlwe.NewOperator(p)
	pOp := o.PlainOperator()

	t.Run("Add", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt1, baseMod[:baseLen])
		randomElementTo(pt2, baseMod[:baseLen])

		pRef := pOp.Add(pt1, pt2)

		c1 := e.Encrypt(pt1, true)
		c2 := e.Encrypt(pt2, true)

		cOut := o.Add(c1, c2)
		pOut := e.Phase(cOut)

		diff := pOp.Sub(pOut, pRef)

		noiseBound := 2 * noiseParams.Bound()
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

			c := e.Encrypt(pt, true)

			cOut := o.AddElement(c, s)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := p.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt1, baseMod[:baseLen])
			randomElementTo(pt2, baseMod[:baseLen])

			ptRef := pOp.Add(pt1, pt2)

			c := e.Encrypt(pt1, false)
			cOut := o.AddElement(c, pt2)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := p.NoiseParams().Bound()
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

			c := e.Encrypt(pt, true)

			cOut := o.SubElement(c, s)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := p.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt1, baseMod[:baseLen])
			randomElementTo(pt2, baseMod[:baseLen])

			ptRef := pOp.Sub(pt1, pt2)

			c := e.Encrypt(pt1, false)
			cOut := o.SubElement(c, pt2)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := p.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("Neg", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt, baseMod[:baseLen])
		ptRef := pOp.Neg(pt)

		c := e.Encrypt(pt, true)

		cOut := o.Neg(c)
		pOut := e.Phase(cOut)
		diff := pOp.Sub(pOut, ptRef)

		noiseBound := p.NoiseParams().Bound()
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("MulElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt, baseMod[:baseLen])
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pOp.Mul(pt, s)

			c := e.Encrypt(pt, true)

			cOut := o.MulElement(c, s)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, pRef)

			noiseBound := 2 * p.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("PolyMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			sP := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt1, baseMod[:baseLen])
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt2, pt2)

			pRef := pOp.Mul(pt1, pt2)
			pOp.InvNTTTo(pRef, pRef)

			c := e.Encrypt(pt1, true)

			cOut := o.MulElement(c, pt2)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, pRef)

			noiseBound := p.NoiseParams().Bound() * float64(rP.ExpFactor())
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("MulAddElement", func(t *testing.T) {
		t.Run("ScalarMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt1, baseMod[:baseLen])
			randomElementTo(pt2, baseMod[:baseLen])

			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt2, s)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddElementTo(c1, c2, s)
			pOut := e.Phase(c1)
			diff := pOp.Sub(pOut, pRef)

			noiseBound := 3 * p.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("PolyMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			sP := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()
			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt1, baseMod[:baseLen])
			randomElementTo(pt2, baseMod[:baseLen])
			pOp.FwdNTTTo(pt1, pt1)
			pOp.FwdNTTTo(pt2, pt2)

			pt3 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt3, pt3)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt3, pt2)
			pOp.InvNTTTo(pRef, pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddElementTo(c1, c2, pt3)
			pOut := e.Phase(c1)
			diff := pOp.Sub(pOut, pRef)

			noiseBound := p.NoiseParams().Bound() * float64(rP.ExpFactor()+1)
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("DivByAux", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
		auxLen := len(auxMod)

		pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
		randomElementTo(pt, p.FullModulus()[:baseLen+auxLen])

		ptRef := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		sc := crt.NewScaler(baseMod[:baseLen], p.FullModulus()[:baseLen+auxLen])
		sc.ScaleTo(ptRef.Value, pt.Value)

		c := e.Encrypt(pt, true)
		cOut := o.DivByAuxModulus(c, true)
		pOut := e.Phase(cOut)
		diff := pOp.Sub(pOut, ptRef)

		// Assume that the auxiliary modulus is large enough.
		noiseBound := float64(rP.ExpFactor()) + 1
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("Scale", func(t *testing.T) {
		t.Run("l < oldLen", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-3)) + 3)
			newLen := int(rSrc.SampleN(uint64(oldLen-2)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			randomElementTo(pt, baseMod[:oldLen])

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := float64(rP.ExpFactor()) + 1
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("l > oldLen", func(t *testing.T) {
			newLen := int(rSrc.SampleN(uint64(len(baseMod)-3)) + 3)
			oldLen := int(rSrc.SampleN(uint64(newLen-2)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			randomElementTo(pt, baseMod[:oldLen])

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := p.NoiseParams().Bound()
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
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			diff := pOp.Sub(pOut, ptRef)

			noiseBound := p.NoiseParams().Bound()
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("GadgetProduct", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
		auxLen := len(auxMod)

		s := crt.TernarySamplerParameters{
			Positive:      float64(1) / float64(3),
			Negative:      float64(1) / float64(3),
			HammingWeight: 0,
		}.Sampler()

		pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(pt1, baseMod[:baseLen])
		pOp.FwdNTTTo(pt1, pt1)

		pt2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
		pOp.FwdNTTTo(pt2, pt2)

		pt2Mod := pt2.WithModLen(baseLen, 0)
		ptRef := pOp.Mul(pt1, pt2Mod)
		pOp.InvNTTTo(ptRef, ptRef)

		c := e.GadgetEncrypt(pt2, true)

		cOut := o.GadgetProd(pt1, c, true)
		pOut := e.Phase(cOut)
		diff := pOp.Sub(pOut, ptRef)

		noiseBound := float64(rP.ExpFactor() + 1)
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("Relinearisation", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		c1 := e.Encrypt(pt, true)
		c2 := e.Encrypt(pt, true)
		rlk := e.NewRelinKey()

		vec := rlwe.NewVectorCustom(rP.Rank(), baseLen, 0, 3, false)

		pOp.MulTo(vec.Value[0], c1.Body, c2.Body)
		pOp.MulTo(vec.Value[1], c1.Mask, c2.Body)
		pOp.MulAddTo(vec.Value[1], c1.Body, c2.Mask)
		pOp.MulTo(vec.Value[2], c1.Mask, c2.Mask)

		cOut := o.Relin(vec, rlk, false)
		pOut := e.Phase(cOut)

		noiseBound := p.NoiseParams().Bound()*p.NoiseParams().Bound()*float64(rP.ExpFactor()) + float64(rP.ExpFactor()) + 1
		assert.True(t, checkBound(pOp.AsBig(pOut), noiseBound))
	})

	t.Run("KeySwitch", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		eNew := rlwe.NewEncryptor(p)
		ksk := e.NewKeySwitchKey(eNew.SecretKey())

		p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		c := eNew.Encrypt(p, true)

		cOut := o.KeySwitch(c, ksk, true)
		pOut := e.Phase(cOut)

		noiseBound := float64(rP.ExpFactor() + 1)
		assert.True(t, checkBound(pOp.AsBig(pOut), noiseBound))
	})

	t.Run("Automorphism", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		var idx int
		for {
			idx = int(rSrc.SampleN(uint64(rP.CycloOrder() - 1)))
			if pOp.CanAut(idx) {
				break
			}
		}
		atk := e.NewAutomorphismKey(idx)

		p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(p, baseMod[:baseLen])
		pRef := pOp.Aut(p, atk.Idx)

		c := e.Encrypt(p, true)

		cOut := o.Aut(c, atk, true)
		pOut := e.Phase(cOut)
		diff := pOp.Sub(pOut, pRef)

		noiseBound := float64(rP.ExpFactor() + 1)
		assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
	})

	t.Run("ExternalProduct", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
		auxLen := len(auxMod)

		s := crt.TernarySamplerParameters{
			Positive:      float64(1) / float64(3),
			Negative:      float64(1) / float64(3),
			HammingWeight: 0,
		}.Sampler()

		p1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
		randomElementTo(p1, baseMod[:baseLen])
		pOp.FwdNTTTo(p1, p1)

		p2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
		pOp.FwdNTTTo(p2, p2)

		p2Mod := p2.WithModLen(baseLen, 0)
		pRef := pOp.Mul(p1, p2Mod)
		pOp.InvNTTTo(pRef, pRef)

		c := e.Encrypt(p1, true)
		r := e.RGSWEncrypt(p2, true)

		cOut := o.ExtProd(c, r, true)
		pOut := e.Phase(cOut)
		diff := pOp.Sub(pOut, pRef)

		noiseBound := p.NoiseParams().Bound() * float64(rP.ExpFactor()+1)
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
