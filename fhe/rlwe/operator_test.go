package rlwe_test

import (
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = csprng.NewUniformSamplerWithSeed(nil)
)

func testOperator(t *testing.T, ringParams dft.RingParameters) {
	q, qAux := rlwe.FindNTTPrimes(ringParams, float64(rSrc.SampleN(100, 800)), float64(rSrc.SampleN(100, 800)))
	params := rlwe.ParametersLiteral{
		RingParams:  ringParams,
		BaseModulus: q,
		AuxModulus:  qAux,

		GadgetParams: rlwe.RNSGadgetParametersLiteral{
			ChunkSize: rSrc.SampleN(1, len(qAux)),
		},
		SecretKeyParams: crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		},
		NoiseParams: crt.TernarySamplerParameters{
			Positive: 0.25,
			Negative: 0.25,
		},
	}.Compile()

	sk := rlwe.NewSecretKey(params, true, dft.FormCoeff)
	params.SecretKeyParams().Sampler().SampleTo(sk.Value, params.FullMod())
	params.Op().FwdNTTTo(sk.Value, sk.Value)

	enc := rlwe.NewEncryptorWithKey(params, sk)
	op := rlwe.NewOperator(params)
	pOp := op.PlainOp()
	us := crt.UniformSamplerParameters{}.Sampler()
	ts := params.SecretKeyParams().Sampler()
	noiseBound := params.NoiseParams().Bound()
	expFac := expansionFactor(ringParams)

	baseLen := rSrc.SampleN(1, len(q)+1)
	baseCvt := crt.NewBigConverter(params.Op().WithModSlice(len(qAux), len(qAux)+baseLen))

	pt0 := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
	pt1 := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
	us.SampleTo(pt0.Value, q[:baseLen])
	us.SampleTo(pt1.Value, q[:baseLen])
	pt0NTT := pOp.FwdNTT(pt0)
	pt1NTT := pOp.FwdNTT(pt1)

	ct0 := enc.Encrypt[*rlwe.Ciphertext](pt0, dft.FormNTT)
	ct1 := enc.Encrypt[*rlwe.Ciphertext](pt1, dft.FormNTT)

	t.Run("Add", func(t *testing.T) {
		t.Run("Ciphertext", func(t *testing.T) {
			ctOut := op.Add(ct0, ct1)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Add(pt0, pt1)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), 2*noiseBound))
		})

		t.Run("Scalar", func(t *testing.T) {
			s := rlwe.NewScalarCustom(baseLen, 0)
			us.SampleTo(s.Value, q[:baseLen])

			ctOut := op.Add(ct0, s)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Add(pt0, s)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			ct := op.InvNTT(ct0)
			ctOut := op.Add(ct, pt1)
			assert.Equal(t, dft.FormCoeff, ctOut.Form())
			op.FwdNTTTo(ctOut, ctOut)

			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Add(pt0, pt1)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
		})
	})

	t.Run("Sub", func(t *testing.T) {
		t.Run("Ciphertext", func(t *testing.T) {
			ctOut := op.Sub(ct0, ct1)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Sub(pt0, pt1)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), 2*noiseBound))
		})

		t.Run("Scalar", func(t *testing.T) {
			s := rlwe.NewScalarCustom(baseLen, 0)
			us.SampleTo(s.Value, q[:baseLen])

			ctOut := op.Sub(ct0, s)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Sub(pt0, s)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			ct := op.InvNTT(ct0)
			ctOut := op.Sub(ct, pt1)
			assert.Equal(t, dft.FormCoeff, ctOut.Form())
			op.FwdNTTTo(ctOut, ctOut)

			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Sub(pt0, pt1)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
		})
	})

	t.Run("Neg", func(t *testing.T) {
		ctOut := op.Neg(ct0)
		ptOut := enc.Phase(ctOut)
		ptOutRef := pOp.Neg(pt0)

		diff := pOp.Sub(ptOut, ptOutRef)
		assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
	})

	t.Run("Mul", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, q[:baseLen]), 0)

			ctOut := op.Mul(ct0, s)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Mul(pt0, s)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), 2*noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			pt := rlwe.NewElementFrom(ts.Sample(params.Rank(), q[:baseLen]), 0)
			pOp.FwdNTTTo(pt, pt)

			ctOut := op.Mul(ct0, pt)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pOp.Mul(pt0NTT, pt)
			pOp.InvNTTTo(ptOutRef, ptOutRef)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), expFac*noiseBound))
		})
	})

	t.Run("MulAdd", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, q[:baseLen]), 0)

			ctOut := ct0.Copy()
			op.MulAddTo(ctOut, ct1, s)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pt0.Copy()
			pOp.MulAddTo(ptOutRef, pt1, s)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), 3*noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			pt := rlwe.NewElementFrom(ts.Sample(params.Rank(), q[:baseLen]), 0)
			pOp.FwdNTTTo(pt, pt)

			ctOut := ct0.Copy()
			op.MulAddTo(ctOut, ct1, pt)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pt0NTT.Copy()
			pOp.MulAddTo(ptOutRef, pt1NTT, pt)
			pOp.InvNTTTo(ptOutRef, ptOutRef)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), (expFac+1)*noiseBound))
		})
	})

	t.Run("MulSub", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, q[:baseLen]), 0)

			ctOut := ct0.Copy()
			op.MulSubTo(ctOut, ct1, s)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pt0.Copy()
			pOp.MulSubTo(ptOutRef, pt1, s)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), 3*noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			pt := rlwe.NewElementFrom(ts.Sample(params.Rank(), q[:baseLen]), 0)
			pOp.FwdNTTTo(pt, pt)

			ctOut := ct0.Copy()
			op.MulSubTo(ctOut, ct1, pt)
			ptOut := enc.Phase(ctOut)
			ptOutRef := pt0NTT.Copy()
			pOp.MulSubTo(ptOutRef, pt1NTT, pt)
			pOp.InvNTTTo(ptOutRef, ptOutRef)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), (expFac+1)*noiseBound))
		})
	})

	t.Run("Scale", func(t *testing.T) {
		t.Run("Aux", func(t *testing.T) {
			pt := rlwe.NewPolyCustom(params.Rank(), baseLen, len(qAux), dft.FormCoeff)
			us.SampleTo(pt.Value, params.FullMod()[:len(qAux)+baseLen])

			ct := enc.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
			ctOut := op.Scale(ct, baseLen, dft.FormNTT)
			ptOut := enc.Phase(ctOut)
			ptOutRef := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
			sc := crt.NewVecScaler(q[:baseLen], params.FullMod()[:len(qAux)+baseLen])
			sc.ScaleTo(ptOutRef.Value.Coeffs, pt.Value.Coeffs)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), expFac+1))
		})

		t.Run("Decrease", func(t *testing.T) {
			oldLen := rSrc.SampleN(2, len(q)+1)
			newLen := rSrc.SampleN(1, oldLen)
			cvt := crt.NewBigConverter(params.Op().WithModSlice(len(qAux), len(qAux)+newLen))

			pt := rlwe.NewPolyCustom(params.Rank(), oldLen, 0, dft.FormCoeff)
			us.SampleTo(pt.Value, q[:oldLen])

			ct := enc.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
			ctOut := op.Scale(ct, newLen, dft.FormCoeff)
			assert.Equal(t, dft.FormCoeff, ctOut.Form())
			op.FwdNTTTo(ctOut, ctOut)

			ptOut := enc.Phase(ctOut)
			ptOutRef := rlwe.NewPolyCustom(params.Rank(), newLen, 0, dft.FormCoeff)
			sc := crt.NewVecScaler(q[:newLen], q[:oldLen])
			sc.ScaleTo(ptOutRef.Value.Coeffs, pt.Value.Coeffs)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(cvt.Convert(diff.Value), expFac+1))
		})

		t.Run("Increase", func(t *testing.T) {
			newLen := rSrc.SampleN(2, len(q)+1)
			oldLen := rSrc.SampleN(1, newLen)

			pt := rlwe.NewPolyCustom(params.Rank(), oldLen, 0, dft.FormCoeff)
			us.SampleTo(pt.Value, q[:oldLen])

			ct := enc.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
			ctOut := op.Scale(ct, newLen, dft.FormNTT)
			ptOut := enc.Phase(ctOut)
			ptOutRef := rlwe.NewPolyCustom(params.Rank(), newLen, 0, dft.FormCoeff)
			sc := crt.NewVecScaler(q[:newLen], q[:oldLen])
			sc.ScaleTo(ptOutRef.Value.Coeffs, enc.Phase(ct).Value.Coeffs)

			assert.Equal(t, ptOutRef, ptOut)
		})

		t.Run("Same", func(t *testing.T) {
			ctOut := op.Scale(ct0, baseLen, dft.FormNTT)

			assert.Equal(t, ct0, ctOut)
		})
	})

	t.Run("GadProd", func(t *testing.T) {
		pt := rlwe.NewElementFrom(ts.Sample(params.Rank(), params.FullMod()), len(qAux))
		pOp.FwdNTTTo(pt, pt)
		ctGad := enc.Encrypt[*rlwe.GadgetEncryption](pt, dft.FormNTT)

		ctOut := op.GadProd(pt0NTT, ctGad, dft.FormNTT)
		ptOut := enc.Phase(ctOut)
		ptOutRef := pOp.Mul(pt0NTT, pt.WithModLen(baseLen, 0))
		pOp.InvNTTTo(ptOutRef, ptOutRef)

		diff := pOp.Sub(ptOut, ptOutRef)
		assert.True(t, checkBound(baseCvt.Convert(diff.Value), expFac+1))
	})

	t.Run("Tensor", func(t *testing.T) {
		ctVec := op.Tensor(ct0, ct1)
		skBase := (*rlwe.Element)(sk).WithModLen(baseLen, 0)

		ptOut := pOp.Mul(ctVec.Value[2], skBase)
		pOp.AddTo(ptOut, ptOut, ctVec.Value[1])
		pOp.MulTo(ptOut, ptOut, skBase)
		pOp.AddTo(ptOut, ptOut, ctVec.Value[0])
		pOp.InvNTTTo(ptOut, ptOut)

		phase0NTT := pOp.FwdNTT(enc.Phase(ct0))
		phase1NTT := pOp.FwdNTT(enc.Phase(ct1))
		ptOutRef := pOp.Mul(phase0NTT, phase1NTT)
		pOp.InvNTTTo(ptOutRef, ptOutRef)

		assert.Equal(t, ptOutRef, ptOut)

		t.Run("Relin", func(t *testing.T) {
			rlk := enc.GenRelinKey()
			ctOut := op.GadProd(ctVec.Value[2], (*rlwe.GadgetEncryption)(rlk), dft.FormNTT)
			pOp.AddTo(ctOut.Body, ctOut.Body, ctVec.Value[0])
			pOp.AddTo(ctOut.Mask, ctOut.Mask, ctVec.Value[1])
			ptOut := enc.Phase(ctOut)

			diff := pOp.Sub(ptOut, ptOutRef)
			assert.True(t, checkBound(baseCvt.Convert(diff.Value), expFac+1))
		})
	})

	t.Run("KeySwitch", func(t *testing.T) {
		skIn := rlwe.NewSecretKey(params, true, dft.FormCoeff)
		params.SecretKeyParams().Sampler().SampleTo(skIn.Value, params.FullMod())
		params.Op().FwdNTTTo(skIn.Value, skIn.Value)
		encIn := rlwe.NewEncryptorWithKey(params, skIn)
		ksk := enc.GenKeySwitchKey(skIn)

		ct := encIn.Encrypt[*rlwe.Ciphertext](pt0, dft.FormNTT)
		ctOut := op.KeySwitch(ct, ksk, dft.FormNTT)
		ptOut := enc.Phase(ctOut)

		diff := pOp.Sub(ptOut, pt0)
		assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound+expFac+1))
	})

	t.Run("Aut", func(t *testing.T) {
		var idx int
		for {
			idx = rSrc.SampleN(1, ringParams.CycloIndex())
			if op.CanAut(idx) {
				break
			}
		}
		atk := enc.GenAutKey(idx)

		ctOut := op.Aut(ct0, atk, dft.FormNTT)
		ptOut := enc.Phase(ctOut)
		ptOutRef := pOp.Aut(pt0, idx)

		diff := pOp.Sub(ptOut, ptOutRef)
		assert.True(t, checkBound(baseCvt.Convert(diff.Value), expFac*(noiseBound+expFac+1)))
	})

	t.Run("ExtProd", func(t *testing.T) {
		pt := rlwe.NewElementFrom(ts.Sample(params.Rank(), params.FullMod()), len(qAux))
		pOp.FwdNTTTo(pt, pt)
		ctRGSW := enc.Encrypt[*rlwe.RGSWEncryption](pt, dft.FormNTT)

		ctOut := op.ExtProd(ct0, ctRGSW, dft.FormNTT)
		ptOut := enc.Phase(ctOut)
		ptOutRef := pOp.Mul(pt0NTT, pt.WithModLen(baseLen, 0))
		pOp.InvNTTTo(ptOutRef, ptOutRef)

		diff := pOp.Sub(ptOut, ptOutRef)
		assert.True(t, checkBound(baseCvt.Convert(diff.Value), expFac*noiseBound+expFac+1))
	})
}

func TestCyclotomicOperator(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testOperator(t, dft.NewCyclotomicParameters(N<<1))
	})

	t.Run("type=Any", func(t *testing.T) {
		M := rSrc.SampleN(3, 1<<10)

		testOperator(t, dft.NewCyclotomicParameters(M))
	})
}

func TestAutFixedOperator(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testOperator(t, dft.NewAutFixedParameters(N<<2, N))
	})

	t.Run("type=Prime", func(t *testing.T) {
		M := num.NextPrime(rSrc.SampleN(3, 1<<10), 1)
		primes, exps := num.Factor(M - 1)
		fold := 1
		for i := range primes {
			e := rSrc.SampleN(0, exps[i])
			for j := 0; j < e; j++ {
				fold *= primes[i]
			}
		}
		N := (M - 1) / fold

		testOperator(t, dft.NewAutFixedParameters(M, N))
	})
}
