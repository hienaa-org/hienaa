package rlwe_test

import (
	"math"
	"math/big"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

func checkBound(v []*big.Int, bound float64) bool {
	viAbs := new(big.Int)
	boundBig, _ := big.NewFloat(bound).Int(nil)
	for i := range v {
		viAbs.Abs(v[i])
		if viAbs.Cmp(boundBig) > 0 {
			return false
		}
	}
	return true
}

func expansionFactor(params dft.RingParameters) float64 {
	N, M := params.Rank(), params.CycloIndex()
	switch params.RingType() {
	case dft.TypeCyclic:
		return float64(N)
	case dft.TypeCyclotomic:
		if num.IsPowerOfTwo(M) {
			return float64(N)
		}
	case dft.TypeAutFixed:
		if num.IsPowerOfTwo(M) {
			return 2*float64(N) - 1
		}
		return 2 * (float64(M) - 1)
	}

	f := params.ModPoly()
	rem := make([]float64, N)
	sums := make([]float64, N)
	rem[0] = 1
	for k := 0; k <= 2*N-2; k++ {
		weight := float64(min(k+1, 2*N-1-k))
		for j := range sums {
			sums[j] += weight * math.Abs(rem[j])
		}
		if k == 2*N-2 {
			break
		}

		lc := rem[N-1]
		for j := N - 1; j >= 0; j-- {
			if j > 0 {
				rem[j] = rem[j-1] - lc*float64(f[j])
			} else {
				rem[j] = -lc * float64(f[j])
			}
		}
	}

	var bound float64
	for _, sum := range sums {
		bound = max(bound, sum)
	}
	return bound
}

func testEncryptor(t *testing.T, ringParams dft.RingParameters) {
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
	op := enc.PlainOp()
	us := crt.UniformSamplerParameters{}.Sampler()
	noiseBound := params.NoiseParams().Bound()

	baseLen := rSrc.SampleN(1, len(q)+1)
	baseCvt := crt.NewBigConverter(params.Op().WithModSlice(len(qAux), len(qAux)+baseLen))
	fullCvt := crt.NewBigConverter(params.Op())

	t.Run("SecretKey", func(t *testing.T) {
		t.Run("EncryptZero", func(t *testing.T) {
			ct := rlwe.NewCiphertextCustom(params.Rank(), baseLen, 0, dft.FormNTT)
			enc.EncryptTo(ct, nil, dft.FormNTT)
			ptOut := enc.Phase(ct)

			assert.True(t, checkBound(baseCvt.Convert(ptOut.Value), noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			t.Run("Scalar", func(t *testing.T) {
				pt := rlwe.NewScalarCustom(baseLen, 0)
				us.SampleTo(pt.Value, q[:baseLen])

				ct := enc.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
				ptOut := enc.Phase(ct)

				diff := op.Sub(ptOut, pt)
				assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
			})

			t.Run("Poly", func(t *testing.T) {
				pt := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
				us.SampleTo(pt.Value, q[:baseLen])

				ct := enc.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
				ptOut := enc.Phase(ct)

				diff := op.Sub(ptOut, pt)
				assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
			})
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			t.Run("Scalar", func(t *testing.T) {
				pt := rlwe.NewScalar(params, true)
				us.SampleTo(pt.Value, params.FullMod())

				ct := enc.Encrypt[*rlwe.GadgetEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					ptOut := enc.Phase(ct.Value[i])
					ptOutRef := op.Mul(pt, params.GadgetVec()[i])

					diff := op.Sub(ptOut, ptOutRef)
					assert.True(t, checkBound(fullCvt.Convert(diff.Value), noiseBound))
				}
			})

			t.Run("Poly", func(t *testing.T) {
				pt := rlwe.NewPoly(params, true, dft.FormCoeff)
				us.SampleTo(pt.Value, params.FullMod())

				ct := enc.Encrypt[*rlwe.GadgetEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					ptOut := enc.Phase(ct.Value[i])
					ptOutRef := op.Mul(pt, params.GadgetVec()[i])

					diff := op.Sub(ptOut, ptOutRef)
					assert.True(t, checkBound(fullCvt.Convert(diff.Value), noiseBound))
				}
			})
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			t.Run("Scalar", func(t *testing.T) {
				pt := rlwe.NewScalar(params, true)
				us.SampleTo(pt.Value, params.FullMod())

				ct := enc.Encrypt[*rlwe.RGSWEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					bodyOut := enc.Phase(ct.BodyStack.Value[i])
					bodyOutRef := op.Mul(pt, params.GadgetVec()[i])
					bodyDiff := op.Sub(bodyOut, bodyOutRef)

					maskOut := enc.Phase(ct.MaskStack.Value[i])
					maskOutRef := op.Mul((*rlwe.Element)(sk), params.GadgetVec()[i])
					op.MulTo(maskOutRef, maskOutRef, pt)
					op.InvNTTTo(maskOutRef, maskOutRef)
					maskDiff := op.Sub(maskOut, maskOutRef)

					assert.True(t, checkBound(fullCvt.Convert(bodyDiff.Value), noiseBound))
					assert.True(t, checkBound(fullCvt.Convert(maskDiff.Value), noiseBound))
				}
			})

			t.Run("Poly", func(t *testing.T) {
				pt := rlwe.NewPoly(params, true, dft.FormCoeff)
				us.SampleTo(pt.Value, params.FullMod())
				op.FwdNTTTo(pt, pt)

				ct := enc.Encrypt[*rlwe.RGSWEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					bodyOut := enc.Phase(ct.BodyStack.Value[i])
					bodyOutRef := op.Mul(pt, params.GadgetVec()[i])
					op.InvNTTTo(bodyOutRef, bodyOutRef)
					bodyDiff := op.Sub(bodyOut, bodyOutRef)

					maskOut := enc.Phase(ct.MaskStack.Value[i])
					maskOutRef := op.Mul((*rlwe.Element)(sk), params.GadgetVec()[i])
					op.MulTo(maskOutRef, maskOutRef, pt)
					op.InvNTTTo(maskOutRef, maskOutRef)
					maskDiff := op.Sub(maskOut, maskOutRef)

					assert.True(t, checkBound(fullCvt.Convert(bodyDiff.Value), noiseBound))
					assert.True(t, checkBound(fullCvt.Convert(maskDiff.Value), noiseBound))
				}
			})
		})
	})

	encPub := rlwe.NewEncryptorWithKey(enc.Params(), enc.GenPublicKey())
	skBound := enc.Params().SecretKeyParams().Bound()
	noiseBound = 2*expansionFactor(ringParams)*skBound*noiseBound + noiseBound

	t.Run("PublicKey", func(t *testing.T) {
		t.Run("EncryptZero", func(t *testing.T) {
			ct := rlwe.NewCiphertextCustom(params.Rank(), baseLen, 0, dft.FormNTT)
			encPub.EncryptTo(ct, nil, dft.FormNTT)
			ptOut := enc.Phase(ct)

			assert.True(t, checkBound(baseCvt.Convert(ptOut.Value), noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			t.Run("Scalar", func(t *testing.T) {
				pt := rlwe.NewScalarCustom(baseLen, 0)
				us.SampleTo(pt.Value, q[:baseLen])

				ct := encPub.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
				ptOut := enc.Phase(ct)

				diff := op.Sub(ptOut, pt)
				assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
			})

			t.Run("Poly", func(t *testing.T) {
				pt := rlwe.NewPolyCustom(params.Rank(), baseLen, 0, dft.FormCoeff)
				us.SampleTo(pt.Value, q[:baseLen])

				ct := encPub.Encrypt[*rlwe.Ciphertext](pt, dft.FormNTT)
				ptOut := enc.Phase(ct)

				diff := op.Sub(ptOut, pt)
				assert.True(t, checkBound(baseCvt.Convert(diff.Value), noiseBound))
			})
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			t.Run("Scalar", func(t *testing.T) {
				pt := rlwe.NewScalar(params, true)
				us.SampleTo(pt.Value, params.FullMod())

				ct := encPub.Encrypt[*rlwe.GadgetEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					ptOut := enc.Phase(ct.Value[i])
					ptOutRef := op.Mul(pt, params.GadgetVec()[i])

					diff := op.Sub(ptOut, ptOutRef)
					assert.True(t, checkBound(fullCvt.Convert(diff.Value), noiseBound))
				}
			})

			t.Run("Poly", func(t *testing.T) {
				pt := rlwe.NewPoly(params, true, dft.FormCoeff)
				us.SampleTo(pt.Value, params.FullMod())

				ct := encPub.Encrypt[*rlwe.GadgetEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					ptOut := enc.Phase(ct.Value[i])
					ptOutRef := op.Mul(pt, params.GadgetVec()[i])

					diff := op.Sub(ptOut, ptOutRef)
					assert.True(t, checkBound(fullCvt.Convert(diff.Value), noiseBound))
				}
			})
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			t.Run("Scalar", func(t *testing.T) {
				pt := rlwe.NewScalar(params, true)
				us.SampleTo(pt.Value, params.FullMod())

				ct := encPub.Encrypt[*rlwe.RGSWEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					bodyOut := enc.Phase(ct.BodyStack.Value[i])
					bodyOutRef := op.Mul(pt, params.GadgetVec()[i])
					bodyDiff := op.Sub(bodyOut, bodyOutRef)

					maskOut := enc.Phase(ct.MaskStack.Value[i])
					maskOutRef := op.Mul((*rlwe.Element)(sk), params.GadgetVec()[i])
					op.MulTo(maskOutRef, maskOutRef, pt)
					op.InvNTTTo(maskOutRef, maskOutRef)
					maskDiff := op.Sub(maskOut, maskOutRef)

					assert.True(t, checkBound(fullCvt.Convert(bodyDiff.Value), noiseBound))
					assert.True(t, checkBound(fullCvt.Convert(maskDiff.Value), noiseBound))
				}
			})

			t.Run("Poly", func(t *testing.T) {
				pt := rlwe.NewPoly(params, true, dft.FormCoeff)
				us.SampleTo(pt.Value, params.FullMod())
				op.FwdNTTTo(pt, pt)

				ct := encPub.Encrypt[*rlwe.RGSWEncryption](pt, dft.FormNTT)
				for i := 0; i < ct.GadgetLen(); i++ {
					bodyOut := enc.Phase(ct.BodyStack.Value[i])
					bodyOutRef := op.Mul(pt, params.GadgetVec()[i])
					op.InvNTTTo(bodyOutRef, bodyOutRef)
					bodyDiff := op.Sub(bodyOut, bodyOutRef)

					maskOut := enc.Phase(ct.MaskStack.Value[i])
					maskOutRef := op.Mul((*rlwe.Element)(sk), params.GadgetVec()[i])
					op.MulTo(maskOutRef, maskOutRef, pt)
					op.InvNTTTo(maskOutRef, maskOutRef)
					maskDiff := op.Sub(maskOut, maskOutRef)

					assert.True(t, checkBound(fullCvt.Convert(bodyDiff.Value), noiseBound))
					assert.True(t, checkBound(fullCvt.Convert(maskDiff.Value), noiseBound))
				}
			})
		})
	})
}

func TestCyclotomicEncryptor(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testEncryptor(t, dft.NewCyclotomicParameters(N<<1))
	})

	t.Run("type=Any", func(t *testing.T) {
		M := rSrc.SampleN(3, 1<<10)

		testEncryptor(t, dft.NewCyclotomicParameters(M))
	})
}

func TestAutFixedEncryptor(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testEncryptor(t, dft.NewAutFixedParameters(N<<2, N))
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

		testEncryptor(t, dft.NewAutFixedParameters(M, N))
	})
}
