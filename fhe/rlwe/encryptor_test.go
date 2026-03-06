package rlwe_test

import (
	"math/big"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
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

func randomElementTo(eOut *rlwe.Element, mod []*num.Modulus) {
	for i := range mod {
		for j := range eOut.Value.Coeffs[i] {
			eOut.Value.Coeffs[i][j] = rSrc.SampleN(mod[i].Value())
		}
	}
	eOut.Value.IsNTT = false
}

func testSKEncryptor(t *testing.T, rP dft.RingParameters) {
	baseBits := float64((rSrc.SampleN(4) + 2) * 50)
	auxBits := float64((rSrc.SampleN(2) + 1) * 50)
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
	noiseBound := noiseParams.Bound()

	t.Run("SampleRlwe", func(t *testing.T) {
		baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

		r := e.SampleRLWECustom(baseLen, 0, true)
		pOut := e.Phase(r)
		res := pOp.AsBig(pOut)

		assert.True(t, checkBound(res, noiseBound))
	})

	t.Run("Encrypt", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(1, baseLen, 0, false)
			randomElementTo(pt, baseMod[:baseLen])

			r := e.Encrypt(pt, true)
			pOut := e.Phase(r)

			diff := pOp.Sub(pOut, pt)
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			randomElementTo(pt, baseMod[:baseLen])

			r := e.Encrypt(pt, true)
			pOut := e.Phase(r)

			diff := pOp.Sub(pOut, pt)
			assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
		})
	})

	t.Run("GadgetEncrypt", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(1, baseLen, auxLen, false)
			randomElementTo(pt, p.FullModulus())

			gEnc := e.GadgetEncrypt(pt, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				tar := pOp.Mul(pt, o.Decomposer().GadgetVector()[i])
				diff := pOp.Sub(pOut, tar)
				assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
			}
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			randomElementTo(pt, p.FullModulus())

			gEnc := e.GadgetEncrypt(pt, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				tar := pOp.Mul(pt, o.Decomposer().GadgetVector()[i])

				diff := pOp.Sub(pOut, tar)
				assert.True(t, checkBound(pOp.AsBig(diff), noiseBound))
			}
		})
	})

	t.Run("RGSWEncrypt", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(1, baseLen, auxLen, false)
			randomElementTo(pt, p.FullModulus())

			gsw := e.RGSWEncrypt(pt, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyTar := pOp.Mul(pt, o.Decomposer().GadgetVector()[i])
				bodyDiff := pOp.Sub(pBodyOut, bodyTar)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskTar := pOp.Mul((*rlwe.Element)(e.SecretKey()), o.Decomposer().GadgetVector()[i])
				pOp.MulTo(maskTar, maskTar, pt)
				pOp.InvNTTTo(maskTar, maskTar)
				maskDiff := pOp.Sub(pMaskOut, maskTar)

				assert.True(t, checkBound(pOp.AsBig(bodyDiff), noiseBound))
				assert.True(t, checkBound(pOp.AsBig(maskDiff), noiseBound))
			}
		})

		t.Run("Poly", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			randomElementTo(pt, p.FullModulus())
			pOp.FwdNTTTo(pt, pt)

			gsw := e.RGSWEncrypt(pt, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyTar := pOp.Mul(pt, o.Decomposer().GadgetVector()[i])
				pOp.InvNTTTo(bodyTar, bodyTar)
				bodyDiff := pOp.Sub(pBodyOut, bodyTar)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskTar := pOp.Mul((*rlwe.Element)(e.SecretKey()), o.Decomposer().GadgetVector()[i])
				pOp.MulTo(maskTar, maskTar, pt)
				pOp.InvNTTTo(maskTar, maskTar)
				maskDiff := pOp.Sub(pMaskOut, maskTar)

				assert.True(t, checkBound(pOp.AsBig(bodyDiff), noiseBound))
				assert.True(t, checkBound(pOp.AsBig(maskDiff), noiseBound))
			}
		})
	})
}

func TestSKEncryptor(t *testing.T) {
	t.Run("type=CyclotomicPow2", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 11)
		testSKEncryptor(t, rP)
	})

	t.Run("type=CyclotomicAny", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))
		testSKEncryptor(t, rP)
	})

	t.Run("type=AutFixedPow2", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<12, 1<<10)
		testSKEncryptor(t, rP)
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
		testSKEncryptor(t, rP)
	})
}
