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

func checkBound(v []*big.Int, bound *big.Int) bool {
	viAbs := big.NewInt(0)
	for i := range v {
		viAbs.Abs(v[i])
		if viAbs.Cmp(bound) > 0 {
			return false
		}
	}
	return true
}

func TestSKEncryptor(t *testing.T) {
	t.Run("type=CyclotomicPow2", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 11)

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		e := rlwe.NewEncryptor(p)
		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := e.SampleRlweCustom(modLen, false, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := e.ScalarEncrypt(s, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := e.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := e.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := e.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := e.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
	})

	t.Run("type=CyclotomicAny", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		e := rlwe.NewEncryptor(p)
		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := e.SampleRlweCustom(modLen, false, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := e.ScalarEncrypt(s, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := e.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := e.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := e.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := e.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
	})

	t.Run("type=AutFixedPow2", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<12, 1<<10)

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		e := rlwe.NewEncryptor(p)
		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := e.SampleRlweCustom(modLen, false, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := e.ScalarEncrypt(s, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := e.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := e.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := e.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := e.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
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

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		e := rlwe.NewEncryptor(p)
		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := e.SampleRlweCustom(modLen, false, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := e.ScalarEncrypt(s, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Add(res[i], sBig)
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				for j := 0; j < rP.Rank(); j++ {
					res[j].Add(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := e.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := e.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				for j := 0; j < rP.Rank(); j++ {
					bodyRes[j].Add(bodyRes[j], bodyTarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := e.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := e.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
	})
}

func TestPKEncryptor(t *testing.T) {
	t.Run("type=CyclotomicPow2", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 11)

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		skEnc := rlwe.NewEncryptor(p)

		pk := skEnc.NewPublicKey()
		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.Rank()+1)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := pkEnc.SampleRlweCustom(modLen, false, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := pkEnc.ScalarEncrypt(s, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := pkEnc.Encrypt(p, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := skEnc.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := skEnc.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := pkEnc.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
	})

	t.Run("type=CyclotomicAny", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		skEnc := rlwe.NewEncryptor(p)

		pk := skEnc.NewPublicKey()
		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.Rank()+1)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := pkEnc.SampleRlweCustom(modLen, false, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := pkEnc.ScalarEncrypt(s, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := pkEnc.Encrypt(p, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := skEnc.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := skEnc.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := pkEnc.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
	})

	t.Run("type=AutFixedPow2", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<12, 1<<10)

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		skEnc := rlwe.NewEncryptor(p)

		pk := skEnc.NewPublicKey()
		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.CycloOrder()+1)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := pkEnc.SampleRlweCustom(modLen, false, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := pkEnc.ScalarEncrypt(s, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := pkEnc.Encrypt(p, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := skEnc.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := skEnc.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := pkEnc.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
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

		modBits := float64((rSrc.SampleN(4) + 2) * 50)
		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,
			AuxModulus: auxMod,

			GadgetParams: rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			},
			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
				Center: 0,
				StdDev: 3.2,
			},
		}.Compile()

		skEnc := rlwe.NewEncryptor(p)

		pk := skEnc.NewPublicKey()
		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

		o := rlwe.NewOperator(p)
		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.CycloOrder()+1)))

		t.Run("SampleRlwe", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			r := pkEnc.SampleRlweCustom(modLen, false, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			s := rlwe.NewPlainScalarCustom(modLen, 0)
			for i := 0; i < modLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			r := pkEnc.ScalarEncrypt(s, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Add(res[i], sBig)
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
			eval := o.SubEvaluatorAt(false, modLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			r := pkEnc.Encrypt(p, true)
			pOut := skEnc.Phase(r)
			res := eval.AsBig(pOut.Value)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				for j := 0; j < rP.Rank(); j++ {
					res[j].Add(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := pkEnc.GadgetEncrypt(p, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := skEnc.Phase(gEnc.Value[i])
				res := eval.AsBig(pOut.Value)
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
			for i := 0; i < modLen+auxLen; i++ {
				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
			}
			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := skEnc.SecretKey().Value.Copy()
			if sk.IsNTT {
				eval.InvNTTTo(sk, sk)
			}
			skBig := eval.AsBig(sk)

			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				for j := 0; j < rP.Rank(); j++ {
					bodyRes[j].Add(bodyRes[j], bodyTarVal)
				}

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})

		t.Run("RGSWEncrypt", func(t *testing.T) {
			modLen := len(mod)
			auxLen := len(auxMod)
			eval := o.SubEvaluatorAt(true, modLen+auxLen)

			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
				}
			}
			pBig := eval.AsBig(p.Value)

			modBig := big.NewInt(1)
			for _, modi := range eval.Modulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := skEnc.SecretKey().Value.Copy()
			pNTT := eval.FwdNTT(p.Value)
			eval.MulTo(skMsg, skMsg, pNTT)
			eval.InvNTTTo(skMsg, skMsg)
			skMsgBig := eval.AsBig(skMsg)

			gsw := pkEnc.RGSWEncrypt(p, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
				bodyRes := eval.AsBig(pBodyOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
				maskRes := eval.AsBig(pMaskOut.Value)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					maskRes[j].Sub(maskRes[j], tarVal)
				}

				assert.True(t, checkBound(bodyRes, noiseBound))
				assert.True(t, checkBound(maskRes, noiseBound))
			}
		})
	})
}
