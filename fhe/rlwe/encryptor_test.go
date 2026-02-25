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
		baseMod, auxMod := rlwe.FindNTTPrimes(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams:  rP,
			BaseModulus: baseMod,
			AuxModulus:  auxMod,

			GadgetParams: rlwe.RNSGadgetParametersLiteral{
				ChunkSize: 2,
			},
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
		pOp := o.PlainOperator()
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			r := e.SampleRLWECustom(baseLen, 0, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			s := rlwe.NewElement(1, baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				s.Value.Coeffs[0][i] = rSrc.SampleN(baseMod[i].Value())
			}
			sBig := pOp.AsBig(s)

			r := e.Encrypt(s, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig[0])).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pBig := pOp.AsBig(p)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, baseLen, auxLen, false)
			for i, modi := range p.FullModulus() {
				s.Value.Coeffs[i][0] = rSrc.SampleN(modi.Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			ptBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(pt, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(ptBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, modLen, auxLen, false)
			for i := range p.FullModulus() {
				s.Value.Coeffs[i][0] = rSrc.SampleN(p.FullModulus()[i].Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := (*rlwe.Element)(e.SecretKey()).Copy()
			if sk.IsNTT() {
				pOp.InvNTTTo(sk, sk)
			}
			skBig := pOp.AsBig(sk)

			gsw := e.RGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], new(big.Int).Mul(skBig[j], gadVecBigi[0])), modBig)
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
			modLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			ptBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := (*rlwe.Element)(e.SecretKey()).Copy()
			ptNTT := pOp.FwdNTT(pt)
			pOp.MulTo(skMsg, skMsg, ptNTT)
			pOp.InvNTTTo(skMsg, skMsg)
			skMsgBig := pOp.AsBig(skMsg)

			gsw := e.RGSWEncrypt(pt, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(ptBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi[0]), modBig)
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
		baseMod, auxMod := rlwe.FindNTTPrimes(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams:  rP,
			BaseModulus: baseMod,
			AuxModulus:  auxMod,

			GadgetParams: rlwe.RNSGadgetParametersLiteral{
				ChunkSize: 1,
			},
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
		pOp := o.PlainOperator()
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			r := e.SampleRLWECustom(baseLen, 0, false)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			s := rlwe.NewElement(1, baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}
			sBig := pOp.AsBig(s)

			r := e.Encrypt(s, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig[0])).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pBig := pOp.AsBig(p)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, baseLen, auxLen, false)
			for i, modi := range p.FullModulus() {
				s.Value.Coeffs[i][0] = rSrc.SampleN(modi.Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			modLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			ptBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(pt, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(ptBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			modLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, modLen, auxLen, false)
			for i, modi := range p.FullModulus() {
				s.Value.Coeffs[i][0] = rSrc.SampleN(modi.Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := (*rlwe.Element)(e.SecretKey()).Copy()
			if sk.IsNTT() {
				pOp.InvNTTTo(sk, sk)
			}
			skBig := pOp.AsBig(sk)

			gsw := e.RGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], new(big.Int).Mul(skBig[j], gadVecBigi[0])), modBig)
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
			modLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), modLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < modLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			ptBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := (*rlwe.Element)(e.SecretKey()).Copy()
			pNTT := pOp.FwdNTT(pt)
			pOp.MulTo(skMsg, skMsg, pNTT)
			pOp.InvNTTTo(skMsg, skMsg)
			skMsgBig := pOp.AsBig(skMsg)

			gsw := e.RGSWEncrypt(pt, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(ptBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi[0]), modBig)
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
		baseMod, auxMod := rlwe.FindNTTPrimes(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams:  rP,
			BaseModulus: baseMod,
			AuxModulus:  auxMod,

			GadgetParams: rlwe.RNSGadgetParametersLiteral{
				ChunkSize: 1,
			},
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
		pOp := o.PlainOperator()
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			r := e.SampleRLWECustom(baseLen, 0, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			s := rlwe.NewElement(1, baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}
			sBig := pOp.AsBig(s)

			r := e.Encrypt(s, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig[0])).Cmp(noiseBound) <= 0)
			assert.True(t, checkBound(res[1:], noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pBig := pOp.AsBig(p)

			r := e.Encrypt(p, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, baseLen, auxLen, false)
			for i, modi := range p.FullModulus() {
				s.Value.Coeffs[i][0] = rSrc.SampleN(modi.Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if tarVal.Cmp(halfmodBig) > 0 {
					tarVal.Sub(tarVal, modBig)
				}

				res[0].Sub(res[0], tarVal)

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("GadgetEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			ptBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(pt, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(ptBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, baseLen, auxLen, false)
			for i, modi := range p.FullModulus() {
				s.Value.Coeffs[i][0] = rSrc.SampleN(modi.Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := (*rlwe.Element)(e.SecretKey()).Copy()
			if sk.IsNTT() {
				pOp.InvNTTTo(sk, sk)
			}
			skBig := pOp.AsBig(sk)

			gsw := e.RGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], new(big.Int).Mul(skBig[j], gadVecBigi[0])), modBig)
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
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			pBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := (*rlwe.Element)(e.SecretKey()).Copy()
			pNTT := pOp.FwdNTT(pt)
			pOp.MulTo(skMsg, skMsg, pNTT)
			pOp.InvNTTTo(skMsg, skMsg)
			skMsgBig := pOp.AsBig(skMsg)

			gsw := e.RGSWEncrypt(pt, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi[0]), modBig)
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
		baseMod, auxMod := rlwe.FindNTTPrimes(rP, modBits, auxModBits)
		p := rlwe.ParametersLiteral{
			RingParams:  rP,
			BaseModulus: baseMod,
			AuxModulus:  auxMod,

			GadgetParams: rlwe.RNSGadgetParametersLiteral{
				ChunkSize: 1,
			},
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
		pOp := o.PlainOperator()
		noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))

		t.Run("SampleRlwe", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			r := e.SampleRLWECustom(baseLen, 0, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarEncrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			s := rlwe.NewElement(1, baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}
			sBig := pOp.AsBig(s)

			r := e.Encrypt(s, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Add(res[i], sBig[0])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Encrypt", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pBig := pOp.AsBig(pt)

			r := e.Encrypt(pt, true)
			pOut := e.Phase(r)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pBig[i])
			}
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, baseLen, auxLen, false)
			for i := 0; i < baseLen+auxLen; i++ {
				s.Value.Coeffs[i][0] = rSrc.SampleN(p.FullModulus()[i].Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(s, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
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
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			pBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			gEnc := e.GadgetEncrypt(pt, true)
			for i := 0; i < gEnc.GadgetLen(); i++ {
				pOut := e.Phase(gEnc.Value[i])
				res := pOp.AsBig(pOut)
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					res[j].Sub(res[j], tarVal)
				}

				assert.True(t, checkBound(res, noiseBound))
			}
		})

		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			s := rlwe.NewElement(1, baseLen, auxLen, false)
			for i := 0; i < baseLen+auxLen; i++ {
				s.Value.Coeffs[i][0] = rSrc.SampleN(p.FullModulus()[i].Value())
			}
			sBig := pOp.AsBig(s)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			sk := (*rlwe.Element)(e.SecretKey()).Copy()
			if sk.IsNTT() {
				pOp.InvNTTTo(sk, sk)
			}
			skBig := pOp.AsBig(sk)

			gsw := e.RGSWEncrypt(s, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], gadVecBigi[0]), modBig)
				if bodyTarVal.Cmp(halfmodBig) > 0 {
					bodyTarVal.Sub(bodyTarVal, modBig)
				}

				for j := 0; j < rP.Rank(); j++ {
					bodyRes[j].Add(bodyRes[j], bodyTarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig[0], new(big.Int).Mul(skBig[j], gadVecBigi[0])), modBig)
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
			baseLen := len(baseMod)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen+auxLen; j++ {
					pt.Value.Coeffs[j][i] = rSrc.SampleN(p.FullModulus()[j].Value())
				}
			}
			ptBig := pOp.AsBig(pt)

			modBig := big.NewInt(1)
			for _, modi := range p.FullModulus() {
				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
			}
			halfmodBig := new(big.Int).Rsh(modBig, 1)

			skMsg := (*rlwe.Element)(e.SecretKey()).Copy()
			pNTT := pOp.FwdNTT(pt)
			pOp.MulTo(skMsg, skMsg, pNTT)
			pOp.InvNTTTo(skMsg, skMsg)
			skMsgBig := pOp.AsBig(skMsg)

			gsw := e.RGSWEncrypt(pt, true)
			for i := 0; i < gsw.GadgetLen(); i++ {
				gadVecBigi := pOp.AsBig(o.Decomposer().GadgetVector()[i])

				pBodyOut := e.Phase(gsw.Body.Value[i])
				bodyRes := pOp.AsBig(pBodyOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(ptBig[j], gadVecBigi[0]), modBig)
					if tarVal.Cmp(halfmodBig) > 0 {
						tarVal.Sub(tarVal, modBig)
					}
					bodyRes[j].Sub(bodyRes[j], tarVal)
				}

				pMaskOut := e.Phase(gsw.Mask.Value[i])
				maskRes := pOp.AsBig(pMaskOut)
				for j := 0; j < rP.Rank(); j++ {
					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi[0]), modBig)
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
	// }

	// func TestPKEncryptor(t *testing.T) {
	// 	t.Run("type=CyclotomicPow2", func(t *testing.T) {
	// 		rP := dft.NewCyclotomicParameters(1 << 11)

	// 		modBits := float64((rSrc.SampleN(4) + 2) * 50)
	// 		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
	// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
	// 		p := rlwe.ParametersLiteral{
	// 			RingParams: rP,
	// 			Modulus:    mod,
	// 			AuxModulus: auxMod,

	// 			GadgetParams: rlwe.NewRNSGadgetParameters(1),
	// 			SecretKeyParams: crt.TernarySamplerParameters{
	// 				Positive:      float64(1) / float64(3),
	// 				Negative:      float64(1) / float64(3),
	// 				HammingWeight: 0,
	// 			},
	// 			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
	// 				Center: 0,
	// 				StdDev: 3.2,
	// 			},
	// 		}.Compile()

	// 		skEnc := rlwe.NewEncryptor(p)

	// 		pk := skEnc.NewPublicKey()
	// 		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

	// 		o := rlwe.NewOperator(p)
	// 		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.Rank()+1)))

	// 		t.Run("SampleRlwe", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			r := pkEnc.SampleRlweCustom(modLen, false, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarEncrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
	// 			for i := 0; i < modLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			r := pkEnc.ScalarEncrypt(s, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
	// 			assert.True(t, checkBound(res[1:], noiseBound))
	// 		})

	// 		t.Run("Encrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			r := pkEnc.Encrypt(p, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			for i := 0; i < rP.Rank(); i++ {
	// 				res[i].Sub(res[i], pBig[i])
	// 			}
	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if tarVal.Cmp(halfmodBig) > 0 {
	// 					tarVal.Sub(tarVal, modBig)
	// 				}

	// 				res[0].Sub(res[0], tarVal)

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("GadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.GadgetEncrypt(p, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					res[j].Sub(res[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			sk := skEnc.SecretKey().Value.Copy()
	// 			if sk.IsNTT {
	// 				eval.InvNTTTo(sk, sk)
	// 			}
	// 			skBig := eval.AsBig(sk)

	// 			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if bodyTarVal.Cmp(halfmodBig) > 0 {
	// 					bodyTarVal.Sub(bodyTarVal, modBig)
	// 				}

	// 				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})

	// 		t.Run("RGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			skMsg := skEnc.SecretKey().Value.Copy()
	// 			pNTT := eval.FwdNTT(p.Value)
	// 			eval.MulTo(skMsg, skMsg, pNTT)
	// 			eval.InvNTTTo(skMsg, skMsg)
	// 			skMsgBig := eval.AsBig(skMsg)

	// 			gsw := pkEnc.RGSWEncrypt(p, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					bodyRes[j].Sub(bodyRes[j], tarVal)
	// 				}

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})
	// 	})

	// 	t.Run("type=CyclotomicAny", func(t *testing.T) {
	// 		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))

	// 		modBits := float64((rSrc.SampleN(4) + 2) * 50)
	// 		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
	// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
	// 		p := rlwe.ParametersLiteral{
	// 			RingParams: rP,
	// 			Modulus:    mod,
	// 			AuxModulus: auxMod,

	// 			GadgetParams: rlwe.NewRNSGadgetParameters(1),
	// 			SecretKeyParams: crt.TernarySamplerParameters{
	// 				Positive:      float64(1) / float64(3),
	// 				Negative:      float64(1) / float64(3),
	// 				HammingWeight: 0,
	// 			},
	// 			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
	// 				Center: 0,
	// 				StdDev: 3.2,
	// 			},
	// 		}.Compile()

	// 		skEnc := rlwe.NewEncryptor(p)

	// 		pk := skEnc.NewPublicKey()
	// 		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

	// 		o := rlwe.NewOperator(p)
	// 		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.Rank()+1)))

	// 		t.Run("SampleRlwe", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			r := pkEnc.SampleRlweCustom(modLen, false, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarEncrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
	// 			for i := 0; i < modLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			r := pkEnc.ScalarEncrypt(s, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
	// 			assert.True(t, checkBound(res[1:], noiseBound))
	// 		})

	// 		t.Run("Encrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			r := pkEnc.Encrypt(p, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			for i := 0; i < rP.Rank(); i++ {
	// 				res[i].Sub(res[i], pBig[i])
	// 			}
	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if tarVal.Cmp(halfmodBig) > 0 {
	// 					tarVal.Sub(tarVal, modBig)
	// 				}

	// 				res[0].Sub(res[0], tarVal)

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("GadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.GadgetEncrypt(p, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					res[j].Sub(res[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			sk := skEnc.SecretKey().Value.Copy()
	// 			if sk.IsNTT {
	// 				eval.InvNTTTo(sk, sk)
	// 			}
	// 			skBig := eval.AsBig(sk)

	// 			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if bodyTarVal.Cmp(halfmodBig) > 0 {
	// 					bodyTarVal.Sub(bodyTarVal, modBig)
	// 				}

	// 				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})

	// 		t.Run("RGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			skMsg := skEnc.SecretKey().Value.Copy()
	// 			pNTT := eval.FwdNTT(p.Value)
	// 			eval.MulTo(skMsg, skMsg, pNTT)
	// 			eval.InvNTTTo(skMsg, skMsg)
	// 			skMsgBig := eval.AsBig(skMsg)

	// 			gsw := pkEnc.RGSWEncrypt(p, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					bodyRes[j].Sub(bodyRes[j], tarVal)
	// 				}

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})
	// 	})

	// 	t.Run("type=AutFixedPow2", func(t *testing.T) {
	// 		rP := dft.NewAutFixedParameters(1<<12, 1<<10)

	// 		modBits := float64((rSrc.SampleN(4) + 2) * 50)
	// 		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
	// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
	// 		p := rlwe.ParametersLiteral{
	// 			RingParams: rP,
	// 			Modulus:    mod,
	// 			AuxModulus: auxMod,

	// 			GadgetParams: rlwe.NewRNSGadgetParameters(1),
	// 			SecretKeyParams: crt.TernarySamplerParameters{
	// 				Positive:      float64(1) / float64(3),
	// 				Negative:      float64(1) / float64(3),
	// 				HammingWeight: 0,
	// 			},
	// 			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
	// 				Center: 0,
	// 				StdDev: 3.2,
	// 			},
	// 		}.Compile()

	// 		skEnc := rlwe.NewEncryptor(p)

	// 		pk := skEnc.NewPublicKey()
	// 		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

	// 		o := rlwe.NewOperator(p)
	// 		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.CycloOrder()+1)))

	// 		t.Run("SampleRlwe", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			r := pkEnc.SampleRlweCustom(modLen, false, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarEncrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
	// 			for i := 0; i < modLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			r := pkEnc.ScalarEncrypt(s, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, new(big.Int).Abs(new(big.Int).Sub(res[0], sBig)).Cmp(noiseBound) <= 0)
	// 			assert.True(t, checkBound(res[1:], noiseBound))
	// 		})

	// 		t.Run("Encrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			r := pkEnc.Encrypt(p, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			for i := 0; i < rP.Rank(); i++ {
	// 				res[i].Sub(res[i], pBig[i])
	// 			}
	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if tarVal.Cmp(halfmodBig) > 0 {
	// 					tarVal.Sub(tarVal, modBig)
	// 				}

	// 				res[0].Sub(res[0], tarVal)

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("GadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.GadgetEncrypt(p, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					res[j].Sub(res[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			sk := skEnc.SecretKey().Value.Copy()
	// 			if sk.IsNTT {
	// 				eval.InvNTTTo(sk, sk)
	// 			}
	// 			skBig := eval.AsBig(sk)

	// 			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if bodyTarVal.Cmp(halfmodBig) > 0 {
	// 					bodyTarVal.Sub(bodyTarVal, modBig)
	// 				}

	// 				bodyRes[0].Sub(bodyRes[0], bodyTarVal)

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})

	// 		t.Run("RGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			skMsg := skEnc.SecretKey().Value.Copy()
	// 			pNTT := eval.FwdNTT(p.Value)
	// 			eval.MulTo(skMsg, skMsg, pNTT)
	// 			eval.InvNTTTo(skMsg, skMsg)
	// 			skMsgBig := eval.AsBig(skMsg)

	// 			gsw := pkEnc.RGSWEncrypt(p, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					bodyRes[j].Sub(bodyRes[j], tarVal)
	// 				}

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})
	// 	})

	// 	t.Run("type=AutFixedPrime", func(t *testing.T) {
	// 		M := num.MustNextPrime(int(rSrc.SampleN(1<<12)), 1)
	// 		primes, exps := num.Factor(M - 1)
	// 		fold := 1
	// 		for i := range primes {
	// 			e := int(rSrc.SampleN(uint64(exps[i])))
	// 			for j := 0; j < e; j++ {
	// 				fold *= primes[i]
	// 			}
	// 		}
	// 		N := (M - 1) / fold
	// 		rP := dft.NewAutFixedParameters(M, N)

	// 		modBits := float64((rSrc.SampleN(4) + 2) * 50)
	// 		auxModBits := float64((rSrc.SampleN(2) + 1) * 50)
	// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
	// 		p := rlwe.ParametersLiteral{
	// 			RingParams: rP,
	// 			Modulus:    mod,
	// 			AuxModulus: auxMod,

	// 			GadgetParams: rlwe.NewRNSGadgetParameters(1),
	// 			SecretKeyParams: crt.TernarySamplerParameters{
	// 				Positive:      float64(1) / float64(3),
	// 				Negative:      float64(1) / float64(3),
	// 				HammingWeight: 0,
	// 			},
	// 			NoiseParams: crt.RoundedGaussianSamplerParameters[float64]{
	// 				Center: 0,
	// 				StdDev: 3.2,
	// 			},
	// 		}.Compile()

	// 		skEnc := rlwe.NewEncryptor(p)

	// 		pk := skEnc.NewPublicKey()
	// 		pkEnc := rlwe.NewEncryptorWithPK(p, pk)

	// 		o := rlwe.NewOperator(p)
	// 		noiseBound := big.NewInt(int64(math.Round(3.2*9.2) * float64(rP.CycloOrder()+1)))

	// 		t.Run("SampleRlwe", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			r := pkEnc.SampleRlweCustom(modLen, false, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarEncrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
	// 			for i := 0; i < modLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			r := pkEnc.ScalarEncrypt(s, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			for i := 0; i < rP.Rank(); i++ {
	// 				res[i].Add(res[i], sBig)
	// 			}
	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("Encrypt", func(t *testing.T) {
	// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
	// 			eval := o.SubEvaluatorAt(false, modLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			r := pkEnc.Encrypt(p, true)
	// 			pOut := skEnc.Phase(r)
	// 			res := eval.AsBig(pOut.Value)

	// 			for i := 0; i < rP.Rank(); i++ {
	// 				res[i].Sub(res[i], pBig[i])
	// 			}
	// 			assert.True(t, checkBound(res, noiseBound))
	// 		})

	// 		t.Run("ScalarGadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.ScalarGadgetEncrypt(s, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if tarVal.Cmp(halfmodBig) > 0 {
	// 					tarVal.Sub(tarVal, modBig)
	// 				}

	// 				for j := 0; j < rP.Rank(); j++ {
	// 					res[j].Add(res[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("GadgetEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			gEnc := pkEnc.GadgetEncrypt(p, true)
	// 			for i := 0; i < gEnc.GadgetLen(); i++ {
	// 				pOut := skEnc.Phase(gEnc.Value[i])
	// 				res := eval.AsBig(pOut.Value)
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					res[j].Sub(res[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(res, noiseBound))
	// 			}
	// 		})

	// 		t.Run("ScalarRGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			s := rlwe.NewPlainScalarCustom(modLen, auxLen)
	// 			for i := 0; i < modLen+auxLen; i++ {
	// 				s.Value[i] = rSrc.SampleN(eval.Modulus()[i].Value())
	// 			}
	// 			sBig := crt.AsBigScalar(s.Value, eval.Modulus())

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			sk := skEnc.SecretKey().Value.Copy()
	// 			if sk.IsNTT {
	// 				eval.InvNTTTo(sk, sk)
	// 			}
	// 			skBig := eval.AsBig(sk)

	// 			gsw := pkEnc.ScalarRGSWEncrypt(s, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				bodyTarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, gadVecBigi), modBig)
	// 				if bodyTarVal.Cmp(halfmodBig) > 0 {
	// 					bodyTarVal.Sub(bodyTarVal, modBig)
	// 				}

	// 				for j := 0; j < rP.Rank(); j++ {
	// 					bodyRes[j].Add(bodyRes[j], bodyTarVal)
	// 				}

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(sBig, new(big.Int).Mul(skBig[j], gadVecBigi)), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	// 				assert.True(t, checkBound(bodyRes, noiseBound))
	// 				assert.True(t, checkBound(maskRes, noiseBound))
	// 			}
	// 		})

	// 		t.Run("RGSWEncrypt", func(t *testing.T) {
	// 			modLen := len(mod)
	// 			auxLen := len(auxMod)
	// 			eval := o.SubEvaluatorAt(true, modLen+auxLen)

	// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
	// 			for i := 0; i < rP.Rank(); i++ {
	// 				for j := 0; j < modLen+auxLen; j++ {
	// 					p.Value.Coeffs[j][i] = rSrc.SampleN(eval.Modulus()[j].Value())
	// 				}
	// 			}
	// 			pBig := eval.AsBig(p.Value)

	// 			modBig := big.NewInt(1)
	// 			for _, modi := range eval.Modulus() {
	// 				modBig.Mul(modBig, big.NewInt(int64(modi.Value())))
	// 			}
	// 			halfmodBig := new(big.Int).Rsh(modBig, 1)

	// 			skMsg := skEnc.SecretKey().Value.Copy()
	// 			pNTT := eval.FwdNTT(p.Value)
	// 			eval.MulTo(skMsg, skMsg, pNTT)
	// 			eval.InvNTTTo(skMsg, skMsg)
	// 			skMsgBig := eval.AsBig(skMsg)

	// 			gsw := pkEnc.RGSWEncrypt(p, true)
	// 			for i := 0; i < gsw.GadgetLen(); i++ {
	// 				gadVecBigi := crt.AsBigScalar(o.Decmp.GadgetVector()[i], eval.Modulus())

	// 				pBodyOut := skEnc.Phase(gsw.Body.Value[i])
	// 				bodyRes := eval.AsBig(pBodyOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(pBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					bodyRes[j].Sub(bodyRes[j], tarVal)
	// 				}

	// 				pMaskOut := skEnc.Phase(gsw.Mask.Value[i])
	// 				maskRes := eval.AsBig(pMaskOut.Value)
	// 				for j := 0; j < rP.Rank(); j++ {
	// 					tarVal := new(big.Int).Mod(new(big.Int).Mul(skMsgBig[j], gadVecBigi), modBig)
	// 					if tarVal.Cmp(halfmodBig) > 0 {
	// 						tarVal.Sub(tarVal, modBig)
	// 					}
	// 					maskRes[j].Sub(maskRes[j], tarVal)
	// 				}

	//				assert.True(t, checkBound(bodyRes, noiseBound))
	//				assert.True(t, checkBound(maskRes, noiseBound))
	//			}
	//		})
	//	})
}
