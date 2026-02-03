package rlwe_test

// import (
// 	"math"
// 	"math/big"
// 	"testing"

// 	"github.com/hienaa-org/hienaa/fhe/rlwe"
// 	"github.com/hienaa-org/hienaa/math/crt"
// 	"github.com/hienaa-org/hienaa/math/dft"
// 	"github.com/hienaa-org/hienaa/math/num"
// 	"github.com/hienaa-org/hienaa/math/vec"
// 	"github.com/stretchr/testify/assert"
// )

// func TestOperator(t *testing.T) {
// 	t.Run("type=CyclotomicPow2", func(t *testing.T) {
// 		rP := dft.NewCyclotomicParameters(1 << 11)

// 		modBits := float64(500)
// 		auxModBits := float64(120)
// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
// 		p := rlwe.ParametersLiteral{
// 			RingParams: rP,
// 			Modulus:    mod,
// 			AuxModulus: auxMod,

// 			GadgetParams: rlwe.NewRNSGadgetParameters(2),
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

// 		e := rlwe.NewEncryptor(p)
// 		o := rlwe.NewOperator(p)

// 		t.Run("Add", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			cOut := o.Add(c1, c2)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
// 			for i := 0; i < modLen; i++ {
// 				s.Value[i] = rSrc.SampleN(mod[i].Value())
// 			}

// 			pRef := eval.ScalarAdd(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarAdd(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			eval.FwdNTTTo(p2.Value, p2.Value)
// 			cOut := o.PolyAdd(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Neg", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Neg(p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Neg(c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := eval.ScalarMul(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarMul(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			pRef := eval.Mul(p1.Value, p2.Value)
// 			eval.InvNTTTo(pRef, pRef)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			cOut := o.PolyMul(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.Rank()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := p1.Copy()
// 			eval.ScalarMulAddTo(pRef.Value, p2.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.ScalarMulAddTo(c1, s, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			p3 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p3.Value, p3.Value)

// 			pRef := p1.Copy()
// 			eval.MulAddTo(pRef.Value, p3.Value, p2.Value)
// 			eval.InvNTTTo(pRef.Value, pRef.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.PolyMulAddTo(c1, p3, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.Rank()) + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("DivByAux", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)
// 			eval := o.SubEvaluatorAt(true, auxLen+modLen)
// 			evalMod := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < auxLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(auxMod[j].Value())
// 				}
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j+auxLen][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), modLen)
// 			s := crt.NewScaler(mod[:modLen], eval.Modulus())
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.DivByAux(c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Scale", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-2)) + 2)
// 			newLen := int(rSrc.SampleN(uint64(modLen-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, newLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), newLen)
// 			s := crt.NewScaler(mod[:newLen], mod[:modLen])
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Scale(c, newLen, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("GadgetProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.GadgetEncrypt(p2, true)

// 			cOut := o.GadgetProd(p1, c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank() + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Relinearisation", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c1 := e.Encrypt(p, true)
// 			c2 := e.Encrypt(p, true)
// 			rlk := e.NewRelinKey()

// 			tensor := rlwe.NewTensorCustom(rP.Rank(), modLen, 0, 3, false)

// 			eval.MulTo(tensor.Value[0], c1.Body, c2.Body)
// 			eval.MulTo(tensor.Value[1], c1.Mask, c2.Body)
// 			eval.MulAddTo(tensor.Value[1], c1.Body, c2.Mask)
// 			eval.MulTo(tensor.Value[2], c1.Mask, c2.Mask)

// 			cOut := o.Relin(tensor, rlk, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(2*rP.Rank()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("KeySwitch", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			eNew := rlwe.NewEncryptor(p)
// 			ksk := e.NewKeySwitchKey(eNew.SecretKey())

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c := eNew.Encrypt(p, true)

// 			cOut := o.KeySwitch(c, ksk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Automorphism", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			var idx int
// 			for {
// 				idx = int(rSrc.SampleN(uint64(rP.CycloOrder() - 1)))
// 				if eval.CanAut(idx) {
// 					break
// 				}
// 			}
// 			atk := e.NewAutomorphismKey(idx)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Aut(p.Value, atk.Idx)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Aut(c, atk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ExternalProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p1, true)
// 			r := e.RGSWEncrypt(p2, true)

// 			cOut := o.ExtProd(c, r, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 	})
// 	t.Run("type=CyclotomicAny", func(t *testing.T) {
// 		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))

// 		modBits := float64(500)
// 		auxModBits := float64(120)
// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
// 		p := rlwe.ParametersLiteral{
// 			RingParams: rP,
// 			Modulus:    mod,
// 			AuxModulus: auxMod,

// 			GadgetParams: rlwe.NewRNSGadgetParameters(2),
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

// 		e := rlwe.NewEncryptor(p)
// 		o := rlwe.NewOperator(p)

// 		t.Run("Add", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			cOut := o.Add(c1, c2)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
// 			for i := 0; i < modLen; i++ {
// 				s.Value[i] = rSrc.SampleN(mod[i].Value())
// 			}

// 			pRef := eval.ScalarAdd(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarAdd(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			eval.FwdNTTTo(p2.Value, p2.Value)
// 			cOut := o.PolyAdd(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Neg", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Neg(p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Neg(c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := eval.ScalarMul(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarMul(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			pRef := eval.Mul(p1.Value, p2.Value)
// 			eval.InvNTTTo(pRef, pRef)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			cOut := o.PolyMul(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := p1.Copy()
// 			eval.ScalarMulAddTo(pRef.Value, p2.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.ScalarMulAddTo(c1, s, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			p3 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p3.Value, p3.Value)

// 			pRef := p1.Copy()
// 			eval.MulAddTo(pRef.Value, p3.Value, p2.Value)
// 			eval.InvNTTTo(pRef.Value, pRef.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.PolyMulAddTo(c1, p3, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.CycloOrder()) + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("DivByAux", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)
// 			eval := o.SubEvaluatorAt(true, auxLen+modLen)
// 			evalMod := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < auxLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(auxMod[j].Value())
// 				}
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j+auxLen][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), modLen)
// 			s := crt.NewScaler(mod[:modLen], eval.Modulus())
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.DivByAux(c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Scale", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-2)) + 2)
// 			newLen := int(rSrc.SampleN(uint64(modLen-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, newLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), newLen)
// 			s := crt.NewScaler(mod[:newLen], mod[:modLen])
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Scale(c, newLen, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("GadgetProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.GadgetEncrypt(p2, true)

// 			cOut := o.GadgetProd(p1, c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank() + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Relinearisation", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c1 := e.Encrypt(p, true)
// 			c2 := e.Encrypt(p, true)
// 			rlk := e.NewRelinKey()

// 			tensor := rlwe.NewTensorCustom(rP.Rank(), modLen, 0, 3, false)

// 			eval.MulTo(tensor.Value[0], c1.Body, c2.Body)
// 			eval.MulTo(tensor.Value[1], c1.Mask, c2.Body)
// 			eval.MulAddTo(tensor.Value[1], c1.Body, c2.Mask)
// 			eval.MulTo(tensor.Value[2], c1.Mask, c2.Mask)

// 			cOut := o.Relin(tensor, rlk, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(2*rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("KeySwitch", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			eNew := rlwe.NewEncryptor(p)
// 			ksk := e.NewKeySwitchKey(eNew.SecretKey())

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c := eNew.Encrypt(p, true)

// 			cOut := o.KeySwitch(c, ksk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Automorphism", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			var idx int
// 			for {
// 				idx = int(rSrc.SampleN(uint64(rP.CycloOrder() - 1)))
// 				if eval.CanAut(idx) {
// 					break
// 				}
// 			}
// 			atk := e.NewAutomorphismKey(idx)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Aut(p.Value, atk.Idx)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Aut(c, atk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ExternalProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p1, true)
// 			r := e.RGSWEncrypt(p2, true)

// 			cOut := o.ExtProd(c, r, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 	})

// 	t.Run("type=AutFixedPow2", func(t *testing.T) {
// 		rP := dft.NewAutFixedParameters(1<<13, 1<<11)

// 		modBits := float64(500)
// 		auxModBits := float64(120)
// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
// 		p := rlwe.ParametersLiteral{
// 			RingParams: rP,
// 			Modulus:    mod,
// 			AuxModulus: auxMod,

// 			GadgetParams: rlwe.NewRNSGadgetParameters(2),
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

// 		e := rlwe.NewEncryptor(p)
// 		o := rlwe.NewOperator(p)

// 		t.Run("Add", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			cOut := o.Add(c1, c2)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
// 			for i := 0; i < modLen; i++ {
// 				s.Value[i] = rSrc.SampleN(mod[i].Value())
// 			}

// 			pRef := eval.ScalarAdd(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarAdd(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			eval.FwdNTTTo(p2.Value, p2.Value)
// 			cOut := o.PolyAdd(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Neg", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Neg(p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Neg(c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := eval.ScalarMul(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarMul(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			pRef := eval.Mul(p1.Value, p2.Value)
// 			eval.InvNTTTo(pRef, pRef)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			cOut := o.PolyMul(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := p1.Copy()
// 			eval.ScalarMulAddTo(pRef.Value, p2.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.ScalarMulAddTo(c1, s, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			p3 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p3.Value, p3.Value)

// 			pRef := p1.Copy()
// 			eval.MulAddTo(pRef.Value, p3.Value, p2.Value)
// 			eval.InvNTTTo(pRef.Value, pRef.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.PolyMulAddTo(c1, p3, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.CycloOrder()) + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("DivByAux", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)
// 			eval := o.SubEvaluatorAt(true, auxLen+modLen)
// 			evalMod := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < auxLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(auxMod[j].Value())
// 				}
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j+auxLen][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), modLen)
// 			s := crt.NewScaler(mod[:modLen], eval.Modulus())
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.DivByAux(c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Scale", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-2)) + 2)
// 			newLen := int(rSrc.SampleN(uint64(modLen-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, newLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), newLen)
// 			s := crt.NewScaler(mod[:newLen], mod[:modLen])
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Scale(c, newLen, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("GadgetProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.GadgetEncrypt(p2, true)

// 			cOut := o.GadgetProd(p1, c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.Rank() + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Relinearisation", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c1 := e.Encrypt(p, true)
// 			c2 := e.Encrypt(p, true)
// 			rlk := e.NewRelinKey()

// 			tensor := rlwe.NewTensorCustom(rP.Rank(), modLen, 0, 3, false)

// 			eval.MulTo(tensor.Value[0], c1.Body, c2.Body)
// 			eval.MulTo(tensor.Value[1], c1.Mask, c2.Body)
// 			eval.MulAddTo(tensor.Value[1], c1.Body, c2.Mask)
// 			eval.MulTo(tensor.Value[2], c1.Mask, c2.Mask)

// 			cOut := o.Relin(tensor, rlk, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(2*rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("KeySwitch", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			eNew := rlwe.NewEncryptor(p)
// 			ksk := e.NewKeySwitchKey(eNew.SecretKey())

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c := eNew.Encrypt(p, true)

// 			cOut := o.KeySwitch(c, ksk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Automorphism", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			var idx int
// 			for {
// 				idx = int(rSrc.SampleN(uint64(rP.CycloOrder() - 1)))
// 				if eval.CanAut(idx) {
// 					break
// 				}
// 			}
// 			atk := e.NewAutomorphismKey(idx)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Aut(p.Value, atk.Idx)
// 			pRefBig := eval.AsBig(pRef)
// 			c := e.Encrypt(p, true)

// 			cOut := o.Aut(c, atk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ExternalProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p1, true)
// 			r := e.RGSWEncrypt(p2, true)

// 			cOut := o.ExtProd(c, r, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
// 			assert.True(t, checkBound(res, noiseBound))
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

// 		modBits := float64(500)
// 		auxModBits := float64(120)
// 		mod, auxMod := rlwe.FindNTTPrimesFromBits(rP, modBits, auxModBits)
// 		p := rlwe.ParametersLiteral{
// 			RingParams: rP,
// 			Modulus:    mod,
// 			AuxModulus: auxMod,

// 			GadgetParams: rlwe.NewRNSGadgetParameters(2),
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

// 		e := rlwe.NewEncryptor(p)
// 		o := rlwe.NewOperator(p)

// 		t.Run("Add", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			cOut := o.Add(c1, c2)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := rlwe.NewPlainScalarCustom(modLen, 0)
// 			for i := 0; i < modLen; i++ {
// 				s.Value[i] = rSrc.SampleN(mod[i].Value())
// 			}

// 			pRef := eval.ScalarAdd(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarAdd(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Add(p1.Value, p2.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			eval.FwdNTTTo(p2.Value, p2.Value)
// 			cOut := o.PolyAdd(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Neg", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Neg(p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Neg(c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := eval.ScalarMul(p.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.ScalarMul(s, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMul", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			pRef := eval.Mul(p1.Value, p2.Value)
// 			eval.InvNTTTo(pRef, pRef)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p1, true)

// 			cOut := o.PolyMul(p2, c)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			s := &rlwe.PlainScalar{
// 				Value:  crt.NewScalar(2, eval.Modulus()),
// 				HasAux: false,
// 			}

// 			pRef := p1.Copy()
// 			eval.ScalarMulAddTo(pRef.Value, p2.Value, s.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.ScalarMulAddTo(c1, s, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("PolyMulAdd", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			p2 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 					p2.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			eval.FwdNTTTo(p1.Value, p1.Value)
// 			eval.FwdNTTTo(p2.Value, p2.Value)

// 			p3 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), eval.Modulus()),
// 				HasAux: false,
// 			}
// 			eval.FwdNTTTo(p3.Value, p3.Value)

// 			pRef := p1.Copy()
// 			eval.MulAddTo(pRef.Value, p3.Value, p2.Value)
// 			eval.InvNTTTo(pRef.Value, pRef.Value)
// 			pRefBig := eval.AsBig(pRef.Value)

// 			c1 := e.Encrypt(p1, true)
// 			c2 := e.Encrypt(p2, true)

// 			o.PolyMulAddTo(c1, p3, c2)
// 			pOut := e.Phase(c1)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.CycloOrder()) + 1))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("DivByAux", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)
// 			eval := o.SubEvaluatorAt(true, auxLen+modLen)
// 			evalMod := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, auxLen, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < auxLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(auxMod[j].Value())
// 				}
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j+auxLen][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), modLen)
// 			s := crt.NewScaler(mod[:modLen], eval.Modulus())
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.DivByAux(c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Scale", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-2)) + 2)
// 			newLen := int(rSrc.SampleN(uint64(modLen-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, newLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}

// 			pRef := crt.NewPoly(rP.Rank(), newLen)
// 			s := crt.NewScaler(mod[:newLen], mod[:modLen])
// 			s.ScaleTo(pRef, p.Value)
// 			pRefBig := eval.AsBig(pRef)

// 			c := e.Encrypt(p, true)

// 			cOut := o.Scale(c, newLen, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("GadgetProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.GadgetEncrypt(p2, true)

// 			cOut := o.GadgetProd(p1, c, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Relinearisation", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c1 := e.Encrypt(p, true)
// 			c2 := e.Encrypt(p, true)
// 			rlk := e.NewRelinKey()

// 			tensor := rlwe.NewTensorCustom(rP.Rank(), modLen, 0, 3, false)

// 			eval.MulTo(tensor.Value[0], c1.Body, c2.Body)
// 			eval.MulTo(tensor.Value[1], c1.Mask, c2.Body)
// 			eval.MulAddTo(tensor.Value[1], c1.Body, c2.Mask)
// 			eval.MulTo(tensor.Value[2], c1.Mask, c2.Mask)

// 			cOut := o.Relin(tensor, rlk, false)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64(2 * (rP.CycloOrder() + 1) * rP.CycloOrder()))

// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("KeySwitch", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			eNew := rlwe.NewEncryptor(p)
// 			ksk := e.NewKeySwitchKey(eNew.SecretKey())

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			c := eNew.Encrypt(p, true)

// 			cOut := o.KeySwitch(c, ksk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)

// 			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("Automorphism", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			eval := o.SubEvaluatorAt(false, modLen)

// 			var idx int
// 			for {
// 				idx = int(rSrc.SampleN(uint64(rP.CycloOrder() - 1)))
// 				if eval.CanAut(idx) {
// 					break
// 				}
// 			}
// 			atk := e.NewAutomorphismKey(idx)

// 			p := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			pRef := eval.Aut(p.Value, atk.Idx)
// 			pRefBig := eval.AsBig(pRef)
// 			c := e.Encrypt(p, true)

// 			cOut := o.Aut(c, atk, true)
// 			pOut := e.Phase(cOut)
// 			res := eval.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 		t.Run("ExternalProduct", func(t *testing.T) {
// 			modLen := int(rSrc.SampleN(uint64(len(mod)-1)) + 1)
// 			auxLen := len(auxMod)

// 			evalMod := o.SubEvaluatorAt(false, modLen)
// 			evalAux := o.SubEvaluatorAt(true, len(mod)+auxLen)

// 			s := crt.TernarySamplerParameters{
// 				Positive:      float64(1) / float64(3),
// 				Negative:      float64(1) / float64(3),
// 				HammingWeight: 0,
// 			}.Sampler()

// 			p1 := rlwe.NewPlainPolyCustom(rP.Rank(), modLen, 0, false)
// 			for i := 0; i < rP.Rank(); i++ {
// 				for j := 0; j < modLen; j++ {
// 					p1.Value.Coeffs[j][i] = rSrc.SampleN(mod[j].Value())
// 				}
// 			}
// 			evalMod.FwdNTTTo(p1.Value, p1.Value)

// 			p2 := &rlwe.PlainPoly{
// 				Value:  s.Sample(rP.Rank(), evalAux.Modulus()),
// 				HasAux: true,
// 			}
// 			evalAux.FwdNTTTo(p2.Value, p2.Value)

// 			p2Mod := p2.WithModIdx(vec.Range(auxLen, auxLen+modLen)...)
// 			pRef := evalMod.Mul(p1.Value, p2Mod.Value)
// 			evalMod.InvNTTTo(pRef, pRef)
// 			pRefBig := evalMod.AsBig(pRef)

// 			c := e.Encrypt(p1, true)
// 			r := e.RGSWEncrypt(p2, true)

// 			cOut := o.ExtProd(c, r, true)
// 			pOut := e.Phase(cOut)
// 			res := evalMod.AsBig(pOut.Value)
// 			for i := 0; i < rP.Rank(); i++ {
// 				res[i].Sub(res[i], pRefBig[i])
// 			}

// 			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
// 			assert.True(t, checkBound(res, noiseBound))
// 		})
// 	})
// }
