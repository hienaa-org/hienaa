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

func TestOperator(t *testing.T) {
	t.Run("type=CyclotomicPow2", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 11)

		modBits := float64(500)
		auxModBits := float64(120)
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

		t.Run("Add", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pRef := pOp.Add(pt1, pt2)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			cOut := o.Add(c1, c2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			s := rlwe.NewElement(1, baseLen, 0, false)

			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}

			ptRef := pOp.Add(pt, s)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.AddPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Add(pt1, pt2)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt1, false)
			cOut := o.AddPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Neg", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Neg(pt)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Neg(c)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pOp.Mul(pt, s)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt, true)

			cOut := o.MulPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			sP := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt2, pt2)

			pRef := pOp.Mul(pt1, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt1, true)

			cOut := o.MulPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.Rank()))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt2, s)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, s)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)
			pOp.FwdNTTTo(pt2, pt2)

			pt3 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt3, pt3)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt3, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, pt3)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.Rank()) + 1))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("DivByAux", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < baseLen+auxLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(p.FullModulus()[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			sc := crt.NewScaler(baseMod[:baseLen], append(auxMod, baseMod[:baseLen]...))
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)
			cOut := o.DivByAuxModulus(c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Scale", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-2)) + 2)
			newLen := int(rSrc.SampleN(uint64(oldLen-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			for i := 0; i < oldLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			var noiseBound *big.Int
			if oldLen > newLen {
				noiseBound = big.NewInt(int64(rP.Rank()) + 1)
			} else {
				noiseBound = big.NewInt(13)
				for i := oldLen; i < newLen; i++ {
					noiseBound.Mul(noiseBound, big.NewInt(int64(baseMod[i].Value())))
				}
			}

			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(pt2, pt2)

			pt2Mod := pt2.WithModLen(baseLen, 0)
			ptRef := pOp.Mul(pt1, pt2Mod)
			pOp.InvNTTTo(ptRef, ptRef)
			pRefBig := pOp.AsBig(ptRef)

			c := e.GadgetEncrypt(pt2, true)

			cOut := o.GadgetProd(pt1, c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.Rank() + 1))
			assert.True(t, checkBound(res, noiseBound))
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
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(2*rP.Rank()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("KeySwitch", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			eNew := rlwe.NewEncryptor(p)
			ksk := e.NewKeySwitchKey(eNew.SecretKey())

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			c := eNew.Encrypt(p, true)

			cOut := o.KeySwitch(c, ksk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pRef := pOp.Aut(p, atk.Idx)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p, true)

			cOut := o.Aut(c, atk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					p1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(p1, p1)

			p2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(p2, p2)

			p2Mod := p2.WithModLen(baseLen, 0)
			pRef := pOp.Mul(p1, p2Mod)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p1, true)
			r := e.RGSWEncrypt(p2, true)

			cOut := o.ExtProd(c, r, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})
	})

	t.Run("type=CyclotomicAny", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))

		modBits := float64(500)
		auxModBits := float64(120)
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

		t.Run("Add", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pRef := pOp.Add(pt1, pt2)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			cOut := o.Add(c1, c2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			s := rlwe.NewElement(1, baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}

			ptRef := pOp.Add(pt, s)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.AddPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Add(pt1, pt2)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt1, false)
			cOut := o.AddPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Neg", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Neg(pt)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Neg(c)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pOp.Mul(pt, s)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt, true)

			cOut := o.MulPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			sP := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt2, pt2)

			pRef := pOp.Mul(pt1, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt1, true)

			cOut := o.MulPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt2, s)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, s)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)
			pOp.FwdNTTTo(pt2, pt2)

			pt3 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt3, pt3)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt3, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, pt3)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.CycloOrder()) + 1))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("DivByAux", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < baseLen+auxLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(p.FullModulus()[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			sc := crt.NewScaler(baseMod[:baseLen], append(auxMod, baseMod[:baseLen]...))
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)
			cOut := o.DivByAuxModulus(c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Scale", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-2)) + 2)
			newLen := int(rSrc.SampleN(uint64(oldLen-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			for i := 0; i < oldLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			var noiseBound *big.Int
			if oldLen > newLen {
				noiseBound = big.NewInt(int64(rP.CycloOrder()) + 1)
			} else {
				noiseBound = big.NewInt(13)
				for i := oldLen; i < newLen; i++ {
					noiseBound.Mul(noiseBound, big.NewInt(int64(baseMod[i].Value())))
				}
			}

			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(pt2, pt2)

			pt2Mod := pt2.WithModLen(baseLen, 0)
			ptRef := pOp.Mul(pt1, pt2Mod)
			pOp.InvNTTTo(ptRef, ptRef)
			pRefBig := pOp.AsBig(ptRef)

			c := e.GadgetEncrypt(pt2, true)

			cOut := o.GadgetProd(pt1, c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder() + 1))
			assert.True(t, checkBound(res, noiseBound))
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
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(2*rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("KeySwitch", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			eNew := rlwe.NewEncryptor(p)
			ksk := e.NewKeySwitchKey(eNew.SecretKey())

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			c := eNew.Encrypt(p, true)

			cOut := o.KeySwitch(c, ksk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pRef := pOp.Aut(p, atk.Idx)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p, true)

			cOut := o.Aut(c, atk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					p1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(p1, p1)

			p2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(p2, p2)

			p2Mod := p2.WithModLen(baseLen, 0)
			pRef := pOp.Mul(p1, p2Mod)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p1, true)
			r := e.RGSWEncrypt(p2, true)

			cOut := o.ExtProd(c, r, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})
	})

	t.Run("type=AutFixedPow2", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<13, 1<<11)

		modBits := float64(500)
		auxModBits := float64(120)
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

		t.Run("Add", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pRef := pOp.Add(pt1, pt2)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			cOut := o.Add(c1, c2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			s := rlwe.NewElement(1, baseLen, 0, false)

			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}

			ptRef := pOp.Add(pt, s)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.AddPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Add(pt1, pt2)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt1, false)
			cOut := o.AddPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})
		t.Run("Neg", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Neg(pt)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Neg(c)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pOp.Mul(pt, s)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt, true)

			cOut := o.MulPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})
		t.Run("PolyMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			sP := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt2, pt2)

			pRef := pOp.Mul(pt1, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt1, true)

			cOut := o.MulPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt2, s)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, s)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)
			pOp.FwdNTTTo(pt2, pt2)

			pt3 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt3, pt3)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt3, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, pt3)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.CycloOrder()) + 1))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("DivByAux", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < baseLen+auxLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(p.FullModulus()[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			sc := crt.NewScaler(baseMod[:baseLen], append(auxMod, baseMod[:baseLen]...))
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)
			cOut := o.DivByAuxModulus(c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Scale", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-2)) + 2)
			newLen := int(rSrc.SampleN(uint64(oldLen-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			for i := 0; i < oldLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			var noiseBound *big.Int
			if oldLen > newLen {
				noiseBound = big.NewInt(int64(rP.CycloOrder()) + 1)
			} else {
				noiseBound = big.NewInt(13)
				for i := oldLen; i < newLen; i++ {
					noiseBound.Mul(noiseBound, big.NewInt(int64(baseMod[i].Value())))
				}
			}

			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(pt2, pt2)

			pt2Mod := pt2.WithModLen(baseLen, 0)
			ptRef := pOp.Mul(pt1, pt2Mod)
			pOp.InvNTTTo(ptRef, ptRef)
			pRefBig := pOp.AsBig(ptRef)

			c := e.GadgetEncrypt(pt2, true)

			cOut := o.GadgetProd(pt1, c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.Rank() + 1))
			assert.True(t, checkBound(res, noiseBound))
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
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(2*rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("KeySwitch", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			eNew := rlwe.NewEncryptor(p)
			ksk := e.NewKeySwitchKey(eNew.SecretKey())

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			c := eNew.Encrypt(p, true)

			cOut := o.KeySwitch(c, ksk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(rP.Rank()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pRef := pOp.Aut(p, atk.Idx)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p, true)

			cOut := o.Aut(c, atk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					p1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(p1, p1)

			p2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(p2, p2)

			p2Mod := p2.WithModLen(baseLen, 0)
			pRef := pOp.Mul(p1, p2Mod)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p1, true)
			r := e.RGSWEncrypt(p2, true)

			cOut := o.ExtProd(c, r, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(rP.CycloOrder()) + 1)
			assert.True(t, checkBound(res, noiseBound))
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

		modBits := float64(500)
		auxModBits := float64(120)
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

		t.Run("Add", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pRef := pOp.Add(pt1, pt2)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			cOut := o.Add(c1, c2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			s := rlwe.NewElement(1, baseLen, 0, false)

			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
				s.Value.Coeffs[i][0] = rSrc.SampleN(baseMod[i].Value())
			}

			ptRef := pOp.Add(pt, s)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.AddPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Add(pt1, pt2)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt1, false)
			cOut := o.AddPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Neg", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			ptRef := pOp.Neg(pt)
			ptRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Neg(c)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], ptRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2 * 9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pOp.Mul(pt, s)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt, true)

			cOut := o.MulPlain(c, s)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(2 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("PolyMul", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			sP := crt.TernarySamplerParameters{
				Positive:      float64(1) / float64(3),
				Negative:      float64(1) / float64(3),
				HammingWeight: 0,
			}.Sampler()

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt2, pt2)

			pRef := pOp.Mul(pt1, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(pt1, true)

			cOut := o.MulPlain(c, pt2)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * int64(rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("ScalarMulAdd", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			pt1 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			pt2 := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			s := rlwe.NewElementFrom(crt.NewScalarFrom(2, baseMod[:baseLen]), 0)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt2, s)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, s)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(3 * int64(math.Round(3.2*9.2)))
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
					pt2.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)
			pOp.FwdNTTTo(pt2, pt2)

			pt3 := rlwe.NewElementFrom(sP.Sample(rP.Rank(), baseMod[:baseLen]), 0)
			pOp.FwdNTTTo(pt3, pt3)

			pRef := pt1.Copy()
			pOp.MulAddTo(pRef, pt3, pt2)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c1 := e.Encrypt(pt1, true)
			c2 := e.Encrypt(pt2, true)

			o.MulAddPlainTo(c1, c2, pt3)
			pOut := e.Phase(c1)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64(math.Round(3.2*9.2)) * (int64(rP.CycloOrder()) + 1))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("DivByAux", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)
			auxLen := len(auxMod)

			pt := rlwe.NewElement(rP.Rank(), baseLen, auxLen, false)
			for i := 0; i < baseLen+auxLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(p.FullModulus()[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			sc := crt.NewScaler(baseMod[:baseLen], append(auxMod, baseMod[:baseLen]...))
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)
			cOut := o.DivByAuxModulus(c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("Scale", func(t *testing.T) {
			oldLen := int(rSrc.SampleN(uint64(len(baseMod)-2)) + 2)
			newLen := int(rSrc.SampleN(uint64(oldLen-1)) + 1)

			pt := rlwe.NewElement(rP.Rank(), oldLen, 0, false)
			for i := 0; i < oldLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}

			ptRef := rlwe.NewElement(rP.Rank(), newLen, 0, false)
			sc := crt.NewScaler(baseMod[:newLen], baseMod[:oldLen])
			sc.ScaleTo(ptRef.Value, pt.Value)
			pRefBig := pOp.AsBig(ptRef)

			c := e.Encrypt(pt, true)

			cOut := o.Scale(c, newLen, false)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			var noiseBound *big.Int
			if oldLen > newLen {
				noiseBound = big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
			} else {
				noiseBound = big.NewInt(13)
				for i := oldLen; i < newLen; i++ {
					noiseBound.Mul(noiseBound, big.NewInt(int64(baseMod[i].Value())))
				}
			}

			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					pt1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(pt1, pt1)

			pt2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(pt2, pt2)

			pt2Mod := pt2.WithModLen(baseLen, 0)
			ptRef := pOp.Mul(pt1, pt2Mod)
			pOp.InvNTTTo(ptRef, ptRef)
			pRefBig := pOp.AsBig(ptRef)

			c := e.GadgetEncrypt(pt2, true)

			cOut := o.GadgetProd(pt1, c, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
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
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64(2 * (rP.CycloOrder() + 1) * rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
		})

		t.Run("KeySwitch", func(t *testing.T) {
			baseLen := int(rSrc.SampleN(uint64(len(baseMod)-1)) + 1)

			eNew := rlwe.NewEncryptor(p)
			ksk := e.NewKeySwitchKey(eNew.SecretKey())

			p := rlwe.NewElement(rP.Rank(), baseLen, 0, false)
			c := eNew.Encrypt(p, true)

			cOut := o.KeySwitch(c, ksk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)

			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < rP.Rank(); i++ {
				for j := 0; j < baseLen; j++ {
					p.Value.Coeffs[j][i] = rSrc.SampleN(baseMod[j].Value())
				}
			}
			pRef := pOp.Aut(p, atk.Idx)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p, true)

			cOut := o.Aut(c, atk, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
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
			for i := 0; i < baseLen; i++ {
				for j := 0; j < rP.Rank(); j++ {
					p1.Value.Coeffs[i][j] = rSrc.SampleN(baseMod[i].Value())
				}
			}
			pOp.FwdNTTTo(p1, p1)

			p2 := rlwe.NewElementFrom(s.Sample(rP.Rank(), p.FullModulus()), auxLen)
			pOp.FwdNTTTo(p2, p2)

			p2Mod := p2.WithModLen(baseLen, 0)
			pRef := pOp.Mul(p1, p2Mod)
			pOp.InvNTTTo(pRef, pRef)
			pRefBig := pOp.AsBig(pRef)

			c := e.Encrypt(p1, true)
			r := e.RGSWEncrypt(p2, true)

			cOut := o.ExtProd(c, r, true)
			pOut := e.Phase(cOut)
			res := pOp.AsBig(pOut)
			for i := 0; i < rP.Rank(); i++ {
				res[i].Sub(res[i], pRefBig[i])
			}

			noiseBound := big.NewInt(int64((rP.CycloOrder() + 1) * rP.CycloOrder()))
			assert.True(t, checkBound(res, noiseBound))
		})
	})
}
