package rlwe_test

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/stretchr/testify/assert"
)

func TestPlainScalarOperator(t *testing.T) {
	t.Run("type=Scalar", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(1 << 4)
		mod := dft.MustFindNextNTTPrimes(rP, 60, 4)
		pL := rlwe.ParametersLiteral{
			RingParams: rP,
			Modulus:    mod,

			GadgetParams:    rlwe.NewRNSGadgetParameters(1),
			SecretKeyParams: crt.TernarySamplerParameters{},
			NoiseParams:     crt.TernarySamplerParameters{},
		}
		p := pL.Compile()
		o := rlwe.NewPlainOperator(p)

		bigMod := big.NewInt(1)
		for i := range mod {
			bigMod.Mul(bigMod, new(big.Int).SetUint64(mod[i].Value()))
		}

		rnd := rand.New(rand.NewSource(0))

		t.Run("AddScalar", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1))) + 1

			v0 := new(big.Int).Rand(rnd, bigMod)
			v1 := new(big.Int).Rand(rnd, bigMod)
			vOutRef := new(big.Int).Add(v0, v1)
			vOutRef.Mod(vOutRef, bigMod)
			sOutRef := crt.NewScalar(vOutRef, mod)[:modLen]

			s0 := &rlwe.PlainScalar{Value: crt.NewScalar(v0, mod[:modLen])}
			s1 := &rlwe.PlainScalar{Value: crt.NewScalar(v1, mod[:modLen])}
			sOut := o.AddScalar(s0, s1)

			assert.Equal(t, sOut.Value, sOutRef)
		})
		t.Run("SubScalar", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1))) + 1

			v0 := new(big.Int).Rand(rnd, bigMod)
			v1 := new(big.Int).Rand(rnd, bigMod)
			vOutRef := new(big.Int).Sub(v0, v1)
			vOutRef.Mod(vOutRef, bigMod)
			sOutRef := crt.NewScalar(vOutRef, mod)[:modLen]

			s0 := &rlwe.PlainScalar{Value: crt.NewScalar(v0, mod[:modLen])}
			s1 := &rlwe.PlainScalar{Value: crt.NewScalar(v1, mod[:modLen])}
			sOut := o.SubScalar(s0, s1)

			assert.Equal(t, sOut.Value, sOutRef)
		})
		t.Run("NegScalar", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1))) + 1

			v := new(big.Int).Rand(rnd, bigMod)
			vOutRef := new(big.Int).Neg(v)
			vOutRef.Mod(vOutRef, bigMod)
			sOutRef := crt.NewScalar(vOutRef, mod)[:modLen]

			s := &rlwe.PlainScalar{Value: crt.NewScalar(v, mod[:modLen])}
			sOut := o.NegScalar(s)

			assert.Equal(t, sOut.Value, sOutRef)
		})
		t.Run("MulScalar", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1))) + 1

			v0 := new(big.Int).Rand(rnd, bigMod)
			v1 := new(big.Int).Rand(rnd, bigMod)
			vOutRef := new(big.Int).Mul(v0, v1)
			vOutRef.Mod(vOutRef, bigMod)
			sOutRef := crt.NewScalar(vOutRef, mod)[:modLen]

			s0 := &rlwe.PlainScalar{Value: crt.NewScalar(v0, mod[:modLen])}
			s1 := &rlwe.PlainScalar{Value: crt.NewScalar(v1, mod[:modLen])}
			sOut := o.MulScalar(s0, s1)

			assert.Equal(t, sOut.Value, sOutRef)
		})
		t.Run("MulAddScalarTo", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1))) + 1

			v0 := new(big.Int).Rand(rnd, bigMod)
			v1 := new(big.Int).Rand(rnd, bigMod)
			v2 := new(big.Int).Rand(rnd, bigMod)
			vOutRef := new(big.Int).Add(new(big.Int).Mul(v0, v1), v2)
			vOutRef.Mod(vOutRef, bigMod)
			sOutRef := crt.NewScalar(vOutRef, mod)[:modLen]

			s0 := &rlwe.PlainScalar{Value: crt.NewScalar(v0, mod[:modLen])}
			s1 := &rlwe.PlainScalar{Value: crt.NewScalar(v1, mod[:modLen])}
			sOut := &rlwe.PlainScalar{Value: crt.NewScalar(v2, mod[:modLen])}
			o.MulAddScalarTo(sOut, s0, s1)

			assert.Equal(t, sOut.Value, sOutRef)
		})
		t.Run("MulSubScalarTo", func(t *testing.T) {
			modLen := int(rSrc.SampleN(uint64(len(mod)-1))) + 1

			v0 := new(big.Int).Rand(rnd, bigMod)
			v1 := new(big.Int).Rand(rnd, bigMod)
			v2 := new(big.Int).Rand(rnd, bigMod)
			vOutRef := new(big.Int).Sub(v2, new(big.Int).Mul(v0, v1))
			vOutRef.Mod(vOutRef, bigMod)
			sOutRef := crt.NewScalar(vOutRef, mod)[:modLen]

			s0 := &rlwe.PlainScalar{Value: crt.NewScalar(v0, mod[:modLen])}
			s1 := &rlwe.PlainScalar{Value: crt.NewScalar(v1, mod[:modLen])}
			sOut := &rlwe.PlainScalar{Value: crt.NewScalar(v2, mod[:modLen])}
			o.MulSubScalarTo(sOut, s0, s1)

			assert.Equal(t, sOut.Value, sOutRef)
		})
	})
}

// func TestPlainPolyOperator(t *testing.T) {
// 	t.Run("type=CyclotomicPow2", func(t *testing.T) {
// 		rP := dft.NewCyclotomicParameters(1 << 11)
// 		mod := dft.MustFindNextNTTPrimes(rP, 60, 4)
// 		pL := rlwe.ParametersLiteral{
// 			RingParams: rP,
// 			Modulus:    mod,

// 			GadgetParams:    rlwe.NewRNSGadgetParameters(1),
// 			SecretKeyParams: crt.TernarySamplerParameters{},
// 			NoiseParams:     crt.TernarySamplerParameters{},
// 		}
// 		p := pL.Compile()
// 		o := rlwe.NewPlainOperator(p)
// 		_ = o

// 		t.Run("FwdNTT", func(t *testing.T) {

// 		})
// 		t.Run("InvNTT", func(t *testing.T) {

// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarSub", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulSub", func(t *testing.T) {

// 		})
// 		t.Run("Add", func(t *testing.T) {

// 		})
// 		t.Run("Sub", func(t *testing.T) {

// 		})
// 		t.Run("Mul", func(t *testing.T) {

// 		})
// 		t.Run("MulAddTo", func(t *testing.T) {

// 		})
// 		t.Run("MulSubTo", func(t *testing.T) {

// 		})
// 	})

// 	t.Run("type=CyclotomicAny", func(t *testing.T) {
// 		t.Run("FwdNTT", func(t *testing.T) {

// 		})
// 		t.Run("InvNTT", func(t *testing.T) {

// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarSub", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulSub", func(t *testing.T) {

// 		})
// 		t.Run("Add", func(t *testing.T) {

// 		})
// 		t.Run("Sub", func(t *testing.T) {

// 		})
// 		t.Run("Mul", func(t *testing.T) {

// 		})
// 		t.Run("MulAddTo", func(t *testing.T) {

// 		})
// 		t.Run("MulSubTo", func(t *testing.T) {

// 		})
// 	})

// 	t.Run("type=AutFixedPow2", func(t *testing.T) {
// 		t.Run("FwdNTT", func(t *testing.T) {

// 		})
// 		t.Run("InvNTT", func(t *testing.T) {

// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarSub", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulSub", func(t *testing.T) {

// 		})
// 		t.Run("Add", func(t *testing.T) {

// 		})
// 		t.Run("Sub", func(t *testing.T) {

// 		})
// 		t.Run("Mul", func(t *testing.T) {

// 		})
// 		t.Run("MulAddTo", func(t *testing.T) {

// 		})
// 		t.Run("MulSubTo", func(t *testing.T) {

// 		})
// 	})

// 	t.Run("type=AutFixedPrime", func(t *testing.T) {
// 		t.Run("FwdNTT", func(t *testing.T) {

// 		})
// 		t.Run("InvNTT", func(t *testing.T) {

// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarSub", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulSub", func(t *testing.T) {

// 		})
// 		t.Run("Add", func(t *testing.T) {

// 		})
// 		t.Run("Sub", func(t *testing.T) {

// 		})
// 		t.Run("Mul", func(t *testing.T) {

// 		})
// 		t.Run("MulAddTo", func(t *testing.T) {

// 		})
// 		t.Run("MulSubTo", func(t *testing.T) {

// 		})
// 	})
// }

// func TestOperator(t *testing.T) {
// 	t.Run("type=Ciphertext", func(t *testing.T) {
// 		t.Run("FwdNTT", func(t *testing.T) {

// 		})
// 		t.Run("InvNTT", func(t *testing.T) {

// 		})
// 		t.Run("Add", func(t *testing.T) {

// 		})
// 		t.Run("ScalarAdd", func(t *testing.T) {

// 		})
// 		t.Run("PolyAdd", func(t *testing.T) {

// 		})
// 		t.Run("Sub", func(t *testing.T) {

// 		})
// 		t.Run("ScalarSub", func(t *testing.T) {

// 		})
// 		t.Run("SubPoly", func(t *testing.T) {

// 		})
// 		t.Run("Neg", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMul", func(t *testing.T) {

// 		})
// 		t.Run("PolyMul", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulAdd", func(t *testing.T) {

// 		})
// 		t.Run("PolyMulAdd", func(t *testing.T) {

// 		})
// 		t.Run("ScalarMulSub", func(t *testing.T) {

// 		})
// 		t.Run("PolyMulSub", func(t *testing.T) {

// 		})
// 		t.Run("DivByAux", func(t *testing.T) {

// 		})
// 		t.Run("Scale", func(t *testing.T) {

// 		})
// 		t.Run("GadgetProduct", func(t *testing.T) {

// 		})
// 		t.Run("Relinearisation", func(t *testing.T) {

// 		})
// 		t.Run("KeySwitch", func(t *testing.T) {

// 		})
// 		t.Run("Automorphism", func(t *testing.T) {

// 		})
// 		t.Run("ExternalProduct", func(t *testing.T) {

// 		})
// 	})
// }
