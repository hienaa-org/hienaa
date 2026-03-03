package bfv_test

import (
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/bfv"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc     = csprng.NewUniformSamplerWithSeed(nil)
	skParams = crt.TernarySamplerParameters{
		Positive: 1.0 / 3.0,
		Negative: 1.0 / 3.0,
	}
	noiseParams = crt.RoundedGaussianSamplerParameters[float64]{
		Center: 0,
		StdDev: 3.2,
	}
	baseModBits = float64(200)
	auxModBits  = float64(100)
)

func randMsg(packLen int, msgMod *num.Modulus) []uint64 {
	msg := make([]uint64, packLen)
	for i := range msg {
		msg[i] = uint64(rand.Intn(int(msgMod.Value())))
	}
	return msg
}

func testOperator(t *testing.T, rP dft.RingParameters, q *num.Modulus) {
	baseMod, auxMod := rlwe.FindNTTPrimes(rP, baseModBits, auxModBits)
	p := bfv.ParametersLiteral{
		RLWEParams: rlwe.ParametersLiteral{
			RingParams:      rP,
			BaseModulus:     baseMod,
			AuxModulus:      auxMod,
			SecretKeyParams: skParams,
			NoiseParams:     noiseParams,
		},
		MessageModulus: q,
		IsLevelled:     true,
		EstimType:      bfv.TypeVariance,
	}.Compile()

	o := bfv.NewOperator(p)
	e := bfv.NewEncryptor(p)
	pack := o.Packer()

	t.Run("Encrypt", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg := uint64(rand.Intn(int(q.Value())))

			ct := e.Encrypt(bfv.NewScalarFrom(msg, q), true)
			res := e.Decrypt(ct)

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloOrder()) {
				for i := range res.Coeffs[0] {
					assert.Equal(t, msg, num.Neg(res.Coeffs[0][i], q))
				}
			} else {
				assert.Equal(t, msg, res.Coeffs[0][0])
				for i := 1; i < len(res.Coeffs[0]); i++ {
					assert.Equal(t, uint64(0), res.Coeffs[0][i])
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg := randMsg(pack.PackLen(), q)

			ct := e.Encrypt(pack.Pack(msg), true)
			res := e.Decrypt(ct)
			out := pack.UnPack(res)

			assert.Equal(t, msg, out)
		})
	})

	t.Run("Add", func(t *testing.T) {
		msg0 := randMsg(pack.PackLen(), q)
		msg1 := randMsg(pack.PackLen(), q)
		msgRef := vec.Add(msg0, msg1, q)

		ct0 := e.Encrypt(pack.Pack(msg0), true)
		ct1 := e.Encrypt(pack.Pack(msg1), true)

		ctOut := o.Add(ct0, ct1, true)
		res := e.Decrypt(ctOut)
		out := pack.UnPack(res)

		assert.Equal(t, msgRef, out)
	})

	t.Run("AddPlain", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := uint64(rand.Intn(int(q.Value())))
			msg1 := uint64(rand.Intn(int(q.Value())))
			msgRef := num.Add(msg0, msg1, q)

			ct := e.Encrypt(bfv.NewScalarFrom(msg0, q), true)
			pt := bfv.NewScalarFrom(msg1, q)

			ctOut := o.AddPlain(ct, pt, true)
			res := e.Decrypt(ctOut)

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloOrder()) {
				for i := range res.Coeffs[0] {
					assert.Equal(t, msgRef, num.Neg(res.Coeffs[0][i], q))
				}
			} else {
				assert.Equal(t, msgRef, res.Coeffs[0][0])
				for i := 1; i < len(res.Coeffs[0]); i++ {
					assert.Equal(t, uint64(0), res.Coeffs[0][i])
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen(), q)
			msg1 := randMsg(pack.PackLen(), q)
			msgRef := vec.Add(msg0, msg1, q)

			ct := e.Encrypt(pack.Pack(msg0), false)
			pt := pack.Pack(msg1)

			ctOut := o.AddPlain(ct, pt, true)
			res := e.Decrypt(ctOut)
			out := pack.UnPack(res)

			assert.Equal(t, msgRef, out)
		})
	})

	t.Run("AddElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := uint64(rand.Intn(int(q.Value())))
			msg1 := uint64(rand.Intn(int(q.Value())))
			msgRef := num.Add(msg0, msg1, q)

			el := o.Encoder().Encode(bfv.NewScalarFrom(msg1, q), false, true)
			ct := e.Encrypt(bfv.NewScalarFrom(msg0, q), true)

			ctOut := o.AddElement(ct, el, true)
			res := e.Decrypt(ctOut)

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloOrder()) {
				for i := range res.Coeffs[0] {
					assert.Equal(t, msgRef, num.Neg(res.Coeffs[0][i], q))
				}
			} else {
				assert.Equal(t, msgRef, res.Coeffs[0][0])
				for i := 1; i < len(res.Coeffs[0]); i++ {
					assert.Equal(t, uint64(0), res.Coeffs[0][i])
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen(), q)
			msg1 := randMsg(pack.PackLen(), q)
			msgRef := vec.Add(msg0, msg1, q)

			ct := e.Encrypt(pack.Pack(msg0), false)
			el := o.Encoder().Encode(pack.Pack(msg1), false, true)

			ctOut := o.AddElement(ct, el, true)
			res := e.Decrypt(ctOut)
			out := pack.UnPack(res)

			assert.Equal(t, msgRef, out)
		})
	})

	t.Run("Sub", func(t *testing.T) {
		msg0 := randMsg(pack.PackLen(), q)
		msg1 := randMsg(pack.PackLen(), q)
		msgRef := vec.Sub(msg0, msg1, q)

		ct0 := e.Encrypt(pack.Pack(msg0), true)
		ct1 := e.Encrypt(pack.Pack(msg1), true)

		ctOut := o.Sub(ct0, ct1, true)
		res := e.Decrypt(ctOut)
		out := pack.UnPack(res)

		assert.Equal(t, msgRef, out)
	})

	t.Run("SubPlain", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := uint64(rand.Intn(int(q.Value())))
			msg1 := uint64(rand.Intn(int(q.Value())))
			msgRef := num.Sub(msg0, msg1, q)

			ct := e.Encrypt(bfv.NewScalarFrom(msg0, q), true)
			pt := bfv.NewScalarFrom(msg1, q)

			ctOut := o.SubPlain(ct, pt, true)
			res := e.Decrypt(ctOut)

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloOrder()) {
				for i := range res.Coeffs[0] {
					assert.Equal(t, msgRef, num.Neg(res.Coeffs[0][i], q))
				}
			} else {
				assert.Equal(t, msgRef, res.Coeffs[0][0])
				for i := 1; i < len(res.Coeffs[0]); i++ {
					assert.Equal(t, uint64(0), res.Coeffs[0][i])
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen(), q)
			msg1 := randMsg(pack.PackLen(), q)
			msgRef := vec.Sub(msg0, msg1, q)

			ct := e.Encrypt(pack.Pack(msg0), false)
			pt := pack.Pack(msg1)

			ctOut := o.SubPlain(ct, pt, true)
			res := e.Decrypt(ctOut)
			out := pack.UnPack(res)

			assert.Equal(t, msgRef, out)
		})
	})

	t.Run("SubElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := uint64(rand.Intn(int(q.Value())))
			msg1 := uint64(rand.Intn(int(q.Value())))
			msgRef := num.Sub(msg0, msg1, q)

			el := o.Encoder().Encode(bfv.NewScalarFrom(msg1, q), false, true)
			ct := e.Encrypt(bfv.NewScalarFrom(msg0, q), true)

			ctOut := o.SubElement(ct, el, true)
			res := e.Decrypt(ctOut)

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloOrder()) {
				for i := range res.Coeffs[0] {
					assert.Equal(t, msgRef, num.Neg(res.Coeffs[0][i], q))
				}
			} else {
				assert.Equal(t, msgRef, res.Coeffs[0][0])
				for i := 1; i < len(res.Coeffs[0]); i++ {
					assert.Equal(t, uint64(0), res.Coeffs[0][i])
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen(), q)
			msg1 := randMsg(pack.PackLen(), q)
			msgRef := vec.Sub(msg0, msg1, q)

			ct := e.Encrypt(pack.Pack(msg0), false)
			el := o.Encoder().Encode(pack.Pack(msg1), false, true)

			ctOut := o.SubElement(ct, el, true)
			res := e.Decrypt(ctOut)
			out := pack.UnPack(res)

			assert.Equal(t, msgRef, out)
		})
	})

	t.Run("Mul", func(t *testing.T) {
		msg0 := randMsg(pack.PackLen(), q)
		msg1 := randMsg(pack.PackLen(), q)
		msgRef := vec.Mul(msg0, msg1, q)

		ct0 := e.Encrypt(pack.Pack(msg0), true)
		ct1 := e.Encrypt(pack.Pack(msg1), true)
		rlk := e.NewRelinKey()

		ctOut := o.Mul(ct0, ct1, rlk, true)
		res := e.Decrypt(ctOut)
		out := pack.UnPack(res)

		assert.Equal(t, msgRef, out)
	})

	t.Run("MulPlain", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := uint64(rand.Intn(int(q.Value())))
			msg1 := uint64(rand.Intn(int(q.Value())))
			msgRef := num.Mul(msg0, msg1, q)

			ct := e.Encrypt(bfv.NewScalarFrom(msg0, q), true)
			pt := bfv.NewScalarFrom(msg1, q)

			ctOut := o.MulPlain(ct, pt, true)
			res := e.Decrypt(ctOut)

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloOrder()) {
				for i := range res.Coeffs[0] {
					assert.Equal(t, msgRef, num.Neg(res.Coeffs[0][i], q))
				}
			} else {
				assert.Equal(t, msgRef, res.Coeffs[0][0])
				for i := 1; i < len(res.Coeffs[0]); i++ {
					assert.Equal(t, uint64(0), res.Coeffs[0][i])
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen(), q)
			msg1 := randMsg(pack.PackLen(), q)
			msgRef := vec.Mul(msg0, msg1, q)

			ct := e.Encrypt(pack.Pack(msg0), true)
			pt := pack.Pack(msg1)

			ctOut := o.MulPlain(ct, pt, true)
			res := e.Decrypt(ctOut)
			out := pack.UnPack(res)

			assert.Equal(t, msgRef, out)
		})
	})
}

func TestBFV(t *testing.T) {
	t.Run("type=CyclotomicPow2Mod1", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)

		var q *num.Modulus
		for {
			prime := uint64(4*rSrc.SampleN(1<<20) + 1)
			if num.IsPrime(prime) {
				q = num.NewModulus(prime)
				break
			}
		}

		testOperator(t, rP, q)
	})

	t.Run("type=CyclotomicPow2Mod3", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)

		var q *num.Modulus
		for {
			prime := uint64(4*rSrc.SampleN(1<<20) + 3)
			if num.IsPrime(prime) {
				q = num.NewModulus(prime)
				break
			}
		}

		testOperator(t, rP, q)
	})

	t.Run("type=CyclotomicAny", func(t *testing.T) {
		rP := dft.NewCyclotomicParameters(int(rSrc.SampleN(100) + 1<<11 - 50))
		q := num.NewModulus(dft.MustFindNextNTTPrimes(rP, float64(rSrc.SampleN(1<<20)), 1)[0].Value())

		testOperator(t, rP, q)
	})

	t.Run("type=AutFixedPow2Mod1", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<13, 1<<11)

		var q *num.Modulus
		for {
			prime := uint64(4*rSrc.SampleN(1<<20) + 1)
			if num.IsPrime(prime) {
				q = num.NewModulus(prime)
				break
			}
		}

		testOperator(t, rP, q)
	})

	t.Run("type=AutFixedPow2Mod3", func(t *testing.T) {
		rP := dft.NewAutFixedParameters(1<<13, 1<<11)
		var q *num.Modulus
		for {
			prime := uint64(4*rSrc.SampleN(1<<20) + 3)
			if num.IsPrime(prime) {
				q = num.NewModulus(prime)
				break
			}
		}

		testOperator(t, rP, q)
	})

	t.Run("type=AutFixedPrime", func(t *testing.T) {
		prime := num.MustNextPrime(int(rSrc.SampleN(20)+1), 1)

		var M int
		var N int
		for {
			M = num.MustNextPrime(int(rSrc.SampleN(1<<10)), 1)
			ord := num.Order(uint64(prime), num.NewModulus(M))
			if ord < 10 {
				N = (M - 1) / int(ord)
				break
			}
		}

		q := num.NewModulus(num.Exp(uint64(prime), 3, nil))
		rP := dft.NewAutFixedParameters(M, N)

		testOperator(t, rP, q)
	})
}
