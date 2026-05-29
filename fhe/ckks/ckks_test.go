package ckks_test

import (
	"math"
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/ckks"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
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

func randMsg(packLen int) []complex128 {
	msg := make([]complex128, packLen)
	for i := range msg {
		msg[i] = complex(rand.Float64(), rand.Float64())
	}
	return msg
}

func checkApprox(v0 []complex128, v1 []complex128, tol float64) bool {
	for i := range v0 {
		if math.Abs(real(v0[i])-real(v1[i])) > tol || math.Abs(imag(v0[i])-imag(v1[i])) > tol {
			return false
		}
	}
	return true
}

func testOperator(t *testing.T, rP dft.RingParameters, scFac float64) {
	baseMod, auxMod := rlwe.FindNTTPrimes(rP, baseModBits, auxModBits)
	p := rlwe.ParametersLiteral{
		RingParams:      rP,
		BaseModulus:     baseMod,
		AuxModulus:      auxMod,
		SecretKeyParams: skParams,
		NoiseParams:     noiseParams,
	}.Compile()

	o := ckks.NewOperator(p, scFac)
	enc := ckks.NewEncryptor(p)
	ecd := ckks.NewEncoder(p)
	pack := ckks.NewComplexPacker(p)

	t.Run("Encrypt", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg := rand.Float64()

			ct := enc.Encrypt(ckks.NewScalarFrom(msg, ckks.TypeReal), scFac, true)
			res := enc.Decrypt(ct)

			tol := p.NoiseParams().Bound() / scFac * float64(rP.Rank())
			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloIndex()) {
				// TODO: Implement
			} else {
				assert.True(t, math.Abs(msg-res.Value[0]) < tol)
				for i := 1; i < res.Rank(); i++ {
					assert.True(t, math.Abs(res.Value[i]) < tol)
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg := randMsg(pack.PackLen())
			pt := pack.Pack(msg)

			ct := enc.Encrypt(pt, scFac, true)
			res := enc.Decrypt(ct)
			out := pack.UnPack(res)

			tol := p.NoiseParams().Bound() / scFac * float64(rP.Rank())
			assert.True(t, checkApprox(msg, out, tol))
		})
	})

	t.Run("Add", func(t *testing.T) {
		msg0 := randMsg(pack.PackLen())
		msg1 := randMsg(pack.PackLen())
		msgRef := make([]complex128, pack.PackLen())
		for i := range msg0 {
			msgRef[i] = msg0[i] + msg1[i]
		}

		pt0 := pack.Pack(msg0)
		pt1 := pack.Pack(msg1)

		ct0 := enc.Encrypt(pt0, scFac, true)
		ct1 := enc.Encrypt(pt1, scFac, true)

		ctOut := o.Add(ct0, ct1, true)
		res := enc.Decrypt(ctOut)
		out := pack.UnPack(res)

		tol := p.NoiseParams().Bound() / scFac * float64(rP.Rank()) * 2
		assert.True(t, checkApprox(msgRef, out, tol))
	})

	t.Run("AddPlain", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := rand.Float64()
			msg1 := rand.Float64()
			msgRef := msg0 + msg1

			ct := enc.Encrypt(ckks.NewScalarFrom(msg0, ckks.TypeReal), scFac, true)
			pt := ckks.NewScalarFrom(msg1, ckks.TypeReal)

			ctOut := o.AddPlain(ct, pt, true)
			res := enc.Decrypt(ctOut)
			tol := p.NoiseParams().Bound()/scFac*float64(rP.Rank()) + 1/scFac

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloIndex()) {
				// TODO: Implement
			} else {
				assert.True(t, math.Abs(msgRef-res.Value[0]) < tol)
				for i := 1; i < res.Rank(); i++ {
					assert.True(t, math.Abs(res.Value[i]) < tol)
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen())
			msg1 := randMsg(pack.PackLen())
			msgRef := make([]complex128, pack.PackLen())
			for i := range msg0 {
				msgRef[i] = msg0[i] + msg1[i]
			}

			pt0 := pack.Pack(msg0)
			pt1 := pack.Pack(msg1)

			ct := enc.Encrypt(pt0, scFac, true)
			ctOut := o.AddPlain(ct, pt1, true)
			res := enc.Decrypt(ctOut)
			out := pack.UnPack(res)

			tol := (p.NoiseParams().Bound() + 1) / scFac * float64(rP.Rank())
			assert.True(t, checkApprox(msgRef, out, tol))
		})
	})

	t.Run("AddElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := rand.Float64()
			msg1 := rand.Float64()
			msgRef := msg0 + msg1

			ct := enc.Encrypt(ckks.NewScalarFrom(msg0, ckks.TypeReal), scFac, true)
			el := ecd.Encode(ckks.NewScalarFrom(msg1, ckks.TypeReal), scFac, false, true)

			ctOut := o.AddElement(ct, el, true)
			res := enc.Decrypt(ctOut)
			tol := p.NoiseParams().Bound()/scFac*float64(rP.Rank()) + 1/scFac

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloIndex()) {
				// TODO: Implement
			} else {
				assert.True(t, math.Abs(msgRef-res.Value[0]) < tol)
				for i := 1; i < res.Rank(); i++ {
					assert.True(t, math.Abs(res.Value[i]) < tol)
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen())
			msg1 := randMsg(pack.PackLen())
			msgRef := make([]complex128, pack.PackLen())
			for i := range msg0 {
				msgRef[i] = msg0[i] + msg1[i]
			}

			pt0 := pack.Pack(msg0)
			pt1 := pack.Pack(msg1)

			ct := enc.Encrypt(pt0, scFac, true)
			el := ecd.Encode(pt1, scFac, false, true)

			ctOut := o.AddElement(ct, el, true)
			res := enc.Decrypt(ctOut)
			out := pack.UnPack(res)

			tol := (p.NoiseParams().Bound() + 1) / scFac * float64(rP.Rank())
			assert.True(t, checkApprox(msgRef, out, tol))
		})
	})

	t.Run("Sub", func(t *testing.T) {
		msg0 := randMsg(pack.PackLen())
		msg1 := randMsg(pack.PackLen())
		msgRef := make([]complex128, pack.PackLen())
		for i := range msg0 {
			msgRef[i] = msg0[i] - msg1[i]
		}

		pt0 := pack.Pack(msg0)
		pt1 := pack.Pack(msg1)

		ct0 := enc.Encrypt(pt0, scFac, true)
		ct1 := enc.Encrypt(pt1, scFac, true)

		ctOut := o.Sub(ct0, ct1, true)
		res := enc.Decrypt(ctOut)
		out := pack.UnPack(res)

		tol := p.NoiseParams().Bound() / scFac * float64(rP.Rank())
		assert.True(t, checkApprox(msgRef, out, tol))
	})

	t.Run("SubPlain", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := rand.Float64()
			msg1 := rand.Float64()
			msgRef := msg0 - msg1

			ct := enc.Encrypt(ckks.NewScalarFrom(msg0, ckks.TypeReal), scFac, true)
			pt := ckks.NewScalarFrom(msg1, ckks.TypeReal)

			ctOut := o.SubPlain(ct, pt, true)
			res := enc.Decrypt(ctOut)
			tol := p.NoiseParams().Bound()/scFac*float64(rP.Rank()) + 1/scFac

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloIndex()) {
				// TODO: Implement
			} else {
				assert.True(t, math.Abs(msgRef-res.Value[0]) < tol)
				for i := 1; i < res.Rank(); i++ {
					assert.True(t, math.Abs(res.Value[i]) < tol)
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen())
			msg1 := randMsg(pack.PackLen())
			msgRef := make([]complex128, pack.PackLen())
			for i := range msg0 {
				msgRef[i] = msg0[i] - msg1[i]
			}

			pt0 := pack.Pack(msg0)
			pt1 := pack.Pack(msg1)

			ct := enc.Encrypt(pt0, scFac, true)
			ctOut := o.SubPlain(ct, pt1, true)
			res := enc.Decrypt(ctOut)
			out := pack.UnPack(res)

			tol := (p.NoiseParams().Bound() + 1) / scFac * float64(rP.Rank())
			assert.True(t, checkApprox(msgRef, out, tol))
		})
	})

	t.Run("SubElement", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := rand.Float64()
			msg1 := rand.Float64()
			msgRef := msg0 - msg1

			ct := enc.Encrypt(ckks.NewScalarFrom(msg0, ckks.TypeReal), scFac, true)
			el := ecd.Encode(ckks.NewScalarFrom(msg1, ckks.TypeReal), scFac, false, true)

			ctOut := o.SubElement(ct, el, true)
			res := enc.Decrypt(ctOut)
			tol := p.NoiseParams().Bound()/scFac*float64(rP.Rank()) + 1/scFac

			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloIndex()) {
				// TODO: Implement
			} else {
				assert.True(t, math.Abs(msgRef-res.Value[0]) < tol)
				for i := 1; i < res.Rank(); i++ {
					assert.True(t, math.Abs(res.Value[i]) < tol)
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen())
			msg1 := randMsg(pack.PackLen())
			msgRef := make([]complex128, pack.PackLen())
			for i := range msg0 {
				msgRef[i] = msg0[i] - msg1[i]
			}

			pt0 := pack.Pack(msg0)
			pt1 := pack.Pack(msg1)

			ct := enc.Encrypt(pt0, scFac, true)
			el := ecd.Encode(pt1, scFac, false, true)

			ctOut := o.SubElement(ct, el, true)
			res := enc.Decrypt(ctOut)
			out := pack.UnPack(res)

			tol := (p.NoiseParams().Bound() + 1) / scFac * float64(rP.Rank())
			assert.True(t, checkApprox(msgRef, out, tol))
		})
	})

	t.Run("Mul", func(t *testing.T) {
		msg0 := randMsg(pack.PackLen())
		msg1 := randMsg(pack.PackLen())
		msgRef := make([]complex128, pack.PackLen())
		for i := range msg0 {
			msgRef[i] = msg0[i] * msg1[i]
		}

		pt0 := pack.Pack(msg0)
		pt1 := pack.Pack(msg1)

		ct0 := enc.Encrypt(pt0, scFac, true)
		ct1 := enc.Encrypt(pt1, scFac, true)
		rlk := enc.NewRelinKey()

		ctOut := o.Mul(ct0, ct1, rlk, true)
		res := enc.Decrypt(ctOut)
		out := pack.UnPack(res)

		// TODO: More precise noise bound.
		tol := float64(rP.Rank()) / scFac
		assert.True(t, checkApprox(msgRef, out, tol))
	})

	t.Run("MulPlain", func(t *testing.T) {
		t.Run("Scalar", func(t *testing.T) {
			msg0 := rand.Float64()
			msg1 := rand.Float64()
			msgRef := msg0 * msg1

			ct := enc.Encrypt(ckks.NewScalarFrom(msg0, ckks.TypeReal), scFac, true)
			pt := ckks.NewScalarFrom(msg1, ckks.TypeReal)

			ctOut := o.MulPlain(ct, pt, true)
			res := enc.Decrypt(ctOut)

			// TODO: More precise noise bound.
			tol := float64(rP.Rank()) / scFac
			if rP.RingType() == dft.TypeAutFixed && num.IsPrime(rP.CycloIndex()) {
				// TODO: Implement
			} else {
				assert.True(t, math.Abs(msgRef-res.Value[0]) < tol)
				for i := 1; i < res.Rank(); i++ {
					assert.True(t, math.Abs(res.Value[i]) < tol)
				}
			}
		})

		t.Run("Poly", func(t *testing.T) {
			msg0 := randMsg(pack.PackLen())
			msg1 := randMsg(pack.PackLen())
			msgRef := make([]complex128, pack.PackLen())
			for i := range msg0 {
				msgRef[i] = msg0[i] * msg1[i]
			}

			pt0 := pack.Pack(msg0)
			pt1 := pack.Pack(msg1)
			ct := enc.Encrypt(pt0, scFac, true)
			ctOut := o.MulPlain(ct, pt1, true)
			res := enc.Decrypt(ctOut)
			out := pack.UnPack(res)

			// TODO: More precise noise bound.
			tol := float64(rP.Rank()) / scFac
			assert.True(t, checkApprox(msgRef, out, tol))
		})
	})
}

func TestCKKS(t *testing.T) {
	t.Run("type=CyclotomicPow2", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)
		scFac := float64(rSrc.SampleN(1<<10)<<10) * float64(N)

		testOperator(t, rP, scFac)
	})
}
