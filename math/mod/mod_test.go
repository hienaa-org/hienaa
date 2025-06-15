package mod_test

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = csprng.NewUniformSamplerWithSeed(nil)
)

func TestReduce(t *testing.T) {
	q := mod.NewModulus(rSrc.SampleN(mod.MaxModulus) | 1)

	x64 := rSrc.Sample()
	x128Hi := rSrc.Sample()
	x128Lo := rSrc.Sample()
	x128 := new(big.Int).Lsh(new(big.Int).SetUint64(x128Hi), 64)
	x128.Add(x128, new(big.Int).SetUint64(x128Lo))

	t.Run("Reduce", func(t *testing.T) {
		assert.Equal(t, mod.Reduce(x64, q), x64%q.Value())
	})

	t.Run("Reduce128", func(t *testing.T) {
		x128.Mod(x128, new(big.Int).SetUint64(q.Value()))
		assert.Equal(t, mod.Reduce128(x128Hi, x128Lo, q), x128.Uint64())
	})
}

func TestOps(t *testing.T) {
	qBig, err := rand.Prime(rSrc, mod.MaxModulusBits)
	assert.NoError(t, err)

	q := mod.NewModulus(qBig.Uint64())
	x0 := rSrc.SampleN(q.Value())
	x1 := rSrc.SampleN(q.Value())

	x0Big := new(big.Int).SetUint64(x0)
	x1Big := new(big.Int).SetUint64(x1)

	t.Run("Add", func(t *testing.T) {
		xAdd := mod.Add(x0, x1, q)
		xAddBig := new(big.Int).Add(x0Big, x1Big)
		xAddBig.Mod(xAddBig, qBig)
		assert.Equal(t, xAdd, xAddBig.Uint64())
	})

	t.Run("Sub", func(t *testing.T) {
		xSub := mod.Sub(x0, x1, q)
		xSubBig := new(big.Int).Sub(x0Big, x1Big)
		xSubBig.Mod(xSubBig, qBig)
		assert.Equal(t, xSub, xSubBig.Uint64())
	})

	xMulBig := new(big.Int).Mul(x0Big, x1Big)
	xMulBig.Mod(xMulBig, qBig)

	t.Run("Barrett", func(t *testing.T) {
		xMul := mod.Mul(x0, x1, q)
		assert.Equal(t, xMul, xMulBig.Uint64())
	})

	t.Run("Montgomery", func(t *testing.T) {
		x0M := mod.MForm(x0, q)
		x1M := mod.MForm(x1, q)
		xMulM := mod.MMul(x0M, x1M, q)
		xMul := mod.InvMForm(xMulM, q)
		assert.Equal(t, xMul, xMulBig.Uint64())
	})

	t.Run("Exp", func(t *testing.T) {
		xExp := mod.Exp(x0, x1, q)
		xExpBig := new(big.Int).Exp(x0Big, x1Big, qBig)
		assert.Equal(t, xExp, xExpBig.Uint64())
	})

	t.Run("Inv", func(t *testing.T) {
		xInv := mod.Inv(x0, q)
		xInvBig := new(big.Int).ModInverse(x0Big, qBig)
		assert.Equal(t, xInv, xInvBig.Uint64())
	})
}

func BenchmarkMul(b *testing.B) {
	N := 1 << 15
	q := mod.NewModulus(rSrc.SampleN(mod.MaxModulus) | 1)

	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)
	for i := 0; i < N; i++ {
		v0[i] = rSrc.SampleN(q.Value())
		v1[i] = rSrc.SampleN(q.Value())
	}

	b.Run("Barrett", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < N; j++ {
				vOut[j] = mod.Mul(v0[j], v1[j], q)
			}
		}
	})

	b.Run("BarrettLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < N; j++ {
				vOut[j] = mod.MulLazy(v0[j], v1[j], q)
			}
		}
	})

	b.Run("Montgomery", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < N; j++ {
				vOut[j] = mod.MMul(v0[j], v1[j], q)
			}
		}
	})

	b.Run("MontgomeryLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < N; j++ {
				vOut[j] = mod.MMulLazy(v0[j], v1[j], q)
			}
		}
	})
}
