package mod_test

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = rand.New(rand.NewSource(0))
)

func TestReduce(t *testing.T) {
	q := mod.NewModulus((rSrc.Uint64()>>10)<<1 + 1)

	x64 := rSrc.Uint64()
	x128Hi := rSrc.Uint64()
	x128Lo := rSrc.Uint64()
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
	q := mod.NewModulus((rSrc.Uint64()>>4)<<1 + 1)
	x0 := rSrc.Uint64() % q.Value()
	x1 := rSrc.Uint64() % q.Value()

	x0Big := new(big.Int).SetUint64(x0)
	x1Big := new(big.Int).SetUint64(x1)
	qBig := new(big.Int).SetUint64(q.Value())

	xAddBig := new(big.Int).Add(x0Big, x1Big)
	xAddBig.Mod(xAddBig, qBig)

	xSubBig := new(big.Int).Sub(x0Big, x1Big)
	xSubBig.Mod(xSubBig, qBig)

	xMulBig := new(big.Int).Mul(x0Big, x1Big)
	xMulBig.Mod(xMulBig, qBig)

	t.Run("Add", func(t *testing.T) {
		xAdd := mod.Add(x0, x1, q)
		assert.Equal(t, xAdd, xAddBig.Uint64())
	})

	t.Run("Sub", func(t *testing.T) {
		xSub := mod.Sub(x0, x1, q)
		assert.Equal(t, xSub, xSubBig.Uint64())
	})

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
}

func BenchmarkMul(b *testing.B) {
	N := 1 << 15
	q := mod.NewModulus(rSrc.Uint64() >> 3)

	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)
	for i := 0; i < N; i++ {
		v0[i] = rSrc.Uint64()
		v1[i] = rSrc.Uint64()
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
