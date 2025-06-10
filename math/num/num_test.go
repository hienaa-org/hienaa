package num_test

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = rand.New(rand.NewSource(0))
)

func TestOps(t *testing.T) {
	q := num.NewModulus((rSrc.Uint64()>>4)<<1 + 1)
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
		xAdd := num.Add(x0, x1, q)
		assert.Equal(t, xAdd, xAddBig.Uint64())
	})

	t.Run("Sub", func(t *testing.T) {
		xSub := num.Sub(x0, x1, q)
		assert.Equal(t, xSub, xSubBig.Uint64())
	})

	t.Run("Barrett", func(t *testing.T) {
		xMul := num.BMul(x0, x1, q)
		assert.Equal(t, xMul, xMulBig.Uint64())
	})

	t.Run("Montgomery", func(t *testing.T) {
		x0M := num.MForm(x0, q)
		x1M := num.MForm(x1, q)
		xMulM := num.MMul(x0M, x1M, q)
		xMul := num.InvMForm(xMulM, q)
		assert.Equal(t, xMul, xMulBig.Uint64())
	})
}

func BenchmarkMul(b *testing.B) {
	N := 1 << 15
	q := num.NewModulus(rSrc.Uint64() >> 3)

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
				vOut[j] = num.BMul(v0[j], v1[j], q)
			}
		}
	})

	b.Run("Montgomery", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := 0; j < N; j++ {
				vOut[j] = num.MMul(v0[j], v1[j], q)
			}
		}
	})
}
