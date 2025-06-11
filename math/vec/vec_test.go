package vec

import (
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = rand.New(rand.NewSource(0))
)

func TestOps(t *testing.T) {
	q := num.NewModulus(rSrc.Uint64() >> 3)

	N := 1<<15 + 3
	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)
	vOutCheck := make([]uint64, N)

	for i := 0; i < N; i++ {
		v0[i] = rSrc.Uint64() % q.Value()
		v1[i] = rSrc.Uint64() % q.Value()
	}

	t.Run("Add", func(t *testing.T) {
		AddTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] -= q.Value()
			}
		}
		assert.EqualValues(t, vOut, vOutCheck)
	})

	t.Run("AddLazy", func(t *testing.T) {
		AddLazyTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.EqualValues(t, vOut, vOutCheck)
	})

	t.Run("Sub", func(t *testing.T) {
		SubTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] += q.Value()
			}
		}
		assert.EqualValues(t, vOut, vOutCheck)
	})

	t.Run("SubLazy", func(t *testing.T) {
		SubLazyTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.EqualValues(t, vOut, vOutCheck)
	})
}

func BenchmarkOps(b *testing.B) {
	q := num.NewModulus(rSrc.Uint64() >> 3)

	N := 1 << 15
	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)

	for i := 0; i < N; i++ {
		v0[i] = rSrc.Uint64() % q.Value()
		v1[i] = rSrc.Uint64() % q.Value()
	}

	b.Run("Add", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			AddTo(v0, v1, q, vOut)
		}
	})

	b.Run("AddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			AddLazyTo(v0, v1, vOut)
		}
	})

	b.Run("Sub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			SubTo(v0, v1, q, vOut)
		}
	})

	b.Run("SubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			SubLazyTo(v0, v1, vOut)
		}
	})
}
