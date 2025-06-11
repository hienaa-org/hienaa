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
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("AddLazy", func(t *testing.T) {
		AddLazyTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("Sub", func(t *testing.T) {
		SubTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] += q.Value()
			}
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("SubLazy", func(t *testing.T) {
		SubLazyTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("BMul", func(t *testing.T) {
		BMulTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.BMul(v0[i], v1[i], q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("BMulAdd", func(t *testing.T) {
		copy(vOut, vOutCheck)

		BMulAddTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.BMul(v0[i], v1[i], q), q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("BMulSub", func(t *testing.T) {
		copy(vOut, vOutCheck)

		BMulSubTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.BMul(v0[i], v1[i], q), q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("BMulLazy", func(t *testing.T) {
		BMulLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.BMulLazy(v0[i], v1[i], q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("BMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutCheck)

		BMulAddLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.BMulLazy(v0[i], v1[i], q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("BMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutCheck)

		BMulSubLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += q.Value()<<1 - num.BMulLazy(v0[i], v1[i], q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
	})

	t.Run("MMul", func(t *testing.T) {
		MMulTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMul(v0[i], v1[i], q)
		}
		assert.EqualValues(t, vOutCheck, vOut)
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

	b.Run("BMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BMulTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BMulAddTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulSub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BMulSubTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BMulLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulAddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BMulAddLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulSubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BMulSubLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			MMulTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMulLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			MMulLazyTo(v0, v1, q, vOut)
		}
	})
}
