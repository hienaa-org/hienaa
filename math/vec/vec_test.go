package vec_test

import (
	"math/rand"
	"slices"
	"testing"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = rand.New(rand.NewSource(0))
)

func TestOps(t *testing.T) {
	q := num.NewModulus(((rSrc.Uint64()%num.MaxModulus)>>1)<<1 + 1)

	N := rSrc.Int() % (1 << 16)
	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)
	vOutCheck := make([]uint64, N)
	vOutInit := make([]uint64, N)

	for i := 0; i < N; i++ {
		v0[i] = rSrc.Uint64() % q.Value()
		v1[i] = rSrc.Uint64() % q.Value()
		vOutInit[i] = rSrc.Uint64() % q.Value()
	}

	t.Run("Add", func(t *testing.T) {
		vec.AddTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] -= q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("AddLazy", func(t *testing.T) {
		vec.AddLazyTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("Sub", func(t *testing.T) {
		vec.SubTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] += q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("SubLazy", func(t *testing.T) {
		vec.SubLazyTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarBMul", func(t *testing.T) {
		vec.ScalarBMulTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.BMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarBMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarBMulAddTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.BMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarBMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarBMulSubTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.BMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarBMulLazy", func(t *testing.T) {
		vec.ScalarBMulLazyTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.BMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.BMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarBMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarBMulAddLazyTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.BMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.BMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarBMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarBMulSubLazyTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.BMulLazy(v0[i], q.Value()-v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.BMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("BMul", func(t *testing.T) {
		vec.BMulTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.BMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("BMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.BMulAddTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.BMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("BMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.BMulSubTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.BMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("BMulLazy", func(t *testing.T) {
		vec.BMulLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.BMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.BMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("BMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.BMulAddLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.BMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.BMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("BMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.BMulSubLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.BMulLazy(q.Value()-v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.BMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMul", func(t *testing.T) {
		vec.ScalarMMulTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulAddTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulSubTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulLazy", func(t *testing.T) {
		vec.ScalarMMulLazyTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulAddLazyTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulSubLazyTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(v0[i], q.Value()-v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMul", func(t *testing.T) {
		vec.MMulTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulAddTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulSubTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulLazy", func(t *testing.T) {
		vec.MMulLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulAddLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulSubLazyTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(q.Value()-v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
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
			vec.AddTo(v0, v1, q, vOut)
		}
	})

	b.Run("AddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.AddLazyTo(v0, v1, vOut)
		}
	})

	b.Run("Sub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.SubTo(v0, v1, q, vOut)
		}
	})

	b.Run("SubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.SubLazyTo(v0, v1, vOut)
		}
	})

	b.Run("ScalarBMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarBMulTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarBMulAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarBMulAddTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarBMulSub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarBMulSubTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarBMulLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarBMulLazyTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarBMulAddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarBMulAddLazyTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarBMulSubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarBMulSubLazyTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("BMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.BMulTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.BMulAddTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulSub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.BMulSubTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.BMulLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulAddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.BMulAddLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("BMulSubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.BMulSubLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("ScalarMMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarMMulTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarMMulAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarMMulAddTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarMMulSub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarMMulSubTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarMMulLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarMMulLazyTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarMMulAddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarMMulAddLazyTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("ScalarMMulSubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ScalarMMulSubLazyTo(v0, v1[0], q, vOut)
		}
	})

	b.Run("MMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MMulTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMulAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MMulAddTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMulSub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MMulSubTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMulLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MMulLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMulAddLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MMulAddLazyTo(v0, v1, q, vOut)
		}
	})

	b.Run("MMulSubLazy", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MMulSubLazyTo(v0, v1, q, vOut)
		}
	})
}
