package mod_test

import (
	"fmt"
	"slices"
	"testing"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/stretchr/testify/assert"
)

var (
	benchLogN = []int{12, 13, 14, 15, 16, 17}
)

func TestVecOps(t *testing.T) {
	q := mod.NewModulus(rSrc.SampleN(mod.MaxModulus) | 1)

	N := int(rSrc.SampleN(1 << 16))
	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)
	vOutCheck := make([]uint64, N)
	vOutInit := make([]uint64, N)

	for i := 0; i < N; i++ {
		v0[i] = rSrc.SampleN(q.Value())
		v1[i] = rSrc.SampleN(q.Value())
		vOutInit[i] = rSrc.SampleN(q.Value())
	}

	t.Run("Add", func(t *testing.T) {
		mod.AddVecTo(vOut, v0, v1, q)
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
		mod.AddLazyVecTo(vOut, v0, v1)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("Sub", func(t *testing.T) {
		mod.SubVecTo(vOut, v0, v1, q)
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
		mod.SubLazyVecTo(vOut, v0, v1)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMul", func(t *testing.T) {
		mod.ScalarMulVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMulAddVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMulSubVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMulLazy", func(t *testing.T) {
		mod.ScalarMulLazyVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.SMulLazy(v0[i], v1[0], mod.SForm(v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMulAddLazyVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.SMulLazy(v0[i], v1[0], mod.SForm(v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Add(vOutInit[i], mod.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMulSubLazyVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.SMulLazy(v0[i], q.Value()-v1[0], mod.SForm(q.Value()-v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Sub(vOutInit[i], mod.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMul", func(t *testing.T) {
		mod.ScalarMMulVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMMulAddVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMMulSubVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulLazy", func(t *testing.T) {
		mod.ScalarMMulLazyVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMMulAddLazyVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.MMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Add(vOutInit[i], mod.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMMulSubLazyVecTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.MMulLazy(v0[i], q.Value()-v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Sub(vOutInit[i], mod.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("Mul", func(t *testing.T) {
		mod.MulVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MulAddVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("MulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MulSubVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("MulLazy", func(t *testing.T) {
		mod.MulLazyVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MulAddLazyVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.MulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Add(vOutInit[i], mod.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MulSubLazyVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.MulLazy(q.Value()-v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Sub(vOutInit[i], mod.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMul", func(t *testing.T) {
		mod.MMulVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MMulAddVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MMulSubVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulLazy", func(t *testing.T) {
		mod.MMulLazyVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MMulAddLazyVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.MMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Add(vOutInit[i], mod.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MMulSubLazyVecTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.MMulLazy(q.Value()-v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Sub(vOutInit[i], mod.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	v1S := mod.SFormVec(v1, q)

	t.Run("SMul", func(t *testing.T) {
		mod.SMulVecTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("SMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.SMulAddVecTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("SMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.SMulSubVecTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("SMulLazy", func(t *testing.T) {
		mod.SMulLazyVecTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.SMulLazy(v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.SMul(v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("SMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.SMulAddLazyVecTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.SMulLazy(v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Add(vOutInit[i], mod.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("SMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.SMulSubLazyVecTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += mod.SMulLazy(q.Value()-v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = mod.Sub(vOutInit[i], mod.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("Reduce", func(t *testing.T) {
		for i := 0; i < N; i++ {
			v0[i] = rSrc.Sample()
		}

		mod.ReduceVecTo(vOut, v0, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] % q.Value()
		}
		assert.Equal(t, vOutCheck, vOut)
	})
}

func BenchmarkVecOps(b *testing.B) {
	q := mod.NewModulus(rSrc.SampleN(mod.MaxModulus) | 1)

	for _, logN := range benchLogN {
		N := 1 << logN
		v0 := make([]uint64, N)
		v1 := make([]uint64, N)
		v1S := make([]uint64, N)
		vOut := make([]uint64, N)

		for i := 0; i < N; i++ {
			v0[i] = rSrc.SampleN(q.Value())
			v1[i] = rSrc.SampleN(q.Value())
			v1S[i] = rSrc.SampleN(q.Value())
		}

		b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
			b.Run("Add", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.AddVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("AddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.AddLazyVecTo(vOut, v0, v1)
				}
			})

			b.Run("Sub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SubVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("SubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SubLazyVecTo(vOut, v0, v1)
				}
			})

			b.Run("ScalarMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulAddVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulSubVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulLazyVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulAddLazyVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulSubLazyVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulAddVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulSubVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulLazyVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulAddLazyVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulSubLazyVecTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("Mul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulAddVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulSubVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulLazyVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulAddLazyVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulSubLazyVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MForm", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MFormVecTo(vOut, v0, q)
				}
			})

			b.Run("InvMForm", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.InvMFormVecTo(vOut, v0, q)
				}
			})

			b.Run("MMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulAddVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulSubVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulLazyVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulAddLazyVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulSubLazyVecTo(vOut, v0, v1, q)
				}
			})

			b.Run("SForm", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SFormVecTo(vOut, v0, q)
				}
			})

			b.Run("SMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SMulVecTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SMulAddVecTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SMulSubVecTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SMulLazyVecTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SMulAddLazyVecTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SMulSubLazyVecTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("Reduce", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ReduceVecTo(vOut, v0, q)
				}
			})
		})
	}
}
