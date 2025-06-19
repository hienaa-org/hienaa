package mod_test

import (
	"fmt"
	"slices"
	"testing"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/stretchr/testify/assert"
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
		mod.AddVecTo(v0, v1, q, vOut)
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
		mod.AddLazyVecTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("Sub", func(t *testing.T) {
		mod.SubVecTo(v0, v1, q, vOut)
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
		mod.SubLazyVecTo(v0, v1, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMul", func(t *testing.T) {
		mod.ScalarMulVecTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMulAddVecTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMulSubVecTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMulLazy", func(t *testing.T) {
		mod.ScalarMulLazyVecTo(v0, v1[0], q, vOut)
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

		mod.ScalarMulAddLazyVecTo(v0, v1[0], q, vOut)
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

		mod.ScalarMulSubLazyVecTo(v0, v1[0], q, vOut)
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
		mod.ScalarMMulVecTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMMulAddVecTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.ScalarMMulSubVecTo(v0, v1[0], q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulLazy", func(t *testing.T) {
		mod.ScalarMMulLazyVecTo(v0, v1[0], q, vOut)
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

		mod.ScalarMMulAddLazyVecTo(v0, v1[0], q, vOut)
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

		mod.ScalarMMulSubLazyVecTo(v0, v1[0], q, vOut)
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
		mod.MulVecTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MulAddVecTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("MulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MulSubVecTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), 2*q.Value())
	})

	t.Run("MulLazy", func(t *testing.T) {
		mod.MulLazyVecTo(v0, v1, q, vOut)
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

		mod.MulAddLazyVecTo(v0, v1, q, vOut)
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

		mod.MulSubLazyVecTo(v0, v1, q, vOut)
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
		mod.MMulVecTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MMulAddVecTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Add(vOutCheck[i], mod.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		mod.MMulSubVecTo(v0, v1, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = mod.Sub(vOutCheck[i], mod.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, slices.Max(vOut), q.Value())
	})

	t.Run("MMulLazy", func(t *testing.T) {
		mod.MMulLazyVecTo(v0, v1, q, vOut)
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

		mod.MMulAddLazyVecTo(v0, v1, q, vOut)
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

		mod.MMulSubLazyVecTo(v0, v1, q, vOut)
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

	t.Run("Reduce", func(t *testing.T) {
		for i := 0; i < N; i++ {
			v0[i] = rSrc.Sample()
		}

		mod.ReduceVecTo(v0, q, vOut)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] % q.Value()
		}
		assert.Equal(t, vOutCheck, vOut)
	})
}

func BenchmarkVecOps(b *testing.B) {
	q := mod.NewModulus(rSrc.SampleN(mod.MaxModulus) | 1)

	for _, logN := range []int{12, 13, 14, 15, 16, 17} {
		N := 1 << logN
		v0 := make([]uint64, N)
		v1 := make([]uint64, N)
		vOut := make([]uint64, N)

		for i := 0; i < N; i++ {
			v0[i] = rSrc.SampleN(q.Value())
			v1[i] = rSrc.SampleN(q.Value())
		}

		b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
			b.Run("Add", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.AddVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("AddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.AddLazyVecTo(v0, v1, vOut)
				}
			})

			b.Run("Sub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SubVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("SubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.SubLazyVecTo(v0, v1, vOut)
				}
			})

			b.Run("ScalarMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulAddVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulSubVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulLazyVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulAddLazyVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMulSubLazyVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulAddVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulSubVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulLazyVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulAddLazyVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("ScalarMMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ScalarMMulSubLazyVecTo(v0, v1[0], q, vOut)
				}
			})

			b.Run("Mul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulAddVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulSubVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulLazyVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulAddLazyVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MulSubLazyVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulAddVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulSubVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulLazyVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulAddLazyVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("MMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.MMulSubLazyVecTo(v0, v1, q, vOut)
				}
			})

			b.Run("Reduce", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					mod.ReduceVecTo(v0, q, vOut)
				}
			})
		})
	}
}
