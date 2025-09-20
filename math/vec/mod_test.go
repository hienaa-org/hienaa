package vec_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{12, 13, 14, 15, 16, 17}
)

func TestOps(t *testing.T) {
	q := num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)

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
		vec.AddTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] -= q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("AddLazy", func(t *testing.T) {
		vec.AddLazyTo(vOut, v0, v1)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())
	})

	t.Run("ScalarAdd", func(t *testing.T) {
		vec.ScalarAddTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[0]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] -= q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarAddLazy", func(t *testing.T) {
		vec.ScalarAddLazyTo(vOut, v0, v1[0])
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[0]
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())
	})

	t.Run("Sub", func(t *testing.T) {
		vec.SubTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] += q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("SubLazy", func(t *testing.T) {
		vec.SubLazyTo(vOut, v0, v1)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarSub", func(t *testing.T) {
		vec.ScalarSubTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[0]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] += q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarSubLazy", func(t *testing.T) {
		vec.ScalarSubLazyTo(vOut, v0, v1[0])
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[0]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("Neg", func(t *testing.T) {
		vec.NegTo(vOut, v0, q)
		for i := 0; i < N; i++ {
			if v0[i] == 0 {
				vOutCheck[i] = 0
			} else {
				vOutCheck[i] = q.Value() - v0[i]
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMul", func(t *testing.T) {
		vec.ScalarMulTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMulAddTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMulSubTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMulLazy", func(t *testing.T) {
		vec.ScalarMulLazyTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.SMulLazy(v0[i], v1[0], num.SForm(v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMulAddLazyTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.SMulLazy(v0[i], v1[0], num.SForm(v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMulSubLazyTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.SMulLazy(v0[i], q.Value()-v1[0], num.SForm(q.Value()-v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMul", func(t *testing.T) {
		vec.ScalarMMulTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulAddTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulSubTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("ScalarMMulLazy", func(t *testing.T) {
		vec.ScalarMMulLazyTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.MMul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulAddLazyTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("ScalarMMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.ScalarMMulSubLazyTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(v0[i], q.Value()-v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.MMul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("Mul", func(t *testing.T) {
		vec.MulTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulAddTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulSubTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MulLazy", func(t *testing.T) {
		vec.MulLazyTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulAddLazyTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulSubLazyTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MulLazy(q.Value()-v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMul", func(t *testing.T) {
		vec.MMulTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulAddTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulSubTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MMulLazy", func(t *testing.T) {
		vec.MMulLazyTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.MMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.MMul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulAddLazyTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MMulSubLazyTo(vOut, v0, v1, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.MMulLazy(q.Value()-v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.MMul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	v1S := vec.SForm(v1, q)

	t.Run("SMul", func(t *testing.T) {
		vec.SMulTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("SMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.SMulAddTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("SMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.SMulSubTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("SMulLazy", func(t *testing.T) {
		vec.SMulLazyTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.SMulLazy(v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.SMul(v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("SMulAddLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.SMulAddLazyTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.SMulLazy(v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Add(vOutInit[i], num.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("SMulSubLazy", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.SMulSubLazyTo(vOut, v0, v1, v1S, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] += num.SMulLazy(q.Value()-v0[i], v1[i], v1S[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 3*q.Value())

		for i := 0; i < N; i++ {
			vOut[i] %= q.Value()
			vOutCheck[i] = num.Sub(vOutInit[i], num.SMul(v0[i], v1[i], v1S[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("Reduce", func(t *testing.T) {
		for i := 0; i < N; i++ {
			v0[i] = rSrc.Sample()
		}

		vec.ReduceTo(vOut, v0, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] % q.Value()
		}
		assert.Equal(t, vOutCheck, vOut)
	})
}

func BenchmarkOps(b *testing.B) {
	q := num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)

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
					vec.AddTo(vOut, v0, v1, q)
				}
			})

			b.Run("AddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.AddLazyTo(vOut, v0, v1)
				}
			})

			b.Run("ScalarAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarAddTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarAddLazyTo(vOut, v0, v1[0])
				}
			})

			b.Run("Sub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SubTo(vOut, v0, v1, q)
				}
			})

			b.Run("SubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SubLazyTo(vOut, v0, v1)
				}
			})

			b.Run("ScalarSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarSubTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarSubLazyTo(vOut, v0, v1[0])
				}
			})

			b.Run("Neg", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.NegTo(vOut, v0, q)
				}
			})

			b.Run("ScalarMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMulTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMulAddTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMulSubTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMulLazyTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMulAddLazyTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMulSubLazyTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMMulTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMMulAddTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMMulSubTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMMulLazyTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMMulAddLazyTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("ScalarMMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ScalarMMulSubLazyTo(vOut, v0, v1[0], q)
				}
			})

			b.Run("Mul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MulTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MulAddTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MulSubTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MulLazyTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MulAddLazyTo(vOut, v0, v1, q)
				}
			})

			b.Run("MulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MulSubLazyTo(vOut, v0, v1, q)
				}
			})

			b.Run("MForm", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MFormTo(vOut, v0, q)
				}
			})

			b.Run("InvMForm", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.InvMFormTo(vOut, v0, q)
				}
			})

			b.Run("MMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MMulTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MMulAddTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MMulSubTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MMulLazyTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MMulAddLazyTo(vOut, v0, v1, q)
				}
			})

			b.Run("MMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.MMulSubLazyTo(vOut, v0, v1, q)
				}
			})

			b.Run("SForm", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SFormTo(vOut, v0, q)
				}
			})

			b.Run("SMul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SMulTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SMulAddTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SMulSubTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SMulLazyTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulAddLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SMulAddLazyTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("SMulSubLazy", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.SMulSubLazyTo(vOut, v0, v1, v1S, q)
				}
			})

			b.Run("Reduce", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					vec.ReduceTo(vOut, v0, q)
				}
			})
		})
	}
}
