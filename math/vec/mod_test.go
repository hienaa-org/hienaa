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
	benchLogN = []int{10, 11, 12, 13, 14, 15, 16}
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

	v1M := vec.ToMulForm(v1, q)

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

	t.Run("AddWord", func(t *testing.T) {
		vec.AddTo(vOut, v0, v1, nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), 2*q.Value())
	})

	t.Run("AddScalar", func(t *testing.T) {
		vec.AddScalarTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] + v1[0]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] -= q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("AddScalarWord", func(t *testing.T) {
		vec.AddScalarTo(vOut, v0, v1[0], nil)
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

	t.Run("SubWord", func(t *testing.T) {
		vec.SubTo(vOut, v0, v1, nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("SubScalar", func(t *testing.T) {
		vec.SubScalarTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] - v1[0]
			if vOutCheck[i] >= q.Value() {
				vOutCheck[i] += q.Value()
			}
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("SubScalarWord", func(t *testing.T) {
		vec.SubScalarTo(vOut, v0, v1[0], nil)
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

	t.Run("NegWord", func(t *testing.T) {
		vec.NegTo(vOut, v0, nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] = -v0[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulScalar", func(t *testing.T) {
		vec.MulScalarTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MulAddScalar", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulAddScalarTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MulSubScalar", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulSubScalarTo(vOut, v0, v1[0], q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("MulScalarWord", func(t *testing.T) {
		vec.MulScalarTo(vOut, v0, v1[0], nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] * v1[0]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulAddScalarWord", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulAddScalarTo(vOut, v0, v1[0], nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] += v0[i] * v1[0]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulSubScalarWord", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulSubScalarTo(vOut, v0, v1[0], nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] -= v0[i] * v1[0]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	v1cM := num.MulForm{Float: v1M.Float[0], SForm: v1M.SForm[0]}
	t.Run("FMulScalar", func(t *testing.T) {
		vec.FMulScalarTo(vOut, v0, v1[0], v1cM, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Mul(v0[i], v1[0], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("FMulAddScalar", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.FMulAddScalarTo(vOut, v0, v1[0], v1cM, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("FMulSubScalar", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.FMulSubScalarTo(vOut, v0, v1[0], v1cM, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.Mul(v0[i], v1[0], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
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

	t.Run("MulWord", func(t *testing.T) {
		vec.MulTo(vOut, v0, v1, nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] * v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulAddWord", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulAddTo(vOut, v0, v1, nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] += v0[i] * v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("MulSubWord", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.MulSubTo(vOut, v0, v1, nil)
		for i := 0; i < N; i++ {
			vOutCheck[i] -= v0[i] * v1[i]
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("FMul", func(t *testing.T) {
		vec.FMulTo(vOut, v0, v1, v1M, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Mul(v0[i], v1[i], q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("FMulAdd", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.FMulAddTo(vOut, v0, v1, v1M, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Add(vOutCheck[i], num.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
	})

	t.Run("FMulSub", func(t *testing.T) {
		copy(vOut, vOutInit)
		copy(vOutCheck, vOutInit)

		vec.FMulSubTo(vOut, v0, v1, v1M, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = num.Sub(vOutCheck[i], num.Mul(v0[i], v1[i], q), q)
		}
		assert.Equal(t, vOutCheck, vOut)

		assert.Less(t, vec.Max(vOut), q.Value())
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

	t.Run("Reduce2Q", func(t *testing.T) {
		for i := 0; i < N; i++ {
			v0[i] = rSrc.Sample() % (2 * q.Value())
		}

		vec.Reduce2QTo(vOut, v0, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] % q.Value()
		}
		assert.Equal(t, vOutCheck, vOut)
	})

	t.Run("Reduce4Q", func(t *testing.T) {
		for i := 0; i < N; i++ {
			v0[i] = rSrc.Sample() % (4 * q.Value())
		}

		vec.Reduce4QTo(vOut, v0, q)
		for i := 0; i < N; i++ {
			vOutCheck[i] = v0[i] % q.Value()
		}
		assert.Equal(t, vOutCheck, vOut)
	})
}

func benchmarkOps(b *testing.B, logN int) {
	q := num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)

	N := 1 << logN
	v0 := make([]uint64, N)
	v1 := make([]uint64, N)
	vOut := make([]uint64, N)

	for i := 0; i < N; i++ {
		v0[i] = rSrc.SampleN(q.Value())
		v1[i] = rSrc.SampleN(q.Value())
	}

	v1M := vec.ToMulForm(v1, q)

	b.Run("Add", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.AddTo(vOut, v0, v1, q)
		}
	})

	b.Run("AddWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.AddTo(vOut, v0, v1, nil)
		}
	})

	b.Run("AddScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.AddScalarTo(vOut, v0, v1[0], q)
		}
	})

	b.Run("AddScalarWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.AddScalarTo(vOut, v0, v1[0], nil)
		}
	})

	b.Run("Sub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.SubTo(vOut, v0, v1, q)
		}
	})

	b.Run("SubWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.SubTo(vOut, v0, v1, nil)
		}
	})

	b.Run("SubScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.SubScalarTo(vOut, v0, v1[0], q)
		}
	})

	b.Run("SubScalarWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.SubScalarTo(vOut, v0, v1[0], nil)
		}
	})

	b.Run("Neg", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.NegTo(vOut, v0, q)
		}
	})

	b.Run("MulScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulScalarTo(vOut, v0, v1[0], q)
		}
	})

	b.Run("MulAddScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulAddScalarTo(vOut, v0, v1[0], q)
		}
	})

	b.Run("MulSubScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulSubScalarTo(vOut, v0, v1[0], q)
		}
	})

	b.Run("MulScalarWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulScalarTo(vOut, v0, v1[0], nil)
		}
	})

	b.Run("MulAddScalarWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulAddScalarTo(vOut, v0, v1[0], nil)
		}
	})

	b.Run("MulSubScalarWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulSubScalarTo(vOut, v0, v1[0], nil)
		}
	})

	v1cM := num.ToMulForm(v1[0], q)
	b.Run("FMulScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.FMulScalarTo(vOut, v0, v1[0], v1cM, q)
		}
	})

	b.Run("FMulAddScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.FMulAddScalarTo(vOut, v0, v1[0], v1cM, q)
		}
	})

	b.Run("FMulSubScalar", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.FMulSubScalarTo(vOut, v0, v1[0], v1cM, q)
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

	b.Run("MulWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulTo(vOut, v0, v1, nil)
		}
	})

	b.Run("MulAddWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulAddTo(vOut, v0, v1, nil)
		}
	})

	b.Run("MulSubWord", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.MulSubTo(vOut, v0, v1, nil)
		}
	})

	b.Run("FMul", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.FMulTo(vOut, v0, v1, v1M, q)
		}
	})

	b.Run("FMulAdd", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.FMulAddTo(vOut, v0, v1, v1M, q)
		}
	})

	b.Run("FMulSub", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.FMulSubTo(vOut, v0, v1, v1M, q)
		}
	})

	b.Run("Reduce", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.ReduceTo(vOut, v0, q)
		}
	})

	b.Run("Reduce2Q", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.Reduce2QTo(vOut, v0, q)
		}
	})

	b.Run("Reduce4Q", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			vec.Reduce4QTo(vOut, v0, q)
		}
	})
}

func BenchmarkOps(b *testing.B) {
	for _, logN := range benchLogN {
		b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
			benchmarkOps(b, logN)
		})
	}
}
