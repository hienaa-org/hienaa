package poly_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/poly"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = csprng.NewUniformSamplerWithSeed(nil)
)

func TestTransform(t *testing.T) {

	t.Run("CyclotomicPow2", func(t *testing.T) {
		testCyclotomicPow2(t, mod.NewModulus(65537))
	})
}

func BenchmarkTransform(b *testing.B) {
	b.Run("CyclotomicPow2", func(b *testing.B) {
		benchmarkCyclotomicPow2(b, mod.NewModulus(0x80000000080001))
	})
}

func testCyclotomicPow2(t *testing.T, q *mod.Modulus) {
	N := 32
	ringParams := poly.NewCyclotomicParameters(N << 1)
	ntt := poly.NewTransformer(ringParams, q)

	p := make([]uint64, N)
	for i := 0; i < N; i++ {
		p[i] = rSrc.SampleN(q.Value())
	}

	t.Run("NTT", func(t *testing.T) {
		pCopy := make([]uint64, N)
		copy(pCopy, p)

		ntt.NTTInPlace(pCopy)
		ntt.InvNTTInPlace(pCopy)

		assert.Equal(t, p, pCopy)
	})

	t.Run("Mul", func(t *testing.T) {
		pMul := make([]uint64, N)
		for i := 0; i < N; i++ {
			for j := 0; j < N; j++ {
				if i+j < N {
					pMul[(i+j)%N] = mod.Add(pMul[(i+j)%N], mod.Mul(p[i], p[j], q), q)
				} else {
					pMul[(i+j)%N] = mod.Sub(pMul[(i+j)%N], mod.Mul(p[i], p[j], q), q)
				}
			}
		}

		ntt.NTTInPlace(p)
		mod.MMulVecTo(p, p, q, p)
		ntt.InvNTTInPlace(p)

		assert.Equal(t, pMul, p)
	})
}

func benchmarkCyclotomicPow2(b *testing.B, q *mod.Modulus) {
	for _, logN := range []int{12, 13, 14, 15, 16, 17} {
		N := 1 << logN
		ringParams := poly.NewCyclotomicParameters(N << 1)
		ntt := poly.NewTransformer(ringParams, q)

		p := make([]uint64, ringParams.Degree())
		for i := 0; i < ringParams.Degree(); i++ {
			p[i] = rSrc.SampleN(q.Value())
		}

		b.Run(fmt.Sprintf("logN=%d", logN), func(b *testing.B) {
			b.Run("NTT", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ntt.NTTInPlace(p)
				}
			})
			b.Run("InvNTT", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					ntt.InvNTTInPlace(p)
				}
			})
		})
	}
}
