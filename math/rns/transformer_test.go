package rns_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/rns"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{12, 13, 14, 15, 16, 17}
)

func randPoly(ringParams rns.RingParameters, q *mod.Modulus) *rns.Poly {
	p := rns.NewPoly(ringParams.Degree(), 1)
	for i := 0; i < ringParams.Degree(); i++ {
		p.Coeffs[0][i] = rSrc.SampleN(q.Value())
	}
	return p
}

func cyclotomicPow2Mul(p0, p1 []uint64, q *mod.Modulus) []uint64 {
	N := len(p0)

	pOut := make([]uint64, N)
	for i := range p0 {
		for j := range p1 {
			if i+j < N {
				pOut[(i+j)%N] = mod.Add(pOut[(i+j)%N], mod.Mul(p0[i], p1[j], q), q)
			} else {
				pOut[(i+j)%N] = mod.Sub(pOut[(i+j)%N], mod.Mul(p0[i], p1[j], q), q)
			}
		}
	}
	return pOut
}

func cyclicMul(p0, p1 []uint64, q *mod.Modulus) []uint64 {
	N := len(p0)

	pOut := make([]uint64, N)
	for i := range p0 {
		for j := range p1 {
			pOut[(i+j)%N] = mod.Add(pOut[(i+j)%N], mod.Mul(p0[i], p1[j], q), q)
		}
	}
	return pOut
}

func TestCyclotomicNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 12
		ringParams := rns.NewCyclotomicParameters(N)
		q := rns.FindNextNTTPrimes(ringParams, 45, 1)
		ntt := rns.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q[0])
		p1 := randPoly(ringParams, q[0])

		p0NTT := ntt.NTT(p0)
		p1NTT := ntt.NTT(p1)

		pOut := rns.NewNTTPoly(ringParams.Degree(), 1)
		mod.MMulVecTo(pOut.Coeffs[0], p0NTT.Coeffs[0], p1NTT.Coeffs[0], q[0])
		ntt.InvNTTTo(pOut, pOut)

		assert.Equal(t, cyclotomicPow2Mul(p0.Coeffs[0], p1.Coeffs[0], q[0]), pOut.Coeffs[0])
	})
}

func TestCyclicNTT(t *testing.T) {
	t.Run("type=Pow235", func(t *testing.T) {
		N := int(num.NextProdPower(rSrc.SampleN(1<<12), []uint64{2, 3, 5}))
		ringParams := rns.NewCyclicParameters(N)
		q := rns.FindNextNTTPrimes(ringParams, 45, 1)
		ntt := rns.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q[0])
		p1 := randPoly(ringParams, q[0])

		p0NTT := ntt.NTT(p0)
		p1NTT := ntt.NTT(p1)

		pOut := rns.NewNTTPoly(ringParams.Degree(), 1)
		mod.MMulVecTo(pOut.Coeffs[0], p0NTT.Coeffs[0], p1NTT.Coeffs[0], q[0])
		ntt.InvNTTTo(pOut, pOut)

		assert.Equal(t, cyclicMul(p0.Coeffs[0], p1.Coeffs[0], q[0]), pOut.Coeffs[0])
	})

	t.Run("type=Bluestein", func(t *testing.T) {
		var N int
		for {
			N = int(rSrc.SampleN(1 << 12))
			if !num.IsProdPowerOf(uint64(N), []uint64{2, 3, 5}) {
				break
			}
		}
		ringParams := rns.NewCyclicParameters(N)
		q := rns.FindNextNTTPrimes(ringParams, 45, 1)
		ntt := rns.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q[0])
		p1 := randPoly(ringParams, q[0])

		p0NTT := ntt.NTT(p0)
		p1NTT := ntt.NTT(p1)

		pOut := rns.NewNTTPoly(ringParams.Degree(), 1)
		mod.MMulVecTo(pOut.Coeffs[0], p0NTT.Coeffs[0], p1NTT.Coeffs[0], q[0])
		ntt.InvNTTTo(pOut, pOut)

		assert.Equal(t, cyclicMul(p0.Coeffs[0], p1.Coeffs[0], q[0]), pOut.Coeffs[0])
	})
}

func BenchmarkCyclotomicNTT(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			ringParams := rns.NewCyclotomicParameters(N)
			q := rns.FindNextNTTPrimes(ringParams, 60, 1)
			ntt := rns.NewTransformer(ringParams, q)

			p := randPoly(ringParams, q[0])
			pOut := rns.NewNTTPoly(ringParams.Degree(), len(q))

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.NTTTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InvNTTTo(p, pOut)
					}
				})
			})
		}
	})
}

func BenchmarkCyclicNTT(b *testing.B) {
	b.Run("type=Pow235", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := int(num.NextProdPower(uint64(1<<logN)+1, []uint64{2, 3, 5}))
			ringParams := rns.NewCyclicParameters(N)
			q := rns.FindNextNTTPrimes(ringParams, 60, 1)
			ntt := rns.NewTransformer(ringParams, q)

			p := randPoly(ringParams, q[0])
			pOut := rns.NewNTTPoly(ringParams.Degree(), len(q))

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.NTTTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InvNTTTo(p, pOut)
					}
				})
			})
		}
	})

	b.Run("type=Bluestein", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := (1 << logN) + 1
			ringParams := rns.NewCyclicParameters(N)
			q := rns.FindNextNTTPrimes(ringParams, 60, 1)
			ntt := rns.NewTransformer(ringParams, q)

			p := randPoly(ringParams, q[0])
			pOut := rns.NewNTTPoly(ringParams.Degree(), len(q))

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.NTTTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InvNTTTo(p, pOut)
					}
				})
			})
		}
	})
}
