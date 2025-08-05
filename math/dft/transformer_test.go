package dft_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{12, 13, 14, 15, 16, 17}
)

func randPoly(ringParams dft.RingParameters, q *num.Modulus) []uint64 {
	p := make([]uint64, ringParams.Rank())
	for i := 0; i < ringParams.Rank(); i++ {
		p[i] = rSrc.SampleN(q.Value())
	}
	return p
}

func cyclotomicPow2Mul(p0, p1 []uint64, q *num.Modulus) []uint64 {
	N := len(p0)

	pOut := make([]uint64, N)
	for i := range p0 {
		for j := range p1 {
			if i+j < N {
				pOut[(i+j)%N] = num.Add(pOut[(i+j)%N], num.Mul(p0[i], p1[j], q), q)
			} else {
				pOut[(i+j)%N] = num.Sub(pOut[(i+j)%N], num.Mul(p0[i], p1[j], q), q)
			}
		}
	}
	return pOut
}

func cyclicMul(p0, p1 []uint64, q *num.Modulus) []uint64 {
	N := len(p0)

	pOut := make([]uint64, N)
	for i := range p0 {
		for j := range p1 {
			pOut[(i+j)%N] = num.Add(pOut[(i+j)%N], num.Mul(p0[i], p1[j], q), q)
		}
	}
	return pOut
}

func reduce(p0, p1 []uint64, q *num.Modulus) []uint64 {
	quo := make([]uint64, len(p0)-len(p1)+1)
	rem := make([]uint64, len(p0))
	copy(rem, p0)

	for i := 0; i <= len(p0)-len(p1); i++ {
		if rem[len(rem)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(rem[len(rem)-i-1], num.Inv(p1[len(p1)-1], q), q)

			for j := 0; j < len(p1); j++ {
				rem[len(rem)-i-j-1] = num.Sub(rem[len(rem)-i-j-1], num.Mul(p1[len(p1)-j-1], quo[len(quo)-i-1], q), q)
			}
		}
	}

	return rem
}

func TestCyclotomicNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 12
		ringParams := dft.NewCyclotomicParameters(N << 1)
		qs := dft.FindNearestNTTPrimes(ringParams, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())
		ntt := dft.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q)
		p1 := randPoly(ringParams, q)

		p0NTT := make([]uint64, N)
		copy(p0NTT, p0)
		ntt.ForwardInPlace(p0NTT)

		p1NTT := make([]uint64, N)
		copy(p1NTT, p1)
		ntt.ForwardInPlace(p1NTT)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseInPlace(pOut)

		assert.Equal(t, cyclotomicPow2Mul(p0, p1, q), pOut)
	})

	t.Run("type=Any", func(t *testing.T) {
		M := int(num.NextPrime(1<<12, 1))
		ringParams := dft.NewCyclotomicParameters(M)
		N := ringParams.Rank()
		q := dft.FindNearestNTTPrimes(ringParams, 45, 1)[0]
		ntt := dft.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q)
		p1 := randPoly(ringParams, q)

		p0NTT := make([]uint64, N)
		copy(p0NTT, p0)
		ntt.ForwardInPlace(p0NTT)

		p1NTT := make([]uint64, N)
		copy(p1NTT, p1)
		ntt.ForwardInPlace(p1NTT)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseInPlace(pOut)

		p0Ref := append(p0, make([]uint64, (2*N-1)-len(p0))...)
		p1Ref := append(p1, make([]uint64, (2*N-1)-len(p1))...)
		pOutRef := cyclicMul(p0Ref, p1Ref, q)

		cycloSigned := dftops.CyclotomicPolynomial(ringParams.CycloOrder())
		cyclo := make([]uint64, len(cycloSigned))
		for i := range cycloSigned {
			if cycloSigned[i] < 0 {
				cyclo[i] = uint64(int(q.Value()) + cycloSigned[i])
			} else {
				cyclo[i] = uint64(cycloSigned[i])
			}
		}

		assert.Equal(t, reduce(pOutRef, cyclo, q)[:ringParams.Rank()], pOut)
	})
}

func TestCyclicNTT(t *testing.T) {
	t.Run("type=Pow235", func(t *testing.T) {
		N := int(num.NextProdPower(rSrc.SampleN(1<<12), []uint64{2, 3, 5}))
		ringParams := dft.NewCyclicParameters(N)
		qs := dft.FindNearestNTTPrimes(ringParams, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())
		ntt := dft.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q)
		p1 := randPoly(ringParams, q)

		p0NTT := make([]uint64, N)
		copy(p0NTT, p0)
		ntt.ForwardInPlace(p0NTT)

		p1NTT := make([]uint64, N)
		copy(p1NTT, p1)
		ntt.ForwardInPlace(p1NTT)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseInPlace(pOut)

		assert.Equal(t, cyclicMul(p0, p1, q), pOut)
	})

	t.Run("type=Bluestein", func(t *testing.T) {
		var N int
		for {
			N = int(rSrc.SampleN(1 << 12))
			if !num.IsProdPowerOf(uint64(N), []uint64{2, 3, 5}) {
				break
			}
		}
		ringParams := dft.NewCyclicParameters(N)
		qs := dft.FindNearestNTTPrimes(ringParams, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())
		ntt := dft.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q)
		p1 := randPoly(ringParams, q)

		p0NTT := make([]uint64, N)
		copy(p0NTT, p0)
		ntt.ForwardInPlace(p0NTT)

		p1NTT := make([]uint64, N)
		copy(p1NTT, p1)
		ntt.ForwardInPlace(p1NTT)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseInPlace(pOut)

		assert.Equal(t, cyclicMul(p0, p1, q), pOut)
	})
}

func BenchmarkCyclotomicNTT(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			ringParams := dft.NewCyclotomicParameters(N << 1)
			q := dft.FindNextNTTPrimes(ringParams, 60, 1)[0]
			ntt := dft.NewTransformer(ringParams, q)

			p := randPoly(ringParams, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardInPlace(p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseInPlace(p)
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
			ringParams := dft.NewCyclicParameters(N)
			q := dft.FindNextNTTPrimes(ringParams, 60, 1)[0]
			ntt := dft.NewTransformer(ringParams, q)

			p := randPoly(ringParams, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardInPlace(p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseInPlace(p)
					}
				})
			})
		}
	})

	b.Run("type=Bluestein", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := (1 << logN) + 1
			ringParams := dft.NewCyclicParameters(N)
			q := dft.FindNextNTTPrimes(ringParams, 60, 1)[0]
			ntt := dft.NewTransformer(ringParams, q)

			p := randPoly(ringParams, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardInPlace(p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseInPlace(p)
					}
				})
			})
		}
	})
}
