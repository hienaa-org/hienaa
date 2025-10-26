package dft_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{12, 13, 14, 15, 16, 17}
)

func randPoly(params dft.RingParameters, q *num.Modulus) []uint64 {
	p := make([]uint64, params.Rank())
	for i := 0; i < params.Rank(); i++ {
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

	lcInv := num.Inv(p1[len(p1)-1], q)
	for i := 0; i <= len(p0)-len(p1); i++ {
		if rem[len(rem)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(rem[len(rem)-i-1], lcInv, q)
			vec.ScalarMulSubTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], q)
		}
	}

	return rem[:len(p1)-1]
}

func TestCyclotomicNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10
		rP := dft.NewCyclotomicParameters(2 * N)

		qs := dft.FindNearestNTTPrimes(rP, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())

		ntt := dft.NewTransformer(rP, q)

		p0 := randPoly(rP, q)
		p1 := randPoly(rP, q)

		p0NTT := make([]uint64, N)
		ntt.ForwardTo(p0NTT, p0)

		p1NTT := make([]uint64, N)
		ntt.ForwardTo(p1NTT, p1)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseTo(pOut, pOut)

		assert.Equal(t, cyclotomicPow2Mul(p0, p1, q), pOut)
	})

	t.Run("type=Any", func(t *testing.T) {
		sqrtN := int(math.Sqrt(math.Exp2(10)))
		m0 := num.NextPrime(sqrtN, 1)
		m1 := num.NextPrime(m0, 2)
		M := m0 * m1
		rP := dft.NewCyclotomicParameters(M)
		N := rP.Rank()

		qs := dft.FindNearestNTTPrimes(rP, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())

		ntt := dft.NewTransformer(rP, q)

		p0 := randPoly(rP, q)
		p1 := randPoly(rP, q)

		p0NTT := make([]uint64, N)
		ntt.ForwardTo(p0NTT, p0)

		p1NTT := make([]uint64, N)
		ntt.ForwardTo(p1NTT, p1)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseTo(pOut, pOut)

		p0Ref := append(p0, make([]uint64, (2*N-1)-len(p0))...)
		p1Ref := append(p1, make([]uint64, (2*N-1)-len(p1))...)
		pOutRef := cyclicMul(p0Ref, p1Ref, q)

		cyclo := vec.Reduce(dft.CyclotomicPolynomial(rP.CycloOrder()), q)

		assert.Equal(t, reduce(pOutRef, cyclo, q), pOut)
	})
}

func TestCyclicNTT(t *testing.T) {
	t.Run("type=Pow235", func(t *testing.T) {
		N := num.NextProdPower(int(rSrc.SampleN(1<<10)), []int{2, 3, 5})
		rP := dft.NewCyclicParameters(N)

		qs := dft.FindNearestNTTPrimes(rP, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())

		ntt := dft.NewTransformer(rP, q)

		p0 := randPoly(rP, q)
		p1 := randPoly(rP, q)

		p0NTT := make([]uint64, N)
		ntt.ForwardTo(p0NTT, p0)

		p1NTT := make([]uint64, N)
		ntt.ForwardTo(p1NTT, p1)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseTo(pOut, pOut)

		assert.Equal(t, cyclicMul(p0, p1, q), pOut)
	})

	t.Run("type=Any", func(t *testing.T) {
		var N int
		for {
			N = int(rSrc.SampleN(1 << 10))
			if !num.IsProdPowerOf(N, []int{2, 3, 5}) {
				break
			}
		}
		rP := dft.NewCyclicParameters(N)

		qs := dft.FindNearestNTTPrimes(rP, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())

		ntt := dft.NewTransformer(rP, q)

		p0 := randPoly(rP, q)
		p1 := randPoly(rP, q)

		p0NTT := make([]uint64, N)
		ntt.ForwardTo(p0NTT, p0)

		p1NTT := make([]uint64, N)
		ntt.ForwardTo(p1NTT, p1)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseTo(pOut, pOut)

		assert.Equal(t, cyclicMul(p0, p1, q), pOut)
	})
}

func TestAutFixedNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10
		rP := dft.NewAutFixedParameters(4*N, N)

		qs := dft.FindNearestNTTPrimes(rP, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())

		ntt := dft.NewTransformer(rP, q)

		p0 := randPoly(rP, q)
		p1 := randPoly(rP, q)

		p0Ref := make([]uint64, 2*N)
		p1Ref := make([]uint64, 2*N)
		copy(p0Ref[:N], p0)
		copy(p1Ref[:N], p1)

		p0Ref[N] = 0
		p1Ref[N] = 0
		for i := 1; i < N; i++ {
			p0Ref[i+N] = num.Neg(p0[N-i], q)
			p1Ref[i+N] = num.Neg(p1[N-i], q)
		}

		p0NTT := make([]uint64, N)
		ntt.ForwardTo(p0NTT, p0)

		p1NTT := make([]uint64, N)
		ntt.ForwardTo(p1NTT, p1)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseTo(pOut, pOut)

		assert.Equal(t, cyclotomicPow2Mul(p0Ref, p1Ref, q)[:N], pOut)
	})

	t.Run("type=Prime", func(t *testing.T) {
		M := num.NextPrime(int(rSrc.SampleN(1<<10)), 1)
		primes, exps := num.Factor(M - 1)
		fold := 1
		for i := range primes {
			e := int(rSrc.SampleN(uint64(exps[i])))
			for j := 0; j < e; j++ {
				fold *= primes[i]
			}
		}
		N := (M - 1) / fold
		rP := dft.NewAutFixedParameters(M, N)

		qs := dft.FindNearestNTTPrimes(rP, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())

		ntt := dft.NewTransformer(rP, q)

		p0 := randPoly(rP, q)
		p1 := randPoly(rP, q)

		p0Ref := make([]uint64, M)
		p1Ref := make([]uint64, M)

		cycloOrdMod := num.NewModulus(M)
		root := num.Generators(cycloOrdMod)[0]
		idx := uint64(1)
		for i := 0; i < fold; i++ {
			for j := 0; j < N; j++ {
				p0Ref[idx], p1Ref[idx] = p0[j], p1[j]
				idx = num.Mul(idx, root, cycloOrdMod)
			}
		}

		p0NTT := make([]uint64, N)
		ntt.ForwardTo(p0NTT, p0)

		p1NTT := make([]uint64, N)
		ntt.ForwardTo(p1NTT, p1)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseTo(pOut, pOut)

		pOutRef := make([]uint64, N)
		pOutRefLong := cyclicMul(p0Ref, p1Ref, q)
		for i := 1; i < M; i++ {
			pOutRefLong[i] = num.Sub(pOutRefLong[i], pOutRefLong[0], q)
		}
		pOutRefLong[0] = 0

		idx = 1
		for i := 0; i < N; i++ {
			pOutRef[i] = pOutRefLong[idx]
			idx = num.Mul(idx, root, cycloOrdMod)
		}

		assert.Equal(t, pOutRef, pOut)
	})
}

func BenchmarkCyclotomicNTT(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(2 * N)

			q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)[0]

			ntt := dft.NewTransformer(rP, q)

			p := randPoly(rP, q)
			pOut := randPoly(rP, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseTo(pOut, p)
					}
				})
			})
		}
	})

	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range benchLogN {
			sqrtN := int(math.Sqrt(math.Exp2(float64(logN))))
			m0 := num.NextPrime(sqrtN, 1)
			m1 := num.NextPrime(m0, 2)
			M := m0 * m1
			rP := dft.NewCyclotomicParameters(M)

			q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)[0]

			ntt := dft.NewTransformer(rP, q)

			p := randPoly(rP, q)
			pOut := randPoly(rP, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseTo(pOut, p)
					}
				})
			})
		}
	})
}

func BenchmarkCyclicNTT(b *testing.B) {
	b.Run("type=Pow235", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := num.NextProdPower((1<<logN)+1, []int{2, 3, 5})
			rP := dft.NewCyclicParameters(N)

			q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)[0]

			ntt := dft.NewTransformer(rP, q)

			p := randPoly(rP, q)
			pOut := randPoly(rP, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseTo(pOut, p)
					}
				})
			})
		}
	})

	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := (1 << logN) + 1
			rP := dft.NewCyclicParameters(N)

			q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)[0]

			ntt := dft.NewTransformer(rP, q)

			p := randPoly(rP, q)
			pOut := randPoly(rP, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseTo(pOut, p)
					}
				})
			})
		}
	})
}

func BenchmarkAutFixedNTT(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewAutFixedParameters(4*N, N)

			q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)[0]

			ntt := dft.NewTransformer(rP, q)

			p := randPoly(rP, q)
			pOut := randPoly(rP, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseTo(pOut, p)
					}
				})
			})
		}
	})

	b.Run("type=Prime", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			M := num.NextPrime(1, N)
			rP := dft.NewAutFixedParameters(M, N)

			q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)[0]

			ntt := dft.NewTransformer(rP, q)

			p := randPoly(rP, q)
			pOut := randPoly(rP, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("NTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.ForwardTo(pOut, p)
					}
				})
				b.Run("InvNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						ntt.InverseTo(pOut, p)
					}
				})
			})
		}
	})
}
