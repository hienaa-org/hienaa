package crt_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	reducerBenchLogN = []int{12, 13, 14}
)

func randTernaryPoly(rank int) []int64 {
	p := make([]int64, rank)
	for i := 0; i < rank-1; i++ {
		p[i] = int64(rSrc.SampleN(3)) - 1
	}
	p[rank-1] = 1
	return p
}

func randPoly(rank int, q []*num.Modulus) *crt.Element {
	p := crt.NewPoly(rank, len(q))
	for i := 0; i < rank; i++ {
		for j := range q {
			p.Coeffs[j][i] = rSrc.SampleN(q[j].Value())
		}
	}

	return p
}

func TestReducer(t *testing.T) {
	t.Run("type=Cyclotomic", func(t *testing.T) {
		sqrtN := int(math.Sqrt(math.Exp2(10)))
		m0 := num.MustNextPrime(sqrtN, 1)
		m1 := num.MustNextPrime(m0, 2)
		M := m0 * m1
		rP := dft.NewCyclotomicParameters(M)
		N := rP.Rank()

		degNext := num.NextProdPower(N, []int{2})
		diffDegNext := num.NextProdPower(2*(M-N)-1, []int{2})
		ambParams := dft.NewCyclicParameters(max(degNext, diffDegNext))

		q := dft.MustFindPrevNTTPrimes(ambParams, 40, 1)
		q = append(q, num.NewModulus(num.MustNextPrime(q[0].Value(), 2)))

		cycloReducer := crt.NewCyclotomicReducer(rP, q)
		longDivReducer := crt.NewLongDivReducer(M, q, dft.CyclotomicPolynomial(M))

		p := randPoly(M, q)
		pOut := cycloReducer.Reduce(p)
		pOutRef := longDivReducer.Reduce(p)

		assert.Equal(t, pOutRef, pOut)
	})

	t.Run("type=Any", func(t *testing.T) {
		N := int(rSrc.SampleN(1<<10)) + 4
		maxRank := N + int(rSrc.SampleN(1<<10)) + 1

		degNext := num.NextProdPower(N, []int{2})
		diffDegNext := num.NextProdPower(2*(maxRank-N)-1, []int{2})
		ambParams := dft.NewCyclicParameters(max(degNext, diffDegNext))

		q := dft.MustFindPrevNTTPrimes(ambParams, 40, 1)
		q = append(q, num.NewModulus(num.MustNextPrime(q[0].Value(), 2)))

		modPoly := randTernaryPoly(N + 1)
		reducer := crt.NewReducer(maxRank, q, modPoly)
		longDivReducer := crt.NewLongDivReducer(maxRank, q, modPoly)

		p := randPoly(maxRank, q)
		pOut := reducer.Reduce(p)
		pOutRef := longDivReducer.Reduce(p)

		assert.Equal(t, pOutRef, pOut)
	})
}

func BenchmarkReducer(b *testing.B) {
	b.Run("type=Cyclotomic", func(b *testing.B) {
		for _, logN := range reducerBenchLogN {
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				sqrtN := int(math.Sqrt(math.Exp2(float64(logN))))
				m0 := num.MustNextPrime(sqrtN, 1)
				m1 := num.MustNextPrime(m0, 2)
				M := m0 * m1
				rP := dft.NewCyclotomicParameters(M)
				N := rP.Rank()

				degNext := num.NextProdPower(N, []int{2})
				diffDegNext := num.NextProdPower(2*(M-N)-1, []int{2})
				ambParams := dft.NewCyclicParameters(max(degNext, diffDegNext))

				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.MustFindPrevNTTPrimes(ambParams, num.MaxModulusBits, 1)

					reducer := crt.NewCyclotomicReducer(rP, q)

					p := randPoly(M, q)
					pOut := crt.NewPoly(N, len(q))

					b.ResetTimer()
					for i := 0; i < b.N; i++ {
						reducer.ReduceTo(pOut, p)
					}
				})

				b.Run("Mod=Any", func(b *testing.B) {
					q := []*num.Modulus{num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)}

					reducer := crt.NewCyclotomicReducer(rP, q)

					p := randPoly(M, q)
					pOut := crt.NewPoly(N, len(q))

					b.ResetTimer()
					for i := 0; i < b.N; i++ {
						reducer.ReduceTo(pOut, p)
					}
				})
			})
		}
	})

	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range reducerBenchLogN {
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				N := 1 << logN
				maxRank := 2 * N

				modPoly := randTernaryPoly(N + 1)

				degNext := num.NextProdPower(N, []int{2})
				diffDegNext := num.NextProdPower(2*(maxRank-N)-1, []int{2})
				ambParams := dft.NewCyclicParameters(max(degNext, diffDegNext))

				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.MustFindPrevNTTPrimes(ambParams, num.MaxModulusBits, 1)

					reducer := crt.NewReducer(maxRank, q, modPoly)

					p := randPoly(maxRank, q)
					pOut := crt.NewPoly(N, len(q))

					b.ResetTimer()
					for i := 0; i < b.N; i++ {
						reducer.ReduceTo(pOut, p)
					}
				})

				b.Run("Mod=Any", func(b *testing.B) {
					q := []*num.Modulus{num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)}

					reducer := crt.NewReducer(maxRank, q, modPoly)

					p := randPoly(maxRank, q)
					pOut := crt.NewPoly(N, 1)

					b.ResetTimer()
					for i := 0; i < b.N; i++ {
						reducer.ReduceTo(pOut, p)
					}
				})
			})
		}
	})
}
