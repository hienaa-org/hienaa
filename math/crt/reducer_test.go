package crt_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
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

func randPoly(rank int, q []*num.Modulus) *crt.Poly {
	p := crt.NewPoly(rank, len(q))
	for i := 0; i < rank; i++ {
		for j := range q {
			p.Coeffs[j][i] = rSrc.SampleN(q[j].Value())
		}
	}

	return p
}

func reduce(p0, p1 []uint64, q *num.Modulus) []uint64 {
	quo := make([]uint64, len(p0)-len(p1)+1)
	rem := make([]uint64, len(p0))
	copy(rem, p0)

	for i := 0; i <= len(p0)-len(p1); i++ {
		if rem[len(rem)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(rem[len(rem)-i-1], num.Inv(p1[len(p1)-1], q), q)
			vec.ScalarMulSubTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], q)
		}
	}

	return rem[:len(p1)-1]
}

func TestReducer(t *testing.T) {
	N := int(rSrc.SampleN(1<<10)) + 4
	maxRank := N + int(rSrc.SampleN(1<<10)) + 1
	q0 := num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)

	degNext := num.NextProdPower(N, []int{2})
	diffDegNext := num.NextProdPower(2*(maxRank-N)-1, []int{2})
	redParams := dft.NewCyclicParameters(max(degNext, diffDegNext))
	q1 := dft.FindNearestNTTPrimes(redParams, 60, 1)[0]

	q := []*num.Modulus{q0, q1}

	modPoly := randTernaryPoly(N + 1)
	reducer := crt.NewReducer(maxRank, q, modPoly)

	p0 := randPoly(maxRank, q)
	pOut := reducer.Reduce(p0)

	assert.Equal(t, reduce(p0.Coeffs[0], vec.Reduce(modPoly, q0), q0), pOut.Coeffs[0])
	assert.Equal(t, reduce(p0.Coeffs[1], vec.Reduce(modPoly, q1), q1), pOut.Coeffs[1])
}

func BenchmarkReducer(b *testing.B) {
	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range reducerBenchLogN {
			N := 1 << logN
			maxRank := 2 * N
			q := []*num.Modulus{num.NewModulus(rSrc.SampleN(num.MaxModulus) | 1)}

			modPoly := randTernaryPoly(N + 1)
			reducer := crt.NewReducer(maxRank, q, modPoly)

			pOut := crt.NewPoly(N, 1)
			p := randPoly(maxRank, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					reducer.ReduceTo(pOut, p)
				}
			})
		}
	})

	b.Run("type=NTT", func(b *testing.B) {
		for _, logN := range reducerBenchLogN {
			N := 1 << logN
			maxRank := 2 * N

			degNext := num.NextProdPower(N, []int{2})
			diffDegNext := num.NextProdPower(2*(maxRank-N)-1, []int{2})
			redParams := dft.NewCyclicParameters(max(degNext, diffDegNext))
			q := dft.FindNearestNTTPrimes(redParams, 60, 1)

			modPoly := randTernaryPoly(N + 1)
			reducer := crt.NewReducer(maxRank, q, modPoly)

			pOut := crt.NewPoly(N, len(q))
			p := randPoly(maxRank, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					reducer.ReduceTo(pOut, p)
				}
			})
		}
	})
}
