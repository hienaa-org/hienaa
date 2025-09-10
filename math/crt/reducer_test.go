package crt_test

import (
	"fmt"
	"testing"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{12, 13, 14, 15, 16, 17}
)

func randTernaryPoly(rank int) []int64 {
	p := make([]int64, rank)
	for i := 0; i < rank-1; i++ {
		p[i] = int64(rSrc.SampleN(3)) - 1
	}
	p[rank-1] = 1
	return p
}

func randPoly(rank int, q *num.Modulus) *crt.Poly {
	p := crt.NewPoly(rank, 1)
	for i := 0; i < rank; i++ {
		p.Coeffs[0][i] = rSrc.SampleN(q.Value())
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

			for j := 0; j < len(p1); j++ {
				rem[len(rem)-i-j-1] = num.Sub(rem[len(rem)-i-j-1], num.Mul(p1[len(p1)-j-1], quo[len(quo)-i-1], q), q)
			}
		}
	}

	return rem
}

func TestReducer(t *testing.T) {
	t.Run("type=Any", func(t *testing.T) {
		deg := int(rSrc.SampleN(1<<12)) + 4
		maxRank := deg + int(rSrc.SampleN(1<<12)) + 1
		q := num.NewModulus(rSrc.SampleN(1 << 61))

		quoPoly := randTernaryPoly(deg + 1)
		reducer := crt.NewReducer(maxRank, q, quoPoly)

		p := randPoly(maxRank, q)
		pOut := reducer.Reduce(p)

		modPoly := make([]uint64, len(quoPoly))
		for i := range modPoly {
			modPoly[i] = num.Reduce(uint64(quoPoly[i])+q.Value(), q)
		}

		assert.Equal(t, reduce(p.Coeffs[0], modPoly, q), pOut.Coeffs[0])
	})

	t.Run("type=NTT", func(t *testing.T) {
		deg := int(rSrc.SampleN(1<<12)) + 4
		maxRank := deg + int(rSrc.SampleN(1<<12)) + 1

		degNext := num.NextProdPower(deg, []int{2})
		diffDegNext := num.NextProdPower(2*(maxRank-deg)-1, []int{2})
		q := dft.FindNearestNTTPrimes(dft.NewCyclicParameters(max(degNext, diffDegNext)), 60, 1)[0]

		quoPoly := randTernaryPoly(deg + 1)
		reducer := crt.NewReducer(maxRank, q, quoPoly)

		p := randPoly(maxRank, q)
		pOut := reducer.Reduce(p)

		modPoly := make([]uint64, len(quoPoly))
		for i := range modPoly {
			modPoly[i] = num.Reduce(uint64(quoPoly[i])+q.Value(), q)
		}

		assert.Equal(t, reduce(p.Coeffs[0], modPoly, q), pOut.Coeffs[0])
	})
}

func BenchmarkReducer(b *testing.B) {
	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range benchLogN {
			deg := 1 << logN
			maxRank := deg + 1<<logN
			q := num.NewModulus(rSrc.SampleN(1 << 61))

			quoPoly := randTernaryPoly(deg)
			reducer := crt.NewReducer(maxRank, q, quoPoly)

			p := randPoly(maxRank, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					reducer.Reduce(p)
				}
			})
		}
	})

	b.Run("type=NTT", func(b *testing.B) {
		for _, logN := range benchLogN {
			deg := 1 << logN
			maxRank := deg + 1<<logN

			degNext := num.NextProdPower(deg, []int{2})
			diffDegNext := num.NextProdPower(2*(maxRank-deg)-1, []int{2})
			q := dft.FindNearestNTTPrimes(dft.NewCyclicParameters(max(degNext, diffDegNext)), 60, 1)[0]

			quoPoly := randTernaryPoly(deg + 1)
			reducer := crt.NewReducer(maxRank, q, quoPoly)

			p := randPoly(maxRank, q)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					reducer.Reduce(p)
				}
			})
		}
	})
}
