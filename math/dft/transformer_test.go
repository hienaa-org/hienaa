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

func mulReduce(p0, p1, pMod []uint64, q *num.Modulus) []uint64 {
	pMul := mul(p0, p1, q)
	return reduce(pMul, pMod, q)
}

func mul(p0, p1 []uint64, q *num.Modulus) []uint64 {
	pOut := make([]uint64, len(p0)+len(p1)-1)
	for i := range p0 {
		for j := range p1 {
			pOut[i+j] = num.Add(pOut[i+j], num.Mul(p0[i], p1[j], q), q)
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
			vec.MulSubScalarTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], q)
		}
	}

	return rem[:len(p1)-1]
}

func expandAutFixedPoly(params dft.RingParameters, p []uint64, q *num.Modulus) []uint64 {
	var pFull []uint64
	if num.IsPowerOfTwo(params.CycloIndex()) {
		pFull = make([]uint64, params.Rank()<<1)
		copy(pFull[:params.Rank()], p)
		pFull[params.Rank()] = 0
		for i := 1; i < params.Rank(); i++ {
			pFull[i+params.Rank()] = num.Neg(p[params.Rank()-i], q)
		}
	} else {
		pFull = make([]uint64, params.CycloIndex())
		cycloIdxMod := num.NewModulus(params.CycloIndex())
		root := num.Generators(cycloIdxMod)[0]
		idx := uint64(1)
		for i := 0; i < (params.CycloIndex()-1)/params.Rank(); i++ {
			for j := 0; j < params.Rank(); j++ {
				pFull[idx] = p[j]
				idx = num.Mul(idx, root, cycloIdxMod)
			}
		}
		for i := 0; i < params.CycloIndex()-1; i++ {
			pFull[i] = num.Sub(pFull[i], pFull[params.CycloIndex()-1], q)
		}
		pFull = pFull[:params.CycloIndex()-1]
	}
	return pFull
}

func testNTT(t *testing.T, params dft.RingParameters) {
	N := params.Rank()
	qs := dft.MustFindPrevNTTPrimes(params, float64((num.MaxModulusBits>>1)-1), 2)
	q := num.NewModulus(qs[0].Value() * qs[1].Value())

	ntt := dft.NewTransformer(params, q)

	p0 := randPoly(params, q)
	p1 := randPoly(params, q)

	var pMod []uint64
	switch params.RingType() {
	case dft.TypeCyclotomic, dft.TypeAutFixed:
		pMod = vec.Reduce(params.ModulusPoly(), q)
	case dft.TypeCyclic:
		pMod = make([]uint64, N+1)
		pMod[N] = 1
		pMod[0] = q.Value() - 1
	}

	p0NTT := make([]uint64, N)
	p1NTT := make([]uint64, N)
	pOutNTT := make([]uint64, N)
	pOut := make([]uint64, N)

	ntt.ForwardTo(p0NTT, p0)
	ntt.ForwardTo(p1NTT, p1)
	vec.MMulTo(pOutNTT, p0NTT, p1NTT, q)
	ntt.InverseTo(pOut, pOutNTT)

	switch params.RingType() {
	case dft.TypeCyclotomic, dft.TypeCyclic:
		assert.Equal(t, mulReduce(p0, p1, pMod, q), pOut)
	case dft.TypeAutFixed:
		p0Full := expandAutFixedPoly(params, p0, q)
		p1Full := expandAutFixedPoly(params, p1, q)
		pOutFull := expandAutFixedPoly(params, pOut, q)
		assert.Equal(t, mulReduce(p0Full, p1Full, pMod, q), pOutFull)
	}

}

func TestCyclotomicNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testNTT(t, dft.NewCyclotomicParameters(N<<1))
	})

	t.Run("type=Any", func(t *testing.T) {
		sqrtN := int(math.Sqrt(math.Exp2(10)))
		m0 := num.MustNextPrime(sqrtN, 1)
		m1 := num.MustNextPrime(m0, 2)
		M := m0 * m1

		testNTT(t, dft.NewCyclotomicParameters(M))
	})
}

func TestCyclicNTT(t *testing.T) {
	t.Run("type=Pow235", func(t *testing.T) {
		N := num.NextProdPower(int(rSrc.SampleN(1<<10)), []int{2, 3, 5})

		testNTT(t, dft.NewCyclicParameters(N))
	})

	t.Run("type=Any", func(t *testing.T) {
		var N int
		for {
			N = int(rSrc.SampleN(1 << 10))
			if !num.IsProdPowerOf(N, []int{2, 3, 5}) {
				break
			}
		}

		testNTT(t, dft.NewCyclicParameters(N))
	})
}

func TestAutFixedNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testNTT(t, dft.NewAutFixedParameters(N<<2, N))
	})

	t.Run("type=Prime", func(t *testing.T) {
		M := num.MustNextPrime(int(rSrc.SampleN(1<<10)), 1)
		primes, exps := num.Factor(M - 1)
		fold := 1
		for i := range primes {
			e := int(rSrc.SampleN(uint64(exps[i])))
			for j := 0; j < e; j++ {
				fold *= primes[i]
			}
		}
		N := (M - 1) / fold

		testNTT(t, dft.NewAutFixedParameters(M, N))
	})
}

func benchmarkNTT(b *testing.B, params dft.RingParameters) {
	q := dft.MustFindPrevNTTPrimes(params, float64(50), 1)[0]
	ntt := dft.NewTransformer(params, q)

	p := randPoly(params, q)
	pNTT := randPoly(params, q)

	b.Run("FwdNTT", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ntt.ForwardTo(pNTT, p)
		}
	})

	b.Run("FwdNTTInPlace", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ntt.ForwardTo(p, p)
		}
	})

	b.Run("InvNTT", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ntt.InverseTo(p, pNTT)
		}
	})

	b.Run("InvNTTInPlace", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ntt.InverseTo(pNTT, pNTT)
		}
	})
}

func BenchmarkCyclotomicNTT(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkNTT(b, dft.NewCyclotomicParameters(1<<(logN+1)))
			})
		}
	})

	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range benchLogN {
			sqrtN := int(math.Sqrt(math.Exp2(float64(logN))))
			m0 := num.MustNextPrime(sqrtN, 1)
			m1 := num.MustNextPrime(m0, 2)
			M := m0 * m1

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkNTT(b, dft.NewCyclotomicParameters(M))
			})
		}
	})
}

func BenchmarkCyclicNTT(b *testing.B) {
	b.Run("type=Pow235", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := num.NextProdPower((1<<logN)+1, []int{2, 3, 5})
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkNTT(b, dft.NewCyclicParameters(N))
			})
		}
	})

	b.Run("type=Any", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := (1 << logN) + 1

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkNTT(b, dft.NewCyclicParameters(N))
			})
		}
	})
}

func BenchmarkAutFixedNTT(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkNTT(b, dft.NewAutFixedParameters(N<<2, N))
			})
		}
	})

	b.Run("type=Prime", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			M := num.MustNextPrime(1, N)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkNTT(b, dft.NewAutFixedParameters(M, N))
			})
		}
	})
}
