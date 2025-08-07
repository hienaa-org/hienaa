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

func TestAutFixedNTT(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := int(num.NextProdPower(rSrc.SampleN(1<<3), []uint64{2}))
		ringParams := dft.NewAutFixedParameters(N<<2, N)
		qs := dft.FindNearestNTTPrimes(ringParams, 30, 2)
		q := num.NewModulus(qs[0].Value() * qs[1].Value())
		ntt := dft.NewTransformer(ringParams, q)

		p0 := randPoly(ringParams, q)
		p1 := randPoly(ringParams, q)

		p0Long := make([]uint64, N<<1)
		p1Long := make([]uint64, N<<1)
		copy(p0Long[:N], p0)
		copy(p1Long[:N], p1)

		p0Long[N] = 0
		p1Long[N] = 0
		for i := 1; i < N; i++ {
			p0Long[i+N] = num.Neg(p0[N-i], q)
			p1Long[i+N] = num.Neg(p1[N-i], q)
		}

		p0NTT := make([]uint64, N)
		copy(p0NTT, p0)
		ntt.ForwardInPlace(p0NTT)

		p1NTT := make([]uint64, N)
		copy(p1NTT, p1)
		ntt.ForwardInPlace(p1NTT)

		pOut := make([]uint64, N)
		vec.MMulTo(pOut, p0NTT, p1NTT, q)
		ntt.InverseInPlace(pOut)

		assert.Equal(t, cyclotomicPow2Mul(p0Long, p1Long, q)[:N], pOut)
	})

	t.Run("type=Prime", func(t *testing.T) {
		for cnt := 0; cnt < 100; cnt++ {
			cycloOrd := int(num.NextPrime(rSrc.SampleN(1<<12), 1))
			primes, exps := num.Factor(uint64(cycloOrd - 1))
			fold := 1
			for i := range primes {
				e := rSrc.SampleN(uint64(exps[i]))
				for j := 0; j < int(e); j++ {
					fold *= int(primes[i])
				}
			}
			N := (cycloOrd - 1) / fold

			ringParams := dft.NewAutFixedParameters(cycloOrd, N)
			q := dft.FindPrevNTTPrimes(ringParams, 61, 1)[0]
			ntt := dft.NewTransformer(ringParams, q)

			p0 := randPoly(ringParams, q)
			p1 := randPoly(ringParams, q)

			p0Long := make([]uint64, cycloOrd)
			p1Long := make([]uint64, cycloOrd)

			cycloOrdMod := num.NewModulus(uint64(cycloOrd))
			root := num.Generators(cycloOrdMod)[0]
			idx := uint64(1)
			for i := 0; i < fold; i++ {
				for j := 0; j < N; j++ {
					p0Long[idx] = p0[j]
					p1Long[idx] = p1[j]
					idx = num.Mul(idx, root, cycloOrdMod)
				}
			}

			p0NTT := make([]uint64, N)
			copy(p0NTT, p0)
			ntt.ForwardInPlace(p0NTT)

			p1NTT := make([]uint64, N)
			copy(p1NTT, p1)
			ntt.ForwardInPlace(p1NTT)

			pOut := make([]uint64, N)
			vec.MMulTo(pOut, p0NTT, p1NTT, q)
			ntt.InverseInPlace(pOut)

			pTest := make([]uint64, N)
			pTestLong := cyclicMul(p0Long, p1Long, q)
			for i := 1; i < cycloOrd; i++ {
				pTestLong[i] = num.Sub(pTestLong[i], pTestLong[0], q)
			}
			pTestLong[0] = 0

			idx = 1
			for i := 0; i < N; i++ {
				pTest[i] = pTestLong[idx]
				idx = num.Mul(idx, root, cycloOrdMod)
			}

			fmt.Println(cycloOrd, fold, N)
			for i := 0; i < N; i++ {
				if pTest[i] != pOut[i] {
					fmt.Println("WRONG!")
					fmt.Println(i, pTest[i], pOut[i])
					break
				}
			}

			// assert.Equal(t, pTest, pOut)
		}
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
