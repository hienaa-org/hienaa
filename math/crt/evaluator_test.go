package crt_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/hienaa-org/hienaa/math/crt"
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

func reduce(p []uint64, q *num.Modulus, modPoly []int64) []uint64 {
	reducer := crt.NewLongDivReducer(len(p), []*num.Modulus{q}, modPoly)
	return reducer.Reduce(&crt.Poly{Coeffs: [][]uint64{p}}).Coeffs[0]
}

func TestCyclotomicEvaluator(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10
		rP := dft.NewCyclotomicParameters(N << 1)

		q := dft.FindPrevNTTPrimes(rP, 40, 1)
		q = append(q, num.NewModulus(num.NextPrime(q[0].Value(), 2)))

		pev := crt.NewPolyEvaluator(rP, q)

		p0 := randPoly(rP.Rank(), q)
		p1 := randPoly(rP.Rank(), q)
		pOut := randPoly(rP.Rank(), q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pOutRef[i] = cyclotomicPow2Mul(p0.Coeffs[i], p1.Coeffs[i], q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclotomicPow2Mul(p0.Coeffs[i], p1.Coeffs[i], q[i])
				vec.AddTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclotomicPow2Mul(p0.Coeffs[i], p1.Coeffs[i], q[i])
				vec.SubTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("Aut", func(t *testing.T) {
			idx := rP.CycloOrder() - 1
			idxInv := int(num.Inv(uint64(idx), num.NewModulus(rP.CycloOrder())))

			pev.AutTo(pOut, p0, idx)
			pev.NTTTo(pOutNTT, pOut)
			pev.AutTo(pOutNTT, pOutNTT, idxInv)
			pev.InvNTTTo(pOut, pOutNTT)

			assert.Equal(t, p0, pOut)
		})
	})

	t.Run("type=Any", func(t *testing.T) {
		M := int(rSrc.SampleN(1 << 10))
		rP := dft.NewCyclotomicParameters(M)
		N := rP.Rank()

		q := dft.FindPrevNTTPrimes(rP, 40, 1)
		q = append(q, num.NewModulus(num.NextPrime(q[0].Value(), 2)))

		cycloSigned := dft.CyclotomicPolynomial(rP.CycloOrder())

		pev := crt.NewPolyEvaluator(rP, q)

		p0 := randPoly(N, q)
		p1 := randPoly(N, q)
		pOut := randPoly(N, q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		p0Ref := make([][]uint64, len(q))
		p1Ref := make([][]uint64, len(q))
		for i := range q {
			p0Ref[i] = make([]uint64, 2*N-1)
			copy(p0Ref[i], p0.Coeffs[i])
			p1Ref[i] = make([]uint64, 2*N-1)
			copy(p1Ref[i], p1.Coeffs[i])
		}

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pMulRef := cyclicMul(p0Ref[i], p1Ref[i], q[i])
				pOutRef[i] = reduce(pMulRef, q[i], cycloSigned)
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0Ref[i], p1Ref[i], q[i])
				vec.AddTo(pOutRef[i], pOutRef[i], reduce(pMulRef, q[i], cycloSigned), q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0Ref[i], p1Ref[i], q[i])
				vec.SubTo(pOutRef[i], pOutRef[i], reduce(pMulRef, q[i], cycloSigned), q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("Aut", func(t *testing.T) {
			idx := rP.CycloOrder() - 1
			idxInv := int(num.Inv(uint64(idx), num.NewModulus(rP.CycloOrder())))

			pev.AutTo(pOut, p0, idx)
			pev.NTTTo(pOut, pOut)
			pev.AutTo(pOut, pOut, idxInv)
			pev.InvNTTTo(pOut, pOut)

			assert.Equal(t, p0.Coeffs, pOut.Coeffs)
		})
	})
}

func TestCyclicEvaluator(t *testing.T) {
	t.Run("type=Pow235", func(t *testing.T) {
		N := num.NextProdPower(int(rSrc.SampleN(1<<10)), []int{2, 3, 5})
		rP := dft.NewCyclicParameters(N)

		q := dft.FindPrevNTTPrimes(rP, 40, 1)
		q = append(q, num.NewModulus(num.NextPrime(q[0].Value(), 2)))

		pev := crt.NewPolyEvaluator(rP, q)

		p0 := randPoly(rP.Rank(), q)
		p1 := randPoly(rP.Rank(), q)
		pOut := randPoly(rP.Rank(), q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pOutRef[i] = cyclicMul(p0.Coeffs[i], p1.Coeffs[i], q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0.Coeffs[i], p1.Coeffs[i], q[i])
				vec.AddTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0.Coeffs[i], p1.Coeffs[i], q[i])
				vec.SubTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})
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

		q := dft.FindPrevNTTPrimes(rP, 40, 1)
		q = append(q, num.NewModulus(num.NextPrime(q[0].Value(), 2)))

		pev := crt.NewPolyEvaluator(rP, q)

		p0 := randPoly(N, q)
		p1 := randPoly(N, q)
		pOut := randPoly(N, q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pOutRef[i] = cyclicMul(p0.Coeffs[i], p1.Coeffs[i], q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0.Coeffs[i], p1.Coeffs[i], q[i])
				vec.AddTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0.Coeffs[i], p1.Coeffs[i], q[i])
				vec.SubTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})
	})
}

func TestAutFixedEvaluator(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10
		rP := dft.NewAutFixedParameters(4*N, N)

		q := dft.FindPrevNTTPrimes(rP, 40, 1)
		q = append(q, num.NewModulus(num.NextPrime(q[0].Value(), 2)))

		pev := crt.NewPolyEvaluator(rP, q)

		p0 := randPoly(rP.Rank(), q)
		p1 := randPoly(rP.Rank(), q)
		pOut := randPoly(rP.Rank(), q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		p0Ref := make([][]uint64, len(q))
		p1Ref := make([][]uint64, len(q))
		for i := range q {
			p0Ref[i] = make([]uint64, 2*N)
			copy(p0Ref[i], p0.Coeffs[i])
			p1Ref[i] = make([]uint64, 2*N)
			copy(p1Ref[i], p1.Coeffs[i])
			for j := 1; j < N; j++ {
				p0Ref[i][j+N] = num.Neg(p0Ref[i][N-j], q[i])
				p1Ref[i][j+N] = num.Neg(p1Ref[i][N-j], q[i])
			}
		}

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pOutRef[i] = cyclotomicPow2Mul(p0Ref[i], p1Ref[i], q[i])[:N]
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclotomicPow2Mul(p0Ref[i], p1Ref[i], q[i])[:N]
				vec.AddTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclotomicPow2Mul(p0Ref[i], p1Ref[i], q[i])[:N]
				vec.SubTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("Aut", func(t *testing.T) {
			idx := rP.CycloOrder() - 3
			idxInv := int(num.Inv(uint64(idx), num.NewModulus(rP.CycloOrder())))

			pev.AutTo(pOut, p0, idx)
			pev.NTTTo(pOutNTT, pOut)
			pev.AutTo(pOutNTT, pOutNTT, idxInv)
			pev.InvNTTTo(pOut, pOutNTT)

			assert.Equal(t, p0.Coeffs, pOut.Coeffs)
		})
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

		q := dft.FindPrevNTTPrimes(rP, 40, 1)
		q = append(q, num.NewModulus(num.NextPrime(q[0].Value(), 2)))

		pev := crt.NewPolyEvaluator(rP, q)

		p0 := randPoly(rP.Rank(), q)
		p1 := randPoly(rP.Rank(), q)
		pOut := randPoly(rP.Rank(), q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		MMod := num.NewModulus(M)
		root := num.Generators(MMod)[0]

		p0Ref := make([][]uint64, len(q))
		p1Ref := make([][]uint64, len(q))
		for i := range q {
			p0Ref[i] = make([]uint64, M)
			p1Ref[i] = make([]uint64, M)
			idx := uint64(1)
			for k := 0; k < fold; k++ {
				for j := 0; j < N; j++ {
					p0Ref[i][idx], p1Ref[i][idx] = p0.Coeffs[i][j], p1.Coeffs[i][j]
					idx = num.Mul(idx, root, MMod)
				}
			}
		}

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pMulRef := cyclicMul(p0Ref[i], p1Ref[i], q[i])
				for j := 1; j < M; j++ {
					pMulRef[j] = num.Sub(pMulRef[j], pMulRef[0], q[i])
				}
				pMulRef[0] = 0

				pOutRef[i] = make([]uint64, N)
				idx := uint64(1)
				for j := 0; j < N; j++ {
					pOutRef[i][j] = pMulRef[idx]
					idx = num.Mul(idx, root, MMod)
				}
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0Ref[i], p1Ref[i], q[i])
				for j := 1; j < M; j++ {
					pMulRef[j] = num.Sub(pMulRef[j], pMulRef[0], q[i])
				}
				pMulRef[0] = 0

				idx := uint64(1)
				for j := 0; j < N; j++ {
					pOutRef[i][j] = num.Add(pOutRef[i][j], pMulRef[idx], q[i])
					idx = num.Mul(idx, root, MMod)
				}
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := cyclicMul(p0Ref[i], p1Ref[i], q[i])
				for j := 1; j < M; j++ {
					pMulRef[j] = num.Sub(pMulRef[j], pMulRef[0], q[i])
				}
				pMulRef[0] = 0

				idx := uint64(1)
				for j := 0; j < N; j++ {
					pOutRef[i][j] = num.Sub(pOutRef[i][j], pMulRef[idx], q[i])
					idx = num.Mul(idx, root, MMod)
				}
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("Aut", func(t *testing.T) {
			idx := rP.CycloOrder() - 3
			idxInv := int(num.Inv(uint64(idx), num.NewModulus(rP.CycloOrder())))

			pev.AutTo(pOut, p0, idx)
			pev.NTTTo(pOutNTT, pOut)
			pev.AutTo(pOutNTT, pOutNTT, idxInv)
			pev.InvNTTTo(pOut, pOutNTT)

			assert.Equal(t, p0.Coeffs, pOut.Coeffs)
		})
	})
}

func TestAnyEvaluator(t *testing.T) {
	t.Run("type=Any", func(t *testing.T) {
		N := 1 << 10

		q := []*num.Modulus{num.NewModulus(num.NextPrime(1<<60+1, 2))}

		modPolySigned := randTernaryPoly(N + 1)

		pev := crt.NewPolyEvaluatorWithModPoly(q, modPolySigned)

		p0 := randPoly(N, q)
		p1 := randPoly(N, q)
		pOut := randPoly(N, q)

		p0NTT := pev.NTT(p0)
		p1NTT := pev.NTT(p1)
		pOutNTT := pev.NTT(pOut)

		p0Ref := make([][]uint64, len(q))
		p1Ref := make([][]uint64, len(q))
		for i := range q {
			p0Ref[i] = make([]uint64, 2*N-1)
			copy(p0Ref[i], p0.Coeffs[i])
			p1Ref[i] = make([]uint64, 2*N-1)
			copy(p1Ref[i], p1.Coeffs[i])
		}

		t.Run("Mul", func(t *testing.T) {
			pev.MulTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			pOutRef := make([][]uint64, len(q))
			for i := range q {
				pOutRef[i] = reduce(cyclicMul(p0Ref[i], p1Ref[i], q[i]), q[i], modPolySigned)
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulAdd", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := reduce(cyclicMul(p0Ref[i], p1Ref[i], q[i]), q[i], modPolySigned)
				vec.AddTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})

		t.Run("MulSub", func(t *testing.T) {
			pOutRef := pOut.Copy().Coeffs

			pev.NTTTo(pOutNTT, pOut)
			pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
			pev.InvNTTTo(pOut, pOutNTT)

			for i := range q {
				pMulRef := reduce(cyclicMul(p0Ref[i], p1Ref[i], q[i]), q[i], modPolySigned)
				vec.SubTo(pOutRef[i], pOutRef[i], pMulRef, q[i])
			}

			assert.Equal(t, pOutRef, pOut.Coeffs)
		})
	})
}

func BenchmarkCyclotomicEvaluator(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(2 * N)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-1)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-1)
						}
					})
				})

				b.Run("Mod=Any", func(b *testing.B) {
					qv := uint64(1)<<60 + 1
					for {
						if !dft.IsNTTFriendly(rP, num.NewModulus(qv)) {
							break
						}
						qv = num.NextPrime(qv, 2)
					}
					q := []*num.Modulus{num.NewModulus(qv)}

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-1)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-1)
						}
					})
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

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-1)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-1)
						}
					})
				})

				b.Run("Mod=Any", func(b *testing.B) {
					qv := uint64(1)<<60 + 1
					for {
						if !dft.IsNTTFriendly(rP, num.NewModulus(qv)) {
							break
						}
						qv = num.NextPrime(qv, 2)
					}
					q := []*num.Modulus{num.NewModulus(qv)}

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-1)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-1)
						}
					})
				})
			})
		}
	})
}

func BenchmarkCyclicEvaluator(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclicParameters(N)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})
				})

				b.Run("Mod=Any", func(b *testing.B) {
					qv := uint64(1)<<60 + 1
					for {
						if !dft.IsNTTFriendly(rP, num.NewModulus(qv)) {
							break
						}
						qv = num.NextPrime(qv, 2)
					}
					q := []*num.Modulus{num.NewModulus(qv)}

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})
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

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})
				})

				b.Run("Mod=Any", func(b *testing.B) {
					qv := uint64(1)<<60 + 1
					for {
						if !dft.IsNTTFriendly(rP, num.NewModulus(qv)) {
							break
						}
						qv = num.NextPrime(qv, 2)
					}
					q := []*num.Modulus{num.NewModulus(qv)}

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})
				})
			})
		}
	})
}

func BenchmarkAutFixedEvaluator(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewAutFixedParameters(4*N, N)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-3)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-3)
						}
					})
				})

				b.Run("Mod=Any", func(b *testing.B) {
					qv := uint64(1)<<60 + 1
					for {
						if !dft.IsNTTFriendly(rP, num.NewModulus(qv)) {
							break
						}
						qv = num.NextPrime(qv, 2)
					}
					q := []*num.Modulus{num.NewModulus(qv)}

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-3)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-3)
						}
					})
				})
			})
		}
	})

	b.Run("type=Prime", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			M := num.NextPrime(1, N)
			rP := dft.NewAutFixedParameters(M, N)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Mod=NTT", func(b *testing.B) {
					q := dft.FindPrevNTTPrimes(rP, num.MaxModulusBits, 1)

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-3)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-3)
						}
					})
				})

				b.Run("Mod=Any", func(b *testing.B) {
					qv := uint64(1)<<60 + 1
					for {
						if !dft.IsNTTFriendly(rP, num.NewModulus(qv)) {
							break
						}
						qv = num.NextPrime(qv, 2)
					}
					q := []*num.Modulus{num.NewModulus(qv)}

					pev := crt.NewPolyEvaluator(rP, q)

					p0 := randPoly(rP.Rank(), q)
					p1 := randPoly(rP.Rank(), q)
					pOut := randPoly(rP.Rank(), q)

					p0NTT := pev.NTT(p0)
					p1NTT := pev.NTT(p1)
					pOutNTT := pev.NTT(pOut)

					b.Run("Add", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AddTo(pOut, p0, p1)
						}
					})

					b.Run("Sub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.SubTo(pOut, p0, p1)
						}
					})

					b.Run("Neg", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NegTo(pOut, p0)
						}
					})

					b.Run("NTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.NTTTo(p0NTT, p0)
						}
					})

					b.Run("InvNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.InvNTTTo(p0, p0NTT)
						}
					})

					b.Run("Mul", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulAdd", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("MulSub", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
						}
					})

					b.Run("Aut", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOut, p0, rP.CycloOrder()-3)
						}
					})

					b.Run("AutNTT", func(b *testing.B) {
						for i := 0; i < b.N; i++ {
							pev.AutTo(pOutNTT, p0NTT, rP.CycloOrder()-3)
						}
					})
				})
			})
		}
	})
}

func BenchmarkAnyEvaluator(b *testing.B) {
	for _, logN := range benchLogN {
		b.Run(fmt.Sprintf("type=Any/LogN=%v/Mod=Any", logN), func(b *testing.B) {
			N := 1 << logN
			modPoly := randTernaryPoly(N + 1)

			q := []*num.Modulus{num.NewModulus(1<<60 + 1)}

			pev := crt.NewPolyEvaluatorWithModPoly(q, modPoly)

			p0 := randPoly(N, q)
			p1 := randPoly(N, q)
			pOut := randPoly(N, q)

			p0NTT := pev.NTT(p0)
			p1NTT := pev.NTT(p1)
			pOutNTT := pev.NTT(pOut)

			b.Run("Add", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.AddTo(pOut, p0, p1)
				}
			})

			b.Run("Sub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.SubTo(pOut, p0, p1)
				}
			})

			b.Run("Neg", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.NegTo(pOut, p0)
				}
			})

			b.Run("NTT", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.NTTTo(p0NTT, p0)
				}
			})

			b.Run("InvNTT", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.InvNTTTo(p0, p0NTT)
				}
			})

			b.Run("Mul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.MulTo(pOutNTT, p0NTT, p1NTT)
				}
			})

			b.Run("MulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.MulAddTo(pOutNTT, p0NTT, p1NTT)
				}
			})

			b.Run("MulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					pev.MulSubTo(pOutNTT, p0NTT, p1NTT)
				}
			})
		})
	}
}
