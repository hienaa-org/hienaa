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

func mulReduce(p0, p1, pMod [][]uint64, q []*num.Modulus) [][]uint64 {
	pOut := make([][]uint64, len(q))
	for i := range q {
		pOut[i] = reduce(mul(p0[i], p1[i], q[i]), pMod[i], q[i])
	}
	return pOut
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
	if num.IsPowerOfTwo(params.CycloOrder()) {
		pFull = make([]uint64, params.Rank()<<1)
		copy(pFull[:params.Rank()], p)
		pFull[params.Rank()] = 0
		for i := 1; i < params.Rank(); i++ {
			pFull[i+params.Rank()] = num.Neg(p[params.Rank()-i], q)
		}
	} else {
		pFull = make([]uint64, params.CycloOrder())
		cycloOrdMod := num.NewModulus(params.CycloOrder())
		root := num.Generators(cycloOrdMod)[0]
		idx := uint64(1)
		for i := 0; i < (params.CycloOrder()-1)/params.Rank(); i++ {
			for j := 0; j < params.Rank(); j++ {
				pFull[idx] = p[j]
				idx = num.Mul(idx, root, cycloOrdMod)
			}
		}
		for i := 0; i < params.CycloOrder()-1; i++ {
			pFull[i] = num.Sub(pFull[i], pFull[params.CycloOrder()-1], q)
		}
		pFull = pFull[:params.CycloOrder()-1]
	}
	return pFull
}

func testOperator(t *testing.T, params dft.RingParameters, modPoly []int64) {
	var op crt.Operator
	var q []*num.Modulus
	if params.RingType() != dft.TypeOther {
		q = dft.MustFindPrevNTTPrimes(params, 40, 1)
		q = append(q, num.NewModulus(num.MustNextPrime(q[0].Value(), 2)))
		op = crt.NewOperator(params, q)
	} else {
		q = []*num.Modulus{num.NewModulus(num.MustNextPrime(1<<40, 1))}
		op = crt.NewOperatorWithModPoly(q, modPoly)
	}

	p0 := randPoly(params.Rank(), q)
	p1 := randPoly(params.Rank(), q)

	pMod := make([][]uint64, len(q))
	switch params.RingType() {
	case dft.TypeCyclotomic, dft.TypeAutFixed:
		cycloPoly := params.ModulusPoly()
		for i := range q {
			pMod[i] = vec.Reduce(cycloPoly, q[i])
		}
	case dft.TypeCyclic:
		for i := range q {
			pMod[i] = make([]uint64, params.Rank()+1)
			pMod[i][params.Rank()] = 1
			pMod[i][0] = q[i].Value() - 1
		}
	case dft.TypeOther:
		for i := range q {
			pMod[i] = vec.Reduce(modPoly, q[i])
		}
	}

	p0NTT := op.FwdNTT(p0)
	p1NTT := op.FwdNTT(p1)

	t.Run("Mul", func(t *testing.T) {
		pOutNTT := op.Mul(p0NTT, p1NTT)
		pOut := op.InvNTT(pOutNTT)

		switch params.RingType() {
		case dft.TypeCyclotomic, dft.TypeCyclic, dft.TypeOther:
			assert.Equal(t, mulReduce(p0.Coeffs, p1.Coeffs, pMod, q), pOut.Coeffs)
		case dft.TypeAutFixed:
			p0Full := make([][]uint64, len(q))
			p1Full := make([][]uint64, len(q))
			pOutFull := make([][]uint64, len(q))
			for i := range q {
				p0Full[i] = expandAutFixedPoly(params, p0.Coeffs[i], q[i])
				p1Full[i] = expandAutFixedPoly(params, p1.Coeffs[i], q[i])
				pOutFull[i] = expandAutFixedPoly(params, pOut.Coeffs[i], q[i])
			}
			assert.Equal(t, mulReduce(p0Full, p1Full, pMod, q), pOutFull)
		}
	})

	t.Run("MulAdd", func(t *testing.T) {
		pOutNTT := randPoly(params.Rank(), q)
		pOutNTT.IsNTT = true
		pOutNTTRef := pOutNTT.Copy()

		op.MulAddTo(pOutNTT, p0NTT, p1NTT)
		op.AddTo(pOutNTTRef, pOutNTTRef, op.Mul(p0NTT, p1NTT))

		assert.Equal(t, pOutNTTRef.Coeffs, pOutNTT.Coeffs)
	})

	t.Run("MulSub", func(t *testing.T) {
		pOutNTT := randPoly(params.Rank(), q)
		pOutNTT.IsNTT = true
		pOutNTTRef := pOutNTT.Copy()

		op.MulSubTo(pOutNTT, p0NTT, p1NTT)
		op.SubTo(pOutNTTRef, pOutNTTRef, op.Mul(p0NTT, p1NTT))

		assert.Equal(t, pOutNTTRef.Coeffs, pOutNTT.Coeffs)
	})

	switch params.RingType() {
	case dft.TypeCyclotomic, dft.TypeAutFixed:
		t.Run("Aut", func(t *testing.T) {
			var idx uint64
			for {
				idx = rSrc.SampleN(uint64(params.CycloOrder()))
				if op.CanAut(int(idx)) {
					break
				}
			}
			idxInv := num.Inv(idx, num.NewModulus(params.CycloOrder()))

			pOut := op.Aut(p0, int(idx))
			op.FwdNTTTo(pOut, pOut)
			op.AutTo(pOut, pOut, int(idxInv))
			op.InvNTTTo(pOut, pOut)

			assert.Equal(t, p0.Coeffs, pOut.Coeffs)
		})
	}
}

func TestCyclotomicOperator(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testOperator(t, dft.NewCyclotomicParameters(N<<1), nil)
	})

	t.Run("type=Any", func(t *testing.T) {
		M := int(rSrc.SampleN(1 << 10))

		testOperator(t, dft.NewCyclotomicParameters(M), nil)
	})
}

func TestCyclicOperator(t *testing.T) {
	t.Run("type=Pow235", func(t *testing.T) {
		N := num.NextProdPower(int(rSrc.SampleN(1<<10)), []int{2, 3, 5})

		testOperator(t, dft.NewCyclicParameters(N), nil)
	})

	t.Run("type=Any", func(t *testing.T) {
		var N int
		for {
			N = int(rSrc.SampleN(1 << 10))
			if !num.IsProdPowerOf(N, []int{2, 3, 5}) {
				break
			}
		}

		testOperator(t, dft.NewCyclicParameters(N), nil)
	})
}

func TestAutFixedOperator(t *testing.T) {
	t.Run("type=Pow2", func(t *testing.T) {
		N := 1 << 10

		testOperator(t, dft.NewAutFixedParameters(N<<2, N), nil)
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

		testOperator(t, dft.NewAutFixedParameters(M, N), nil)
	})
}

func TestAnyOperator(t *testing.T) {
	t.Run("type=Any", func(t *testing.T) {
		N := 1 << 10
		modPolySigned := randTernaryPoly(N + 1)

		testOperator(t, dft.NewOtherParameters(modPolySigned), modPolySigned)
	})
}

func benchmarkOperator(b *testing.B, params dft.RingParameters, modPoly []int64) {
	var modTypes []string
	if params.RingType() != dft.TypeOther {
		modTypes = []string{"NTT", "Any"}
	} else {
		modTypes = []string{"Any"}
	}

	for _, modType := range modTypes {
		var q []*num.Modulus
		if params.RingType() != dft.TypeOther {
			if modType == "NTT" {
				q = dft.MustFindPrevNTTPrimes(params, num.MaxModulusBits, 1)
			} else {
				q = make([]*num.Modulus, 1)
				for {
					q[0] = num.NewModulus(rSrc.SampleN(num.MaxModulus))
					if !dft.IsNTTFriendly(params, q[0]) {
						break
					}
				}
			}
		} else {
			q = []*num.Modulus{num.NewModulus(1<<60 + 1)}
		}

		b.Run(fmt.Sprintf("Mod=%v", modType), func(b *testing.B) {
			var op crt.Operator
			if params.RingType() != dft.TypeOther {
				op = crt.NewOperator(params, q)
			} else {
				op = crt.NewOperatorWithModPoly(q, modPoly)
			}

			p0 := randPoly(params.Rank(), q)
			p1 := randPoly(params.Rank(), q)
			pOut := randPoly(params.Rank(), q)

			p0NTT := op.FwdNTT(p0)
			p1NTT := op.FwdNTT(p1)
			pOutNTT := op.FwdNTT(pOut)

			b.Run("Add", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.AddTo(pOut, p0, p1)
				}
			})

			b.Run("Sub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.SubTo(pOut, p0, p1)
				}
			})

			b.Run("Neg", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.NegTo(pOut, p0)
				}
			})

			b.Run("FwdNTT", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.FwdNTTTo(p0NTT, p0)
				}
			})

			b.Run("InvNTT", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.InvNTTTo(p0, p0NTT)
				}
			})

			b.Run("Mul", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.MulTo(pOutNTT, p0NTT, p1NTT)
				}
			})

			b.Run("MulAdd", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.MulAddTo(pOutNTT, p0NTT, p1NTT)
				}
			})

			b.Run("MulSub", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					op.MulSubTo(pOutNTT, p0NTT, p1NTT)
				}
			})

			switch params.RingType() {
			case dft.TypeCyclotomic, dft.TypeAutFixed:
				var idx uint64
				for {
					idx = rSrc.SampleN(uint64(params.CycloOrder()))
					if op.CanAut(int(idx)) {
						break
					}
				}

				b.Run("Aut", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						op.AutTo(pOut, p0, int(idx))
					}
				})

				b.Run("AutNTT", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						op.AutTo(pOutNTT, p0NTT, int(idx))
					}
				})
			}
		})
	}
}

func BenchmarkCyclotomicOperator(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkOperator(b, dft.NewCyclotomicParameters(N<<1), nil)
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
				benchmarkOperator(b, dft.NewCyclotomicParameters(M), nil)
			})
		}
	})
}

func BenchmarkCyclicOperator(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkOperator(b, dft.NewCyclicParameters(N), nil)
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
				benchmarkOperator(b, dft.NewCyclicParameters(M), nil)
			})
		}
	})
}

func BenchmarkAutFixedOperator(b *testing.B) {
	b.Run("type=Pow2", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkOperator(b, dft.NewAutFixedParameters(N<<2, N), nil)
			})
		}
	})

	b.Run("type=Prime", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			M := num.MustNextPrime(1, N)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				benchmarkOperator(b, dft.NewAutFixedParameters(M, N), nil)
			})
		}
	})
}

func BenchmarkAnyOperator(b *testing.B) {
	for _, logN := range benchLogN {
		N := 1 << logN
		modPoly := randTernaryPoly(N + 1)

		b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
			benchmarkOperator(b, dft.NewOtherParameters(modPoly), modPoly)
		})
	}
}
