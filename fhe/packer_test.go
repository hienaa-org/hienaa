package fhe_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc      = csprng.NewUniformSamplerWithSeed(nil)
	benchLogN = []int{10, 11, 12, 13, 14, 15, 16}
)

func randMsg(packLen int, q *num.Modulus) []uint64 {
	msg := make([]uint64, packLen)
	for i := range msg {
		msg[i] = rSrc.SampleN(q.Value())
	}
	return msg
}

func testPackerInt(t *testing.T, rP dft.RingParameters, q *num.Modulus) {
	eval := crt.NewOperator(rP, []*num.Modulus{q})
	packer := pack.NewIntPacker(rP, q)
	packLen := packer.PackLen()

	t.Run("Pack", func(t *testing.T) {
		m0 := randMsg(packLen, q)
		m1 := randMsg(packLen, q)

		p1 := &crt.Element{Coeffs: [][]uint64{packer.Pack(m0)}, IsNTT: false}
		p2 := &crt.Element{Coeffs: [][]uint64{packer.Pack(m1)}, IsNTT: false}

		eval.FwdNTTTo(p1, p1)
		eval.FwdNTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1.Coeffs[0])
		assert.Equal(t, vec.Mul(m0, m1, q), res)
	})

	t.Run("Rotate", func(t *testing.T) {
		m := randMsg(packLen, q)

		rotIdx := make([]int, len(packer.Cube()))
		for i := range rotIdx {
			rotIdx[i] = int(rSrc.SampleN(uint64(packer.Cube()[i]-1))) + 1
		}
		autIdx := packer.RotIdxToAutIdx(rotIdx)

		mRef := make([]uint64, packLen)
		outDigits := make([]int, len(rotIdx))
		for i := 0; i < packLen; i++ {
			idxIn := i
			for j := 0; j < len(outDigits); j++ {
				outDigits[j] = (idxIn + rotIdx[j]) % packer.Cube()[j]
				idxIn /= packer.Cube()[j]
			}

			idxOut := outDigits[len(rotIdx)-1]
			for j := len(outDigits) - 2; j >= 0; j-- {
				idxOut *= packer.Cube()[j]
				idxOut += outDigits[j]
			}
			idxOut = idxOut % packLen

			mRef[idxOut] = m[i]
		}

		p := &crt.Element{Coeffs: [][]uint64{packer.Pack(m)}, IsNTT: false}
		pOut := eval.Aut(p, autIdx)
		res := packer.UnPack(pOut.Coeffs[0])

		assert.Equal(t, mRef, res)
	})
}

func TestPackerInt(t *testing.T) {
	t.Run("type=CyclotomicPow2Mod1", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 1)
			if num.IsPrime(prime) && rP.CycloIndex()/int(num.Order(prime, num.NewModulus(rP.CycloIndex()))) >= 8 {
				break
			}
		}
		q := num.NewModulus(prime * prime)
		testPackerInt(t, rP, q)
	})

	t.Run("type=CyclotomicPow2Mod3", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 3)
			if num.IsPrime(prime) && rP.CycloIndex()/int(num.Order(prime, num.NewModulus(rP.CycloIndex()))) >= 8 {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		testPackerInt(t, rP, q)
	})

	t.Run("type=CyclotomicAnyNTT", func(t *testing.T) {
		sqrtN := int(math.Sqrt(math.Exp2(10)))
		m0 := num.MustNextPrime(sqrtN, 1)
		m1 := num.MustNextPrime(m0, 2)
		M := m0 * m1
		rP := dft.NewCyclotomicParameters(M)

		q := dft.MustFindNextNTTPrimes(rP, 20, 1)[0]

		testPackerInt(t, rP, q)
	})

	t.Run("type=AutFixedPow2Mod1", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(9)+5)
		rP := dft.NewAutFixedParameters(N<<2, N)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 1)
			if num.IsPrime(prime) && rP.CycloIndex()/int(num.Order(prime, num.NewModulus(rP.CycloIndex()))) >= 16 {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		testPackerInt(t, rP, q)
	})

	t.Run("type=AutFixedPow2Mod3", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(9)+5)
		rP := dft.NewAutFixedParameters(N<<2, N)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 3)
			if num.IsPrime(prime) && rP.CycloIndex()/int(num.Order(prime, num.NewModulus(rP.CycloIndex()))) >= 32 {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		testPackerInt(t, rP, q)
	})

	t.Run("type=AutFixedPrime", func(t *testing.T) {
		prime := num.MustNextPrime(int(rSrc.SampleN(5)+1), 1)

		var M int
		var N int
		for {
			M = num.MustNextPrime(int(rSrc.SampleN(1<<12)), 1)
			ord := num.Order(uint64(prime), num.NewModulus(M))
			if M > 10 && ord < 10 {
				N = (M - 1) / int(ord)
				break
			}
		}

		q := num.NewModulus(num.Exp(uint64(prime), 3, nil))
		rP := dft.NewAutFixedParameters(M, N)

		testPackerInt(t, rP, q)
	})
}

func BenchmarkPackerInt(b *testing.B) {
	b.Run("type=CyclotomicPow2Mod1", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(N << 1)
			q := dft.MustFindNextNTTPrimes(rP, 1, 1)[0]
			packer := pack.NewIntPacker(rP, q)
			packLen := packer.PackLen()

			msg := randMsg(packLen, q)
			p := packer.Pack(msg)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Pack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.PackTo(p, msg)
					}
				})
				b.Run("UnPack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.UnPackTo(msg, p)
					}
				})
			})
		}
	})

	b.Run("type=CyclotomicPow2Mod3", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(N << 2)

			var prime uint64
			for {
				prime = uint64(2*uint64(N)*rSrc.SampleN(1<<20) - 1)
				if num.IsPrime(prime) {
					break
				}
			}

			q := num.NewModulus(prime)
			packer := pack.NewIntPacker(rP, q)
			packLen := packer.PackLen()

			msg := randMsg(packLen, q)
			p := packer.Pack(msg)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Pack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.PackTo(p, msg)
					}
				})
				b.Run("UnPack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.UnPackTo(msg, p)
					}
				})
			})
		}
	})

	b.Run("type=CyclotomicAnyNTT", func(b *testing.B) {
		for _, logN := range benchLogN {
			sqrtN := int(math.Sqrt(math.Exp2(float64(logN))))
			m0 := num.MustNextPrime(sqrtN, 1)
			m1 := num.MustNextPrime(m0, 2)
			M := m0 * m1
			rP := dft.NewCyclotomicParameters(M)
			q := dft.MustFindNextNTTPrimes(rP, 1, 1)[0]
			packer := pack.NewIntPacker(rP, q)
			packLen := packer.PackLen()
			msg := randMsg(packLen, q)
			p := packer.Pack(msg)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Pack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.PackTo(p, msg)
					}
				})
				b.Run("UnPack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.UnPackTo(msg, p)
					}
				})
			})
		}
	})

	b.Run("type=AutFixedPow2Mod1", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewAutFixedParameters(N<<2, N)
			q := dft.MustFindNextNTTPrimes(rP, 1, 1)[0]
			packer := pack.NewIntPacker(rP, q)
			packLen := packer.PackLen()
			msg := randMsg(packLen, q)
			p := packer.Pack(msg)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Pack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.PackTo(p, msg)
					}
				})
				b.Run("UnPack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.UnPackTo(msg, p)
					}
				})
			})
		}
	})

	b.Run("type=AutFixedPow2Mod3", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewAutFixedParameters(N<<2, N)

			var prime uint64
			for {
				prime = uint64(4*uint64(N)*rSrc.SampleN(1<<20) - 1)
				if num.IsPrime(prime) {
					break
				}
			}
			q := num.NewModulus(prime)
			packer := pack.NewIntPacker(rP, q)
			packLen := packer.PackLen()
			msg := randMsg(packLen, q)
			p := packer.Pack(msg)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Pack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.PackTo(p, msg)
					}
				})
				b.Run("UnPack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.UnPackTo(msg, p)
					}
				})
			})
		}
	})

	b.Run("type=AutFixedPrime", func(b *testing.B) {
		for _, logN := range benchLogN {
			prime := 2

			var N int
			M := 1 << logN
			for {
				M = num.MustNextPrime(M, 1)
				ord := num.Order(uint64(prime), num.NewModulus(M))
				N = (M - 1) / int(ord)
				if N > 1<<logN-1000 && N < 1<<logN+1000 {
					break
				}
			}

			q := num.NewModulus(num.Exp(uint64(prime), 3, nil))
			rP := dft.NewAutFixedParameters(M, N)

			packer := pack.NewIntPacker(rP, q)
			packLen := packer.PackLen()

			msg := randMsg(packLen, q)
			p := packer.Pack(msg)

			b.Run(fmt.Sprintf("LogN=%v", logN), func(b *testing.B) {
				b.Run("Pack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.PackTo(p, msg)
					}
				})
				b.Run("UnPack", func(b *testing.B) {
					for i := 0; i < b.N; i++ {
						packer.UnPackTo(msg, p)
					}
				})
			})
		}
	})
}
