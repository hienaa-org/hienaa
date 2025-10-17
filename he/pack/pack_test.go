package pack

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
	benchLogN = []int{12, 13, 14}
)

func randMsg(packLen int, q *num.Modulus) []uint64 {
	msg := make([]uint64, packLen)
	for i := range msg {
		msg[i] = rSrc.SampleN(q.Value())
	}
	return msg
}

func TestPackerInt(t *testing.T) {
	t.Run("type=CyclotomicPow2Mod1", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 1)
			if num.IsPrime(prime) {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		eval := crt.NewPolyEvaluator(rP, []*num.Modulus{q})
		packer := NewPackerInt(rP, q)
		packLen := packer.PackLen()

		m0 := randMsg(packLen, q)
		m1 := randMsg(packLen, q)
		p1 := packer.Pack(m0)
		p2 := packer.Pack(m1)

		eval.NTTTo(p1, p1)
		eval.NTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1)

		assert.Equal(t, vec.Mul(m0, m1, q), res)
	})

	t.Run("type=CyclotomicPow2Mod3", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(10)+5)
		rP := dft.NewCyclotomicParameters(N << 1)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 3)
			if num.IsPrime(prime) {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		eval := crt.NewPolyEvaluator(rP, []*num.Modulus{q})
		packer := NewPackerInt(rP, q)
		packLen := packer.PackLen()

		m0 := randMsg(packLen, q)
		m1 := randMsg(packLen, q)
		p1 := packer.Pack(m0)
		p2 := packer.Pack(m1)

		eval.NTTTo(p1, p1)
		eval.NTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1)

		assert.Equal(t, vec.Mul(m0, m1, q), res)
	})

	t.Run("type=CyclotomicAnyNTT", func(t *testing.T) {
		sqrtN := int(math.Sqrt(math.Exp2(10)))
		m0 := num.NextPrime(sqrtN, 1)
		m1 := num.NextPrime(m0, 2)
		M := m0 * m1
		rP := dft.NewCyclotomicParameters(M)

		q := dft.FindNextNTTPrimes(rP, 20, 1)[0]

		eval := crt.NewPolyEvaluator(rP, []*num.Modulus{q})
		packer := NewPackerInt(rP, q)
		packLen := packer.PackLen()

		msg0 := randMsg(packLen, q)
		msg1 := randMsg(packLen, q)
		p1 := packer.Pack(msg0)
		p2 := packer.Pack(msg1)

		eval.NTTTo(p1, p1)
		eval.NTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1)

		assert.Equal(t, vec.Mul(msg0, msg1, q), res)
	})

	t.Run("type=AutFixedPow2Mod1", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(9)+5)
		rP := dft.NewAutFixedParameters(N<<2, N)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 1)
			if num.IsPrime(prime) {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		eval := crt.NewPolyEvaluator(rP, []*num.Modulus{q})
		packer := NewPackerInt(rP, q)
		packLen := packer.PackLen()

		m0 := randMsg(packLen, q)
		m1 := randMsg(packLen, q)
		p1 := packer.Pack(m0)
		p2 := packer.Pack(m1)

		eval.NTTTo(p1, p1)
		eval.NTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1)

		assert.Equal(t, vec.Mul(m0, m1, q), res)
	})

	t.Run("type=AutFixedPow2Mod3", func(t *testing.T) {
		N := 1 << int(rSrc.SampleN(9)+5)
		rP := dft.NewAutFixedParameters(N<<2, N)

		var prime uint64
		for {
			prime = uint64(4*rSrc.SampleN(1<<20) + 3)
			if num.IsPrime(prime) {
				break
			}
		}
		q := num.NewModulus(prime * prime)

		eval := crt.NewPolyEvaluator(rP, []*num.Modulus{q})
		packer := NewPackerInt(rP, q)
		packLen := packer.PackLen()

		m0 := randMsg(packLen, q)
		m1 := randMsg(packLen, q)
		p1 := packer.Pack(m0)
		p2 := packer.Pack(m1)

		eval.NTTTo(p1, p1)
		eval.NTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1)

		assert.Equal(t, vec.Mul(m0, m1, q), res)
	})

	t.Run("type=AutFixedPrime", func(t *testing.T) {
		prime := num.NextPrime(int(rSrc.SampleN(20)+1), 1)

		var M int
		var N int
		for {
			M = num.NextPrime(int(rSrc.SampleN(1<<10)), 1)
			ord := num.Order(uint64(prime), num.NewModulus(M))
			if ord < 10 {
				N = (M - 1) / int(ord)
				break
			}
		}

		q := num.NewModulus(num.Exp(uint64(prime), 3, nil))
		rP := dft.NewAutFixedParameters(M, N)

		eval := crt.NewPolyEvaluator(rP, []*num.Modulus{q})
		packer := NewPackerInt(rP, q)
		packLen := packer.PackLen()

		msg0 := randMsg(packLen, q)
		msg1 := randMsg(packLen, q)
		p1 := packer.Pack(msg0)
		p2 := packer.Pack(msg1)

		eval.NTTTo(p1, p1)
		eval.NTTTo(p2, p2)
		eval.MulTo(p1, p1, p2)
		eval.InvNTTTo(p1, p1)

		res := packer.UnPack(p1)

		assert.Equal(t, vec.Mul(msg0, msg1, q), res)
	})
}

func BenchmarkPackerInt(b *testing.B) {
	b.Run("type=CyclotomicPow2Mod1", func(b *testing.B) {
		for _, logN := range benchLogN {
			N := 1 << logN
			rP := dft.NewCyclotomicParameters(N << 1)
			q := dft.FindNextNTTPrimes(rP, 1, 1)[0]
			packer := NewPackerInt(rP, q)
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
			packer := NewPackerInt(rP, q)
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
			m0 := num.NextPrime(sqrtN, 1)
			m1 := num.NextPrime(m0, 2)
			M := m0 * m1
			rP := dft.NewCyclotomicParameters(M)
			q := dft.FindNextNTTPrimes(rP, 1, 1)[0]
			packer := NewPackerInt(rP, q)
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
			q := dft.FindNextNTTPrimes(rP, 1, 1)[0]
			packer := NewPackerInt(rP, q)
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
			packer := NewPackerInt(rP, q)
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
				M = num.NextPrime(M, 1)
				ord := num.Order(uint64(prime), num.NewModulus(M))
				N = (M - 1) / int(ord)
				if N > 1<<logN-1000 && N < 1<<logN+1000 {
					break
				}
			}

			q := num.NewModulus(num.Exp(uint64(prime), 3, nil))
			rP := dft.NewAutFixedParameters(M, N)

			packer := NewPackerInt(rP, q)
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
