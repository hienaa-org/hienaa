// Package pack implements SIMD packing/unpacking of integers and complex numbers.
package pack

import (
	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// packerInt is the interface for packing/unpacking integers.
type PackerInt interface {
	// Params returns the ring parameters.
	Params() dft.RingParameters
	// Modulus returns the modulus used for the packing/unpacking.
	Modulus() *num.Modulus
	// PackLen returns the length of the packing/unpacking.
	PackLen() int
	// SafeCopy returns a thread-safe copy.
	SafeCopy() PackerInt
	// Pack packs the input integer vector.
	Pack(vIn []uint64) []uint64
	// PackTo packs the input integer vector to the polynomial p.
	PackTo(vOut []uint64, vIn []uint64)
	// UnPack unpacks the polynomial p.
	UnPack(vIn []uint64) []uint64
	// UnPackTo unpacks the polynomial p to the output integer vector.
	UnPackTo(vOut []uint64, vIn []uint64)
}

func NewPackerInt(params dft.RingParameters, mod *num.Modulus) PackerInt {
	switch params.RingType() {
	case dft.Cyclotomic:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			primes, _ := num.Factor(mod.Value())
			isMod1 := true
			for _, prime := range primes {
				if prime%4 != 1 {
					isMod1 = false
					break
				}
			}

			if isMod1 {
				return newCyclotomicPow2Mod1Packer(params, mod)
			} else if len(primes) == 1 {
				return newCyclotomicPow2Mod3Packer(params, mod)
			}
		case dft.IsNTTFriendly(params, mod):
			return newCyclotomicAnyNTTPacker(params, mod)
		}
	case dft.AutFixed:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			primes, _ := num.Factor(mod.Value())
			isMod1 := true
			for _, prime := range primes {
				if prime%4 != 1 {
					isMod1 = false
					break
				}
			}
			if isMod1 {
				return newAutFixedPow2Mod1Packer(params, mod)
			} else if len(primes) == 1 {
				return newAutFixedPow2Mod3Packer(params, mod)
			}
		case num.IsPrime(params.CycloOrder()):
			primes, _ := num.Factor(mod.Value())
			if len(primes) == 1 {
				return newAutFixedPrimePacker(params, mod)
			}
		}
	}

	panic("NewPackerInt: unsupported ring type or parameters")
}

// packerBuffer is a buffer for [PackerInt].
type packerBuffer struct {
	coeffs [][]uint64
}

func newPackerBuffer(dim, rank int) packerBuffer {
	coeffs := make([][]uint64, dim)
	for i := range coeffs {
		coeffs[i] = make([]uint64, rank)
	}
	return packerBuffer{coeffs: coeffs}
}

// pow2Mod3PackerBuffer is a buffer for [pow2Mod3Packer].
type pow2Mod3PackerBuffer struct {
	coeffs []gnum.GaussianInt
}

func newPow2Mod3PackerBuffer(rank int) pow2Mod3PackerBuffer {
	coeffs := make([]gnum.GaussianInt, rank)
	for i := range coeffs {
		coeffs[i] = gnum.GaussianInt{Real: 0, Imag: 0}
	}
	return pow2Mod3PackerBuffer{coeffs: coeffs}
}
