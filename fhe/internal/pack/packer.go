// Package pack implements SIMD packing/unpacking of integers and complex numbers.
package pack

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// IntPacker is the interface for packing/unpacking integers.
type IntPacker interface {
	// Params returns the ring parameters.
	Params() dft.RingParameters
	// Modulus returns the modulus used for the packing/unpacking.
	Modulus() *num.Modulus
	// PackLen returns the length of the packing/unpacking.
	PackLen() int
	// Pack returns the packing of v.
	Pack(v []uint64) []uint64
	// PackTo packs v to vPack.
	PackTo(vPack, v []uint64)
	// UnPack returns the unpacking of vPack.
	UnPack(vPack []uint64) []uint64
	// UnPackTo unpacks vPack to v.
	UnPackTo(v, vPack []uint64)
	// Cube returns the form of the hypercube structure.
	Cube() []int
	// CubeGen returns the corresponding generator for the hypercube structure.
	CubeGen() []uint64
	// RotIdxToAutIdx converts a rotation index to an automorphism index.
	RotIdxToAutIdx(idx []int) int
}

// NewIntPacker creates a new [IntPacker].
func NewIntPacker(params dft.RingParameters, mod *num.Modulus) IntPacker {
	switch params.RingType() {
	case dft.TypeCyclotomic:
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
				return newPow2CyclotomicMod1Packer(params, mod)
			} else if primes[0] != 2 && len(primes) == 1 {
				return newPow2CyclotomicMod3Packer(params, mod)
			} else {
				return newTrivialPacker(params, mod)
			}

		case dft.IsNTTFriendly(params, mod):
			return newAnyCyclotomicPacker(params, mod)

		default:
			return newTrivialPacker(params, mod)
		}

	case dft.TypeAutFixed:
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
				return newPow2AutFixedMod1Packer(params, mod)
			} else if len(primes) == 1 {
				return newPow2AutFixedMod3Packer(params, mod)
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
