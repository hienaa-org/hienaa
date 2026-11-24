package regev

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// packerType is a type for packing algorithm.
type packerType byte

const (
	// typePow2Cyclotomic1Mod4 is a power-of-two cyclotomic ring with modulus 1 mod 4.
	typePow2Cyclotomic1Mod4 packerType = iota
	// typePow2Cyclotomic3Mod4 is a power-of-two cyclotomic ring with modulus 3 mod 4.
	typePow2Cyclotomic3Mod4
	// typeAnyCyclotomic is an arbitrary index cyclotomic ring.
	typeAnyCyclotomic
	// typePow2AutFixed1Mod4 is a power-of-two conjugate invariant ring with modulus 1 mod 4.
	typePow2AutFixed1Mod4
	// typePow2AutFixed3Mod4 is a power-of-two conjugate invariant ring with modulus 3 mod 4.
	typePow2AutFixed3Mod4
	// typePrimeAutFixed is a prime order fixed automorphism subring.
	typePrimeAutFixed
	// typeTrivial is a trivial packer.
	typeTrivial
	// typeNone is an invalid type.
	typeNone
)

// packerTypeOf returns the [packerType] of [dft.RingParameters].
func packerTypeOf(params dft.RingParameters, mod *num.Modulus) packerType {
	switch params.RingType() {
	case dft.TypeCyclotomic:
		if params.Rank() == 1 {
			return typeTrivial
		}

		if num.IsPowerOfTwo(params.CycloIdx()) {
			primes, _ := num.Factor(mod.Value())
			is1Mod4 := true
			for i := range primes {
				if primes[i]%4 != 1 {
					is1Mod4 = false
					break
				}
			}

			if is1Mod4 {
				return typePow2Cyclotomic1Mod4
			} else if len(primes) == 1 && primes[0] != 2 {
				return typePow2Cyclotomic3Mod4
			}
			return typeTrivial
		}

		if dft.IsNTTFriendly(params, mod) {
			return typeAnyCyclotomic
		}
		return typeTrivial

	case dft.TypeAutFixed:
		switch {
		case num.IsPowerOfTwo(params.CycloIdx()):
			primes, _ := num.Factor(mod.Value())
			is1Mod4 := true
			for i := range primes {
				if primes[i]%4 != 1 {
					is1Mod4 = false
					break
				}
			}

			if is1Mod4 {
				return typePow2AutFixed1Mod4
			} else if len(primes) == 1 && primes[0] != 2 {
				return typePow2AutFixed3Mod4
			}
			return typeTrivial

		case num.IsPrime(params.CycloIdx()):
			primes, _ := num.Factor(mod.Value())
			if len(primes) == 1 {
				return typePrimeAutFixed
			}
		}
	}

	return typeNone
}

// Packer packs/unpacks integers in "slots".
//
// Operations usually take two forms: for example,
//   - Pack(v) packs v, allocates a new vector to store the result and returns it.
//   - PackTo(vPack, v) packs v and writes the result to pre-allocated vPack without returning.
//
// Operations panic when inputs are not consistent with the Packer's parameters.
type Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	packerType                packerType
	pow2Cyclotomic1Mod4Packer *pow2Cyclotomic1Mod4Packer
	pow2Cyclotomic3Mod4Packer *pow2Cyclotomic3Mod4Packer
	anyCyclotomicPacker       *anyCyclotomicPacker
	pow2AutFixed1Mod4Packer   *pow2AutFixed1Mod4Packer
	pow2AutFixed3Mod4Packer   *pow2AutFixed3Mod4Packer
	primeAutFixedPacker       *primeAutFixedPacker
	trivialPacker             *trivialPacker
}

// NewPacker creates a new [Packer].
//
// Panics when the ring parameters or modulus are unsupported.
func NewPacker(params dft.RingParameters, mod *num.Modulus) *Packer {
	pck := &Packer{
		params:     params,
		mod:        mod,
		packerType: packerTypeOf(params, mod),
	}

	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		pck.pow2Cyclotomic1Mod4Packer = newPow2Cyclotomic1Mod4Packer(params, mod)
	case typePow2Cyclotomic3Mod4:
		pck.pow2Cyclotomic3Mod4Packer = newPow2Cyclotomic3Mod4Packer(params, mod)
	case typeAnyCyclotomic:
		pck.anyCyclotomicPacker = newAnyCyclotomicPacker(params, mod)
	case typePow2AutFixed1Mod4:
		pck.pow2AutFixed1Mod4Packer = newPow2AutFixed1Mod4Packer(params, mod)
	case typePow2AutFixed3Mod4:
		pck.pow2AutFixed3Mod4Packer = newPow2AutFixed3Mod4Packer(params, mod)
	case typePrimeAutFixed:
		pck.primeAutFixedPacker = newPrimeAutFixedPacker(params, mod)
	case typeTrivial:
		pck.trivialPacker = newTrivialPacker(params)
	default:
		panic("unsupported ring type or modulus")
	}

	return pck
}

// PackLen returns the number of slots.
func (pck *Packer) PackLen() int {
	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		return pck.pow2Cyclotomic1Mod4Packer.packLen
	case typePow2Cyclotomic3Mod4:
		return pck.pow2Cyclotomic3Mod4Packer.packLen
	case typeAnyCyclotomic:
		return pck.anyCyclotomicPacker.packLen
	case typePow2AutFixed1Mod4:
		return pck.pow2AutFixed1Mod4Packer.packLen
	case typePow2AutFixed3Mod4:
		return pck.pow2AutFixed3Mod4Packer.packLen
	case typePrimeAutFixed:
		return pck.primeAutFixedPacker.packLen
	case typeTrivial:
		return pck.trivialPacker.packLen
	default:
		panic("unsupported ring type or modulus")
	}
}

// Cube returns the dimensions of the packing hypercube.
func (pck *Packer) Cube() []int {
	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		return pck.pow2Cyclotomic1Mod4Packer.cube
	case typePow2Cyclotomic3Mod4:
		return pck.pow2Cyclotomic3Mod4Packer.cube
	case typeAnyCyclotomic:
		return pck.anyCyclotomicPacker.cube
	case typePow2AutFixed1Mod4:
		return pck.pow2AutFixed1Mod4Packer.cube
	case typePow2AutFixed3Mod4:
		return pck.pow2AutFixed3Mod4Packer.cube
	case typePrimeAutFixed:
		return pck.primeAutFixedPacker.cube
	case typeTrivial:
		return pck.trivialPacker.cube
	default:
		panic("unsupported ring type or modulus")
	}
}

// CubeGen returns the automorphism generator for each packing dimension.
func (pck *Packer) CubeGen() []uint64 {
	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		return pck.pow2Cyclotomic1Mod4Packer.cubeGen
	case typePow2Cyclotomic3Mod4:
		return pck.pow2Cyclotomic3Mod4Packer.cubeGen
	case typeAnyCyclotomic:
		return pck.anyCyclotomicPacker.cubeGen
	case typePow2AutFixed1Mod4:
		return pck.pow2AutFixed1Mod4Packer.cubeGen
	case typePow2AutFixed3Mod4:
		return pck.pow2AutFixed3Mod4Packer.cubeGen
	case typePrimeAutFixed:
		return pck.primeAutFixedPacker.cubeGen
	case typeTrivial:
		return pck.trivialPacker.cubeGen
	default:
		panic("unsupported ring type or modulus")
	}
}

// Pack returns the packed coefficients of v.
func (pck *Packer) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, pck.params.Rank())
	pck.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (pck *Packer) PackTo(vPack, v []uint64) {
	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		pck.pow2Cyclotomic1Mod4Packer.packTo(vPack, v)
	case typePow2Cyclotomic3Mod4:
		pck.pow2Cyclotomic3Mod4Packer.packTo(vPack, v)
	case typeAnyCyclotomic:
		pck.anyCyclotomicPacker.packTo(vPack, v)
	case typePow2AutFixed1Mod4:
		pck.pow2AutFixed1Mod4Packer.packTo(vPack, v)
	case typePow2AutFixed3Mod4:
		pck.pow2AutFixed3Mod4Packer.packTo(vPack, v)
	case typePrimeAutFixed:
		pck.primeAutFixedPacker.packTo(vPack, v)
	case typeTrivial:
		pck.trivialPacker.packTo(vPack, v)
	}
}

// Unpack returns all slots of vPack.
func (pck *Packer) Unpack(vPack []uint64) []uint64 {
	v := make([]uint64, pck.PackLen())
	pck.UnpackTo(v, vPack)
	return v
}

// UnpackTo unpacks vPack to v.
func (pck *Packer) UnpackTo(v, vPack []uint64) {
	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		pck.pow2Cyclotomic1Mod4Packer.unpackTo(v, vPack)
	case typePow2Cyclotomic3Mod4:
		pck.pow2Cyclotomic3Mod4Packer.unpackTo(v, vPack)
	case typeAnyCyclotomic:
		pck.anyCyclotomicPacker.unpackTo(v, vPack)
	case typePow2AutFixed1Mod4:
		pck.pow2AutFixed1Mod4Packer.unpackTo(v, vPack)
	case typePow2AutFixed3Mod4:
		pck.pow2AutFixed3Mod4Packer.unpackTo(v, vPack)
	case typePrimeAutFixed:
		pck.primeAutFixedPacker.unpackTo(v, vPack)
	case typeTrivial:
		pck.trivialPacker.unpackTo(v, vPack)
	}
}

// RotIdxToAutIdx returns the automorphism index for the rotation indices idx.
// Each rotation index wraps modulo the corresponding dimension in [Packer.Cube].
//
// Panics when the number of indices does not match the packing dimensions.
func (pck *Packer) RotIdxToAutIdx(idx []int) int {
	switch pck.packerType {
	case typePow2Cyclotomic1Mod4:
		return pck.pow2Cyclotomic1Mod4Packer.rotIdxToAutIdx(idx)
	case typePow2Cyclotomic3Mod4:
		return pck.pow2Cyclotomic3Mod4Packer.rotIdxToAutIdx(idx)
	case typeAnyCyclotomic:
		return pck.anyCyclotomicPacker.rotIdxToAutIdx(idx)
	case typePow2AutFixed1Mod4:
		return pck.pow2AutFixed1Mod4Packer.rotIdxToAutIdx(idx)
	case typePow2AutFixed3Mod4:
		return pck.pow2AutFixed3Mod4Packer.rotIdxToAutIdx(idx)
	case typePrimeAutFixed:
		return pck.primeAutFixedPacker.rotIdxToAutIdx(idx)
	case typeTrivial:
		return pck.trivialPacker.rotIdxToAutIdx(idx)
	default:
		panic("unsupported ring type or modulus")
	}
}

type trivialPacker struct {
	rank int

	packLen int
	cube    []int
	cubeGen []uint64
}

// newTrivialPacker creates a new [trivialPacker].
func newTrivialPacker(params dft.RingParameters) *trivialPacker {
	return &trivialPacker{
		rank:    params.Rank(),
		packLen: 1,
		cube:    []int{1},
		cubeGen: []uint64{1},
	}
}

func (pck *trivialPacker) packTo(vPack, v []uint64) {
	if len(vPack) != pck.rank || len(v) != pck.packLen {
		panic("inconsistent inputs")
	}
	vPack[0] = v[0]
	if len(vPack) > 1 {
		clear(vPack[1:])
	}
}

func (pck *trivialPacker) unpackTo(v, vPack []uint64) {
	if len(v) != pck.packLen || len(vPack) != pck.rank {
		panic("inconsistent inputs")
	}
	v[0] = vPack[0]
}

func (pck *trivialPacker) rotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("inconsistent inputs")
	}
	return 1
}
