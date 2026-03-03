package fhe

import (
	"github.com/hienaa-org/hienaa/math/num"
)

// Packer pack/unpacks vector of numbers into [Plaintext] in a SIMD manner.
type Packer[P Plaintext, T num.Number] interface {
	// PackLen returns the length of packable vector.
	PackLen() int
	// Pack packs v.
	Pack(v []T) P
	// PackTo packs v to ptOut.
	PackTo(ptOut P, v []T)
	// UnPack unpacks pt.
	UnPack(pt P) []T
	// UnPack unpacks pt to vOut.
	UnPackTo(vOut []T, pt P)
	// Cube returns the form of the hypercube structure.
	Cube() []int
	// CubeGen returns the corresponding generator for the hypercube structure.
	CubeGen() []uint64
	// RotIdxToAutIdx converts a rotation index to an automorphism index.
	RotIdxToAutIdx(idx []int) int
}
