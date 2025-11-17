package fhe

import (
	"github.com/hienaa-org/hienaa/math/num"
)

// Packer pack/unpacks vector of numbers into [Plaintext] in a SIMD manner.
type Packer[Self any, P Plaintext, T num.Number] interface {
	// PackLen returns the length of packable vector.
	PackLen() int
	// SafeCopy returns a thread-safe copy.
	SafeCopy() Packer[Self, P, T]
	// Pack packs v.
	Pack(v []T) P
	// PackTo packs v to ptOut.
	PackTo(ptOut P, v []T)
	// UnPack unpacks pt.
	UnPack(pt P) []T
	// UnPack unpacks pt to vOut.
	UnPackTo(vOut []T, pt P)
}
