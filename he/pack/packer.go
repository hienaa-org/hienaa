// Package pack implements SIMD packing/unpacking of integers and complex numbers.
package pack

import (
	"github.com/hienaa-org/hienaa/math/crt"
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
	Pack(vIn []uint64) *crt.Poly
	// PackTo packs the input integer vector to the polynomial p.
	PackTo(pOut *crt.Poly, vIn []uint64)
	// UnPack unpacks the polynomial p.
	UnPack(pIn *crt.Poly) []uint64
	// UnPackTo unpacks the polynomial p to the output integer vector.
	UnPackTo(vOut []uint64, pIn *crt.Poly)
}

func NewPackerInt(params dft.RingParameters, mod *num.Modulus) PackerInt {
	switch params.RingType() {
	case dft.Cyclotomic:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			// TODO: Add Bruun NTT.
			switch {
			case mod.Value()%4 == 1:
				return newCyclotomicPow2NTTPacker(params, mod)
			}
		case dft.IsNTTFriendly(params, mod):
			return newCyclotomicAnyNTTPacker(params, mod)
		}
	case dft.AutFixed:
		switch {
		case num.IsPrime(params.CycloOrder()):
			return newAutFixedPrimePacker(params, mod)
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
