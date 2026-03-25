package ckks

import (
	"github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

type RealPacker struct {
	params dft.RingParameters
	pack   pack.ComplexPacker
}

// NewRealPacker creates a new [RealPacker].
func NewRealPacker(params dft.RingParameters) *RealPacker {
	return &RealPacker{
		params: params,
		pack:   pack.NewComplexPacker(params),
	}
}

// ComplexPacker packs a vector of complex numbers into [*Plaintext].
type ComplexPacker struct {
	params rlwe.Parameters
	pack   pack.ComplexPacker
}

// NewComplexPacker creates a new [ComplexPacker].
func NewComplexPacker(params rlwe.Parameters) *ComplexPacker {
	rP := params.RingParams()
	switch rP.RingType() {
	case dft.TypeAutFixed:
		if num.IsPowerOfTwo(rP.CycloOrder()) || ((rP.CycloOrder()-1)/rP.Rank())&1 == 0 {
			panic("unsupported ring type.")
		}
	}

	return &ComplexPacker{
		params: params,
		pack:   pack.NewComplexPacker(params.RingParams()),
	}
}

// Params returns the parameters.
func (p *ComplexPacker) Params() rlwe.Parameters {
	return p.params
}

// PackLen returns the length of packable vector.
func (p *ComplexPacker) PackLen() int {
	return p.pack.PackLen()
}

// Pack packs v.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *ComplexPacker) Pack(v []complex128) *Plaintext {
	ptOut := NewPoly(p.params.Rank(), TypeReal)
	p.PackTo(ptOut, v)
	return ptOut
}

// PackTo packs v to ptOut.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *ComplexPacker) PackTo(ptOut *Plaintext, v []complex128) {
	if p.PackLen()%len(v) != 0 {
		panic("len(v) must divide PackLen")
	} else if ptOut.Rank() != p.params.Rank() {
		panic("inconsistent input(s)")
	}

	p.pack.PackTo(ptOut.Value, v)
}

// UnPack unpacks pt.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *ComplexPacker) UnPack(pt *Plaintext) []complex128 {
	vOut := make([]complex128, p.pack.PackLen())
	p.UnPackTo(vOut, pt)
	return vOut
}

// UnPackTo unpacks pt to vOut.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *ComplexPacker) UnPackTo(vOut []complex128, pt *Plaintext) {
	if p.PackLen()%len(vOut) != 0 {
		panic("len(vOut) must divide PackLen")
	} else if pt.Rank() != p.params.Rank() {
		panic("inconsistent input(s)")
	}

	p.pack.UnPackTo(vOut, pt.Value)
}

// Cube returns the form of the hypercube structure.
func (p *ComplexPacker) Cube() []int {
	return p.pack.Cube()
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *ComplexPacker) CubeGen() []uint64 {
	return p.pack.CubeGen()
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *ComplexPacker) RotIdxToAutIdx(idx []int) int {
	return p.pack.RotIdxToAutIdx(idx)
}
