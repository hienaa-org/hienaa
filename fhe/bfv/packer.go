package bfv

import (
	pack "github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

// Packer packs a vector of uint64 into [Plaintext].
type Packer struct {
	params rlwe.Parameters
	packer pack.IntPacker
}

// NewPacker creates a new [Packer].
func NewPacker(params rlwe.Parameters, msgMod *num.Modulus) *Packer {
	return &Packer{
		params: params,
		packer: pack.NewIntPacker(params.RingParams(), msgMod),
	}
}

// Params returns the parameters.
func (p *Packer) Params() rlwe.Parameters {
	return p.params
}

// PackLen returns the length of packable vector.
func (p *Packer) PackLen() int {
	return p.packer.PackLen()
}

// Pack packs v.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) Pack(v []uint64) Plaintext {
	ptOut := NewPoly(p.params.Rank())
	p.PackTo(ptOut, v)
	return ptOut
}

// PackTo packs v to ptOut.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) PackTo(ptOut Plaintext, v []uint64) {
	if p.PackLen()%len(v) != 0 {
		panic("len(v) must divide PackLen")
	} else if ptOut.Rank() != p.params.Rank() {
		panic("inconsistent input(s)")
	}

	p.packer.PackTo(ptOut, v)
}

// UnPack unpacks pt.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPack(pt Plaintext) []uint64 {
	vOut := make([]uint64, p.packer.PackLen())
	p.UnPackTo(vOut, pt)
	return vOut
}

// UnPack unpacks pt to vOut.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPackTo(vOut []uint64, pt Plaintext) {
	if p.PackLen()%len(vOut) != 0 {
		panic("len(vOut) must divide PackLen")
	} else if pt.Rank() != p.params.Rank() {
		panic("inconsistent input(s)")
	}

	p.packer.UnPackTo(vOut, pt)
}

// Cube returns the form of the hypercube structure.
func (p *Packer) Cube() []int {
	return p.packer.Cube()
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *Packer) CubeGen() []uint64 {
	return p.packer.CubeGen()
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *Packer) RotIdxToAutIdx(idx []int) int {
	return p.packer.RotIdxToAutIdx(idx)
}
