package ckks

import "github.com/hienaa-org/hienaa/fhe/internal/pack"

// Packer packs a vector of complex numbers into [*Plaintext].
type Packer struct {
	params Parameters
	packer pack.RealPacker
}

// NewPacker creates a new [Packer].
func NewPacker(params Parameters) *Packer {
	return &Packer{
		params: params,
		packer: pack.NewRealPacker(params.ringParams),
	}
}

func (p *Packer) Params() Parameters {
	return p.params
}

func (p *Packer) PackLen() int {
	return p.packer.PackLen()
}

// func (p *Packer) Pack(v []complex128) *Plaintext {
// 	ptOut := NewPlaintext(p.params.ringParams.Rank())
// 	p.PackTo(ptOut, v)
// 	return ptOut
// }
