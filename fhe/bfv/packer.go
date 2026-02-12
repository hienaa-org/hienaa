package bfv

import (
	"sync"

	pack "github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Packer packs a vector of uint64 into [*Plaintext].
type Packer struct {
	params Parameters
	packer pack.IntPacker

	crtOp crt.Operator

	embModToPlainMod *crt.Embedder
	embPlainModToMod *crt.Embedder

	pool *sync.Pool
}

func NewPacker(params Parameters) *Packer {
	return &Packer{
		params: params,
		packer: pack.NewIntPacker(params.ringParams, params.messageModulus),

		crtOp: crt.NewOperator(params.ringParams, params.modulus),

		embModToPlainMod: crt.NewEmbedder([]*num.Modulus{params.messageModulus}, params.modulus),
		embPlainModToMod: crt.NewEmbedder(params.modulus, []*num.Modulus{params.messageModulus}),

		pool: &sync.Pool{
			New: func() any {
				return crt.NewPoly(params.ringParams.Rank(), len(params.modulus))
			},
		},
	}
}

// PackLen returns the length of packable vector.
func (p *Packer) PackLen() int {
	return p.packer.PackLen()
}

// Pack packs v.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) Pack(v []uint64) *Plaintext {
	pOut := Plaintext{Value: p.crtOp.NewPoly()}
	p.PackTo(&pOut, v)
	return &pOut
}

// PackTo packs v to ptOut.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) PackTo(ptOut *Plaintext, v []uint64) {
	if p.PackLen()%len(v) != 0 {
		panic("len(v) must divide PackLen")
	} else if ptOut.Value.Rank() != p.crtOp.Params().Rank() || ptOut.Value.ModLen() > len(p.crtOp.Modulus()) {
		panic("inconsistent input(s)")
	}

	ptBuf := p.pool.Get().(*crt.Element)
	defer p.pool.Put(ptBuf)

	ptBuf = ptBuf.WithModIdx(0)
	p.packer.PackTo(ptBuf.Coeffs[0], v)
	p.embPlainModToMod.EmbedTo(ptOut.Value, ptBuf)
}

// UnPack unpacks pt.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPack(pt *Plaintext) []uint64 {
	vOut := make([]uint64, p.packer.PackLen())
	p.UnPackTo(vOut, pt)
	return vOut
}

// UnPack unpacks pt to vOut.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPackTo(vOut []uint64, pt *Plaintext) {
	if p.PackLen()%len(vOut) != 0 {
		panic("len(vOut) must divide PackLen")
	} else if pt.Value.Rank() != p.crtOp.Params().Rank() || pt.Value.ModLen() > len(p.crtOp.Modulus()) {
		panic("inconsistent input(s)")
	}

	ptBuf := p.pool.Get().(*crt.Element)
	defer p.pool.Put(ptBuf)

	ptBuf = ptBuf.WithModIdx(vec.Range(0, pt.Value.ModLen())...)
	ptBuf.CopyFrom(pt.Value)
	if ptBuf.IsNTT {
		p.crtOp.SubOperator(vec.Range(0, pt.Value.ModLen())...).InvNTTTo(ptBuf, ptBuf)
	}
	p.embModToPlainMod.EmbedVecTo(ptBuf.Coeffs[:1], ptBuf.Coeffs)
	p.packer.UnPackTo(vOut, ptBuf.Coeffs[0])
}
