package bgv

import (
	pack "github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type packerBuffer struct {
	poly *crt.Element
}

func newPackerBuffer(rank int, modLen int) packerBuffer {
	return packerBuffer{poly: crt.NewPoly(rank, modLen)}
}

// Packer packs a vector of uint64 into [*Plaintext].
type Packer struct {
	params Parameters
	pack   pack.PackerInt
	eval   crt.Operator
	emb    *crt.Embedder
	buf    packerBuffer
}

func NewPacker(params Parameters) *Packer {
	pack := pack.NewPackerInt(params.ringParams, params.messageModulus)
	eval := crt.NewOperator(params.ringParams, params.modulus)
	emb := crt.NewEmbedder([]*num.Modulus{params.messageModulus}, params.modulus)
	buf := newPackerBuffer(params.ringParams.Rank(), len(params.modulus))
	return &Packer{
		params: params,
		pack:   pack,
		eval:   eval,
		emb:    emb,
		buf:    buf,
	}
}

// PackLen returns the length of packable vector.
func (p *Packer) PackLen() int {
	return p.pack.PackLen()
}

// SafeCopy returns a thread-safe copy.
func (p *Packer) SafeCopy() *Packer {
	return &Packer{
		params: p.params,
		pack:   p.pack.SafeCopy(),
		eval:   p.eval.SafeCopy(),
		emb:    p.emb.SafeCopy(),
		buf:    newPackerBuffer(p.params.ringParams.Rank(), len(p.params.modulus)),
	}
}

// Pack packs v.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) Pack(v []uint64) *Plaintext {
	pOut := Plaintext{Value: p.eval.NewPoly()}
	p.PackTo(&pOut, v)
	return &pOut
}

// PackTo packs v to ptOut.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) PackTo(ptOut *Plaintext, v []uint64) {
	if p.PackLen()%len(v) != 0 {
		panic("PackTo: message length should divide the maximum packing length")
	}
	if ptOut.Value.Rank() != p.eval.Params().Rank() {
		panic("PackTo: output rank should match the parameters of the packer")
	}
	if ptOut.Value.ModLen() > len(p.eval.Modulus()) {
		panic("PackTo: output modulus length is larger than the number of moduli")
	}

	p.pack.PackTo(p.buf.poly.Coeffs[0], v)

	modulus := p.eval.Modulus()
	for i := 0; i < ptOut.Value.ModLen(); i++ {
		vec.ReduceTo(ptOut.Value.Coeffs[i], p.buf.poly.Coeffs[0], modulus[i])
	}

	ptOut.Value.IsNTT = false
}

// UnPack unpacks pt.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPack(pt *Plaintext) []uint64 {
	vOut := make([]uint64, p.pack.PackLen())
	p.UnPackTo(vOut, pt)
	return vOut
}

// UnPack unpacks pt to vOut.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPackTo(vOut []uint64, pt *Plaintext) {
	if p.PackLen()%len(vOut) != 0 {
		panic("UnPackTo: message length should divide the maximum packing length")
	}
	if pt.Value.Rank() != p.eval.Params().Rank() {
		panic("UnPackTo: input rank should match the parameters of the packer")
	}
	if pt.Value.ModLen() > len(p.eval.Modulus()) {
		panic("UnPackTo: input modulus length is larger than the number of moduli")
	}

	ptLen := len(pt.Value.Coeffs)
	copy(p.buf.poly.Coeffs[:ptLen], pt.Value.Coeffs)

	if pt.Value.IsNTT {
		p.eval.InvNTTTo(p.buf.poly, p.buf.poly)
	}

	p.emb.EmbedVecTo(p.buf.poly.Coeffs[:1], p.buf.poly.Coeffs[:ptLen])

	p.pack.UnPackTo(vOut, p.buf.poly.Coeffs[0])
}
