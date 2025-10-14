package pack

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclotomicPow2NTTPacker is a packer for cyclotomic NTT with power-of-two cyclotomic order.
type cyclotomicPow2NTTPacker struct {
	params  dft.RingParameters
	mod     *num.Modulus
	packLen int

	ntt dft.Transformer

	buf packerBuffer
}

// newCyclotomicPow2NTTPacker creates a new [cyclotomicPow2NTTPacker].
func newCyclotomicPow2NTTPacker(params dft.RingParameters, mod *num.Modulus) *cyclotomicPow2NTTPacker {
	primes, _ := num.Factor(mod.Value())
	packLen := int(params.Rank())
	for i := range primes {
		ord := num.Order(primes[i], num.NewModulus(params.CycloOrder()))
		if packLen > int(params.Rank())/int(ord) {
			packLen = int(params.Rank()) / int(ord)
		}
	}
	nttParams := dft.NewCyclotomicParameters(packLen << 1)
	ntt := dft.NewTransformer(nttParams, mod)

	return &cyclotomicPow2NTTPacker{
		params:  params,
		mod:     mod,
		packLen: packLen,

		ntt: ntt,

		buf: newPackerBuffer(2, packLen),
	}
}

// Params returns the parameters of the packer.
func (p *cyclotomicPow2NTTPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *cyclotomicPow2NTTPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *cyclotomicPow2NTTPacker) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *cyclotomicPow2NTTPacker) SafeCopy() PackerInt {
	return &cyclotomicPow2NTTPacker{
		params:  p.params,
		mod:     p.mod,
		packLen: p.packLen,

		ntt: p.ntt.SafeCopy(),

		buf: newPackerBuffer(2, p.packLen),
	}
}

// Pack packs the input vector into a polynomial.
func (p *cyclotomicPow2NTTPacker) Pack(vIn []uint64) *crt.Poly {
	pOut := crt.NewPoly(p.params.Rank(), 1)
	p.PackTo(pOut, vIn)
	return pOut
}

// PackTo packs the input vector into a polynomial.
func (p *cyclotomicPow2NTTPacker) PackTo(pOut *crt.Poly, vIn []uint64) {
	vLen := len(vIn)

	if p.packLen%vLen != 0 {
		panic("packTo: message length should divide the maximum packing length")
	}
	if pOut.Rank() != p.params.Rank() {
		panic("packTo: output rank should match the parameters of the packer")
	}
	if pOut.ModLen() != 1 {
		panic("packTo: output modulus length should be 1")
	}

	for i := 0; i < p.packLen/vLen; i++ {
		copy(p.buf.coeffs[0][vLen*i:vLen*(i+1)], vIn)
	}

	idx := 1
	mask := p.packLen<<1 - 1
	for i := 0; i < p.packLen>>1; i++ {
		pOut.Coeffs[0][idx>>1] = vIn[p.packLen>>1-i-1]
		pOut.Coeffs[0][p.packLen-idx>>1-1] = vIn[p.packLen-i-1]
		idx = (5 * idx) & mask
	}

	copy(p.buf.coeffs[0], pOut.Coeffs[0][:p.packLen])
	vec.RadixReverseInPlace(p.buf.coeffs[0], 2)
	p.ntt.InverseTo(p.buf.coeffs[0], p.buf.coeffs[0])

	clear(pOut.Coeffs[0])
	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		pOut.Coeffs[0][i*skip] = p.buf.coeffs[0][i]
	}
}

// UnPack unpacks the polynomial into a vector.
func (p *cyclotomicPow2NTTPacker) UnPack(pIn *crt.Poly) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, pIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *cyclotomicPow2NTTPacker) UnPackTo(vOut []uint64, pIn *crt.Poly) {
	vLen := len(vOut)

	if p.packLen%vLen != 0 {
		panic("packTo: message length should divide the maximum packing length")
	}
	if pIn.Rank() != p.params.Rank() {
		panic("packTo: input rank should match the parameters of the packer")
	}
	if pIn.ModLen() != 1 {
		panic("packTo: input modulus length should be 1")
	}

	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		p.buf.coeffs[0][i] = pIn.Coeffs[0][i*skip]
	}
	p.ntt.ForwardTo(p.buf.coeffs[0], p.buf.coeffs[0])
	vec.RadixReverseInPlace(p.buf.coeffs[0], 2)

	idx := 1
	mask := p.packLen<<1 - 1
	for i := 0; i < p.packLen>>1; i++ {
		p.buf.coeffs[1][p.packLen>>1-i-1] = p.buf.coeffs[0][idx>>1]
		p.buf.coeffs[1][p.packLen-i-1] = p.buf.coeffs[0][p.packLen-idx>>1-1]
		idx = (5 * idx) & mask
	}
	copy(vOut, p.buf.coeffs[1][:vLen])
}

type cyclotomicAnyNTTPacker struct {
	params  dft.RingParameters
	mod     *num.Modulus
	packLen int

	ntt dft.Transformer
}

func newCyclotomicAnyNTTPacker(params dft.RingParameters, mod *num.Modulus) *cyclotomicAnyNTTPacker {
	packLen := int(params.Rank())
	ntt := dft.NewTransformer(params, mod)

	return &cyclotomicAnyNTTPacker{
		params:  params,
		mod:     mod,
		packLen: packLen,

		ntt: ntt,
	}
}

// Params returns the parameters of the packer.
func (p *cyclotomicAnyNTTPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *cyclotomicAnyNTTPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *cyclotomicAnyNTTPacker) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *cyclotomicAnyNTTPacker) SafeCopy() PackerInt {
	return &cyclotomicAnyNTTPacker{
		params:  p.params,
		mod:     p.mod,
		packLen: p.packLen,

		ntt: p.ntt.SafeCopy(),
	}
}

// Pack packs the input vector into a polynomial.
func (p *cyclotomicAnyNTTPacker) Pack(vIn []uint64) *crt.Poly {
	pOut := crt.NewPoly(p.params.Rank(), 1)
	p.PackTo(pOut, vIn)
	return pOut
}

// PackTo packs the input vector into a polynomial.
func (p *cyclotomicAnyNTTPacker) PackTo(pOut *crt.Poly, vIn []uint64) {
	vLen := len(vIn)

	if p.packLen != vLen {
		panic("packTo: message length should match the packing length")
	}
	if pOut.Rank() != p.params.Rank() {
		panic("packTo: output rank should match the parameters of the packer")
	}
	if pOut.ModLen() != 1 {
		panic("packTo: output modulus length should be 1")
	}

	p.ntt.InverseTo(pOut.Coeffs[0], vIn)
}

// UnPack unpacks the polynomial into a vector.
func (p *cyclotomicAnyNTTPacker) UnPack(pIn *crt.Poly) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, pIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *cyclotomicAnyNTTPacker) UnPackTo(vOut []uint64, pIn *crt.Poly) {
	vLen := len(vOut)

	if p.packLen != vLen {
		panic("packTo: message length should match the packing length")
	}
	if pIn.Rank() != p.params.Rank() {
		panic("packTo: input rank should match the parameters of the packer")
	}
	if pIn.ModLen() != 1 {
		panic("packTo: input modulus length should be 1")
	}

	p.ntt.ForwardTo(vOut, pIn.Coeffs[0])
}
