package pack

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclotomicPow2Mod1Packer is a packer for the power-of-two cyclotomic ring,
// where the modulus is a multiple of prime powers that are 1 mod 4.
type cyclotomicPow2Mod1Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// ntt is the transformer for NTT.
	ntt dft.Transformer

	buf packerBuffer
}

// newCyclotomicPow2Mod1Packer creates a new [cyclotomicPow2Mod1Packer].
func newCyclotomicPow2Mod1Packer(params dft.RingParameters, mod *num.Modulus) *cyclotomicPow2Mod1Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloOrder() >> 1
	for i := range primes {
		ithLogPackLen := bits.TrailingZeros64(primes[i]-1) - 1
		if packLen > (1 << ithLogPackLen) {
			packLen = 1 << ithLogPackLen
		}
	}
	nttParams := dft.NewCyclotomicParameters(packLen << 1)
	ntt := dft.NewTransformer(nttParams, mod)

	return &cyclotomicPow2Mod1Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		ntt:     ntt,

		buf: newPackerBuffer(2, packLen),
	}
}

// Params returns the parameters of the packer.
func (p *cyclotomicPow2Mod1Packer) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *cyclotomicPow2Mod1Packer) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *cyclotomicPow2Mod1Packer) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *cyclotomicPow2Mod1Packer) SafeCopy() PackerInt {
	return &cyclotomicPow2Mod1Packer{
		params: p.params,
		mod:    p.mod,

		packLen: p.packLen,
		ntt:     p.ntt.SafeCopy(),

		buf: newPackerBuffer(2, p.packLen),
	}
}

// Pack packs the input vector into a polynomial.
func (p *cyclotomicPow2Mod1Packer) Pack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.params.Rank())
	p.PackTo(vOut, vIn)
	return vOut
}

// PackTo packs the input vector into a polynomial.
func (p *cyclotomicPow2Mod1Packer) PackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vIn)

	for i := 0; i < p.packLen/vLen; i++ {
		copy(p.buf.coeffs[0][vLen*i:vLen*(i+1)], vIn)
	}

	pow5 := 1
	mask := p.packLen<<1 - 1
	for i := 0; i < p.packLen>>1; i++ {
		p.buf.coeffs[1][pow5>>1] = p.buf.coeffs[0][p.packLen>>1-i-1]
		p.buf.coeffs[1][p.packLen-pow5>>1-1] = p.buf.coeffs[0][p.packLen-i-1]
		pow5 = (5 * pow5) & mask
	}

	copy(p.buf.coeffs[0], p.buf.coeffs[1])
	vec.RadixReverseInPlace(p.buf.coeffs[0], 2)
	vec.MFormTo(p.buf.coeffs[0], p.buf.coeffs[0], p.mod)
	p.ntt.InverseTo(p.buf.coeffs[0], p.buf.coeffs[0])

	clear(vOut)
	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		vOut[i*skip] = p.buf.coeffs[0][i]
	}
}

// UnPack unpacks the polynomial into a vector.
func (p *cyclotomicPow2Mod1Packer) UnPack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, vIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *cyclotomicPow2Mod1Packer) UnPackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vOut)

	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		p.buf.coeffs[0][i] = vIn[i*skip]
	}
	p.ntt.ForwardTo(p.buf.coeffs[0], p.buf.coeffs[0])
	vec.InvMFormTo(p.buf.coeffs[0], p.buf.coeffs[0], p.mod)
	vec.RadixReverseInPlace(p.buf.coeffs[0], 2)

	pow5 := 1
	mask := p.packLen<<1 - 1
	for i := 0; i < p.packLen>>1; i++ {
		p.buf.coeffs[1][p.packLen>>1-i-1] = p.buf.coeffs[0][pow5>>1]
		p.buf.coeffs[1][p.packLen-i-1] = p.buf.coeffs[0][p.packLen-pow5>>1-1]
		pow5 = (5 * pow5) & mask
	}
	copy(vOut, p.buf.coeffs[1][:vLen])
}

// cyclotomicPow2Mod3Packer is a packer for the power-of-two cyclotomic ring,
// where the modulus is a multiple of prime powers that are 3 mod 4.
type cyclotomicPow2Mod3Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// packIdx is the index mapping for the packing.
	packIdx []int

	// tw is the twiddle factor for NTT.
	tw []gnum.GaussianInt
	// twInv is the twiddle factor for InvNTT.
	twInv []gnum.GaussianInt

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	buf pow2Mod3PackerBuffer
}

// newCyclotomicPow2Mod3Packer creates a new [cyclotomicPow2Mod3Packer].
func newCyclotomicPow2Mod3Packer(params dft.RingParameters, mod *num.Modulus) *cyclotomicPow2Mod3Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloOrder() >> 2
	LogPackLen := bits.TrailingZeros64(primes[0]+1) - 1
	if packLen > (1 << LogPackLen) {
		packLen = 1 << LogPackLen
	}

	packIdx := make([]int, packLen<<1)
	mask := packLen<<2 - 1
	for i := 0; i < packLen; i++ {
		packIdx[i] = int(num.Exp(5, uint64(i), nil) & uint64(mask))
		packIdx[i+packLen] = (packIdx[i] * int(primes[0])) & mask

		packIdx[i] >>= 1
		packIdx[i+packLen] >>= 1
	}

	root := gIntNthRoot(packLen<<2, mod)
	tw := make([]gnum.GaussianInt, packLen<<1)
	twInv := make([]gnum.GaussianInt, packLen<<1)
	tw[0], tw[1] = gnum.GaussianInt{Real: 1, Imag: 0}, root
	twInv[0], twInv[1] = gnum.GaussianInt{Real: 1, Imag: 0}, gnum.Inv(root, mod)
	for i := 2; i < 2*packLen; i++ {
		tw[i] = gnum.Mul(tw[i-1], tw[1], mod)
		twInv[i] = gnum.Mul(twInv[i-1], twInv[1], mod)
	}
	bitReverseInPlace(tw)
	bitReverseInPlace(twInv)

	rankInv := num.Inv(uint64(packLen<<1), mod)

	return &cyclotomicPow2Mod3Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		packIdx: packIdx,

		tw:    tw,
		twInv: twInv,

		rankInv: rankInv,

		buf: newPow2Mod3PackerBuffer(packLen << 1),
	}
}

// Params returns the parameters of the packer.
func (p *cyclotomicPow2Mod3Packer) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *cyclotomicPow2Mod3Packer) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *cyclotomicPow2Mod3Packer) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *cyclotomicPow2Mod3Packer) SafeCopy() PackerInt {
	return &cyclotomicPow2Mod3Packer{
		params: p.params,
		mod:    p.mod,

		packLen: p.packLen,
		packIdx: p.packIdx,

		tw:    p.tw,
		twInv: p.twInv,

		rankInv: p.rankInv,

		buf: newPow2Mod3PackerBuffer(p.packLen << 1),
	}
}

// Pack packs the input vector into a polynomial.
func (p *cyclotomicPow2Mod3Packer) Pack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.params.Rank())
	p.PackTo(vOut, vIn)
	return vOut
}

// PackTo packs the input vector into a polynomial.
func (p *cyclotomicPow2Mod3Packer) PackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vIn)

	for i := 0; i < p.packLen; i++ {
		idx1 := p.packIdx[i]
		idx2 := p.packIdx[i+p.packLen]

		p.buf.coeffs[idx1].Real = vIn[i&(vLen-1)]
		p.buf.coeffs[idx1].Imag = 0

		p.buf.coeffs[idx2].Real = vIn[i&(vLen-1)]
		p.buf.coeffs[idx2].Imag = 0
	}

	bitReverseInPlace(p.buf.coeffs)
	invNTTGaloisRingInPlacePow2(p.buf.coeffs, p.twInv, p.mod)

	clear(vOut)
	skip := p.params.Rank() / (p.packLen << 1)
	for i := range p.buf.coeffs {
		vOut[i*skip] = num.Mul(p.buf.coeffs[i].Real, p.rankInv, p.mod)
	}
}

// UnPack unpacks the polynomial into a vector.
func (p *cyclotomicPow2Mod3Packer) UnPack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, vIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *cyclotomicPow2Mod3Packer) UnPackTo(vOut []uint64, vIn []uint64) {
	skip := p.params.Rank() / (p.packLen << 1)
	for i := range p.buf.coeffs {
		p.buf.coeffs[i].Real = vIn[i*skip]
		p.buf.coeffs[i].Imag = 0
	}

	nttGaloisRingInPlacePow2(p.buf.coeffs, p.tw, p.mod)
	bitReverseInPlace(p.buf.coeffs)

	for i := range vOut {
		idx := p.packIdx[i]
		vOut[i] = p.buf.coeffs[idx].Real
	}
}

// cyclotomicAnyNTTPacker is a packer for the cyclotomic ring,
// where the modulus is a multiple of prime powers that are 1 mod 4.
type cyclotomicAnyNTTPacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// ntt is the transformer for NTT.
	ntt dft.Transformer
}

// newCyclotomicAnyNTTPacker creates a new [cyclotomicAnyNTTPacker].
func newCyclotomicAnyNTTPacker(params dft.RingParameters, mod *num.Modulus) *cyclotomicAnyNTTPacker {
	packLen := int(params.Rank())
	ntt := dft.NewTransformer(params, mod)

	return &cyclotomicAnyNTTPacker{
		params: params,
		mod:    mod,

		packLen: packLen,
		ntt:     ntt,
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
func (p *cyclotomicAnyNTTPacker) Pack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.params.Rank())
	p.PackTo(vOut, vIn)
	return vOut
}

// PackTo packs the input vector into a polynomial.
func (p *cyclotomicAnyNTTPacker) PackTo(vOut []uint64, vIn []uint64) {
	vec.MFormTo(vOut, vIn, p.mod)
	p.ntt.InverseTo(vOut, vOut)
}

// UnPack unpacks the polynomial into a vector.
func (p *cyclotomicAnyNTTPacker) UnPack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, vIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *cyclotomicAnyNTTPacker) UnPackTo(vOut []uint64, vIn []uint64) {
	p.ntt.ForwardTo(vOut, vIn)
	vec.InvMFormTo(vOut, vOut, p.mod)
}
