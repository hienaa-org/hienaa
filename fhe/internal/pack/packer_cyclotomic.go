package pack

import (
	"math/bits"
	"sync"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// pow2CyclotomicMod1Packer is a packer for the power-of-two cyclotomic ring,
// where the modulus is a multiple of prime powers that are 1 mod 4.
type pow2CyclotomicMod1Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// ntt is the transformer for NTT.
	ntt dft.Transformer

	// cube is the form of the hypercube structure.
	cube []int
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	pool *sync.Pool
}

// newPow2CyclotomicMod1Packer creates a new [pow2CyclotomicMod1Packer].
func newPow2CyclotomicMod1Packer(params dft.RingParameters, mod *num.Modulus) *pow2CyclotomicMod1Packer {
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

	return &pow2CyclotomicMod1Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		ntt:     ntt,

		cube:    []int{packLen >> 1, 2},
		cubeGen: []uint64{5, uint64(params.CycloOrder() - 1)},

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, packLen)
				return &v
			},
		},
	}
}

// Params returns the ring parameters.
func (p *pow2CyclotomicMod1Packer) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus used for the packing/unpacking.
func (p *pow2CyclotomicMod1Packer) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *pow2CyclotomicMod1Packer) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *pow2CyclotomicMod1Packer) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *pow2CyclotomicMod1Packer) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) > p.packLen {
		panic("input(s) shape not consistent")
	}

	vBufPtr := p.pool.Get().(*[]uint64)
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	vBufPow5Ptr := p.pool.Get().(*[]uint64)
	vBufPow5 := *vBufPow5Ptr
	defer p.pool.Put(vBufPow5Ptr)

	for i := 0; i < p.packLen/len(v); i++ {
		copy(vBuf[i*len(v):(i+1)*len(v)], v)
	}

	pow5 := 1
	mask := p.packLen<<1 - 1
	for i := 0; i < p.packLen>>1; i++ {
		vBufPow5[pow5>>1] = vBuf[p.packLen>>1-i-1]
		vBufPow5[p.packLen-pow5>>1-1] = vBuf[p.packLen-i-1]
		pow5 = (5 * pow5) & mask
	}

	vec.RadixReverseInPlace(vBufPow5, 2)
	vec.MFormTo(vBufPow5, vBufPow5, p.mod)
	p.ntt.InverseTo(vBufPow5, vBufPow5)

	clear(vPack)
	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		vPack[i*skip] = vBufPow5[i]
	}
}

// UnPack returns the unpacking of vPack.
func (p *pow2CyclotomicMod1Packer) UnPack(vPack []uint64) []uint64 {
	v := make([]uint64, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

// UnPackTo unpacks vPack to v.
func (p *pow2CyclotomicMod1Packer) UnPackTo(v, vPack []uint64) {
	if len(v)%p.packLen != 0 || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	vBufPtr := p.pool.Get().(*[]uint64)
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)
	clear(vBuf)

	vBufPow5Ptr := p.pool.Get().(*[]uint64)
	vBufPow5 := *vBufPow5Ptr
	defer p.pool.Put(vBufPow5Ptr)

	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		vBuf[i] = vPack[i*skip]
	}
	p.ntt.ForwardTo(vBuf, vBuf)
	vec.InvMFormTo(vBuf, vBuf, p.mod)
	vec.RadixReverseInPlace(vBuf, 2)

	pow5 := 1
	mask := p.packLen<<1 - 1
	for i := 0; i < p.packLen>>1; i++ {
		vBufPow5[p.packLen>>1-i-1] = vBuf[pow5>>1]
		vBufPow5[p.packLen-i-1] = vBuf[p.packLen-pow5>>1-1]
		pow5 = (5 * pow5) & mask
	}
	copy(v, vBufPow5[:len(v)])
}

// Cube returns the form of the hypercube structure.
func (p *pow2CyclotomicMod1Packer) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *pow2CyclotomicMod1Packer) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *pow2CyclotomicMod1Packer) RotIdxToAutIdx(idx []int) int {
	switch {
	case p.packLen == 1:
		if len(idx) != 1 {
			panic("input(s) shape not consistent")
		}
		return 1

	case p.packLen == 2:
		if len(idx) != 1 {
			panic("input(s) shape not consistent")
		}
		return int(num.Exp(5, uint64(idx[0]), nil)) & (p.params.CycloOrder() - 1)

	default:
		if len(idx) != 2 {
			panic("input(s) shape not consistent")
		}
		autIdx := int(num.Exp(5, uint64(idx[0]), nil)) & (p.params.CycloOrder() - 1)
		if idx[1]&1 == 1 {
			autIdx = p.params.CycloOrder() - autIdx
		}
		return autIdx
	}
}

// pow2CyclotomicMod3Packer is a packer for the power-of-two cyclotomic ring,
// where the modulus is a multiple of prime powers that are 3 mod 4.
type pow2CyclotomicMod3Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// packIdx is the index mapping for the packing.
	packIdx []int

	// cube is the form of the hypercube structure.
	cube []int
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	// tw is the twiddle factor for NTT.
	tw []gnum.GaussianInt
	// twInv is the twiddle factor for InvNTT.
	twInv []gnum.GaussianInt

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	pool *sync.Pool
}

// newPow2CyclotomicMod3Packer creates a new [pow2CyclotomicMod3Packer].
func newPow2CyclotomicMod3Packer(params dft.RingParameters, mod *num.Modulus) *pow2CyclotomicMod3Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloOrder() >> 2
	logPackLen := bits.TrailingZeros64(primes[0]+1) - 1
	if packLen > (1 << logPackLen) {
		packLen = 1 << logPackLen
	}

	packIdx := make([]int, packLen<<1)
	mask := packLen<<2 - 1
	for i := 0; i < packLen; i++ {
		packIdx[i] = int(num.Exp(5, uint64(packLen<<1-i), nil) & uint64(mask))
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
	vec.RadixReverseInPlace(tw, 2)
	vec.RadixReverseInPlace(twInv, 2)

	rankInv := num.Inv(uint64(packLen<<1), mod)

	return &pow2CyclotomicMod3Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		packIdx: packIdx,

		cube:    []int{packLen},
		cubeGen: []uint64{5},

		tw:    tw,
		twInv: twInv,

		rankInv: rankInv,

		pool: &sync.Pool{
			New: func() any {
				v := make([]gnum.GaussianInt, packLen<<1)
				return &v
			},
		},
	}
}

// Params returns the ring parameters.
func (p *pow2CyclotomicMod3Packer) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *pow2CyclotomicMod3Packer) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *pow2CyclotomicMod3Packer) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *pow2CyclotomicMod3Packer) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *pow2CyclotomicMod3Packer) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) > p.packLen {
		panic("input(s) shape not consistent")
	}

	vBufPtr := p.pool.Get().(*[]gnum.GaussianInt)
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	for i := 0; i < p.packLen; i++ {
		idx0, idx1 := p.packIdx[i], p.packIdx[i+p.packLen]

		vBuf[idx0].Real = v[i%len(v)]
		vBuf[idx0].Imag = 0

		vBuf[idx1].Real = v[i%len(v)]
		vBuf[idx1].Imag = 0
	}

	vec.RadixReverseInPlace(vBuf, 2)
	invNTTGaloisRingInPlacePow2(vBuf, p.twInv, p.mod)

	clear(vPack)
	skip := p.params.Rank() / (p.packLen << 1)
	for i := range vBuf {
		vPack[i*skip] = num.Mul(vBuf[i].Real, p.rankInv, p.mod)
	}
}

// UnPack returns the unpacking of vPack.
func (p *pow2CyclotomicMod3Packer) UnPack(vPack []uint64) []uint64 {
	v := make([]uint64, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

// UnPackTo unpacks vPack to v.
func (p *pow2CyclotomicMod3Packer) UnPackTo(v, vPack []uint64) {
	if len(v)%p.packLen != 0 || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	vBufPtr := p.pool.Get().(*[]gnum.GaussianInt)
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)
	clear(vBuf)

	skip := p.params.Rank() / (p.packLen << 1)
	for i := range vBuf {
		vBuf[i].Real = vPack[i*skip]
		vBuf[i].Imag = 0
	}

	nttGaloisRingInPlacePow2(vBuf, p.tw, p.mod)
	vec.RadixReverseInPlace(vBuf, 2)

	for i := range v {
		v[i] = vBuf[p.packIdx[i]].Real
	}
}

// Cube returns the form of the hypercube structure.
func (p *pow2CyclotomicMod3Packer) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *pow2CyclotomicMod3Packer) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *pow2CyclotomicMod3Packer) RotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("input(s) shape not consistent")
	}
	return int(num.Exp(5, uint64(idx[0]), nil)) & (p.params.CycloOrder() - 1)
}

// anyCyclotomicPacker is a packer for arbitrary cyclotomic rings,
// where the modulus is NTT-friendly.
type anyCyclotomicPacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// ntt is the transformer for NTT.
	ntt dft.Transformer

	// cube is the form of the hypercube structure.
	cube []int
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	// cycloOrdMod is the cyclotomic order modulus.
	cycloOrdMod *num.Modulus
}

// newAnyCyclotomicPacker creates a new [anyCyclotomicPacker].
func newAnyCyclotomicPacker(params dft.RingParameters, mod *num.Modulus) *anyCyclotomicPacker {
	packLen := int(params.Rank())
	ntt := dft.NewTransformer(params, mod)

	cycloOrdMod := num.NewModulus(params.CycloOrder())
	primes, exps := num.Factor(cycloOrdMod.Value())
	cubeGen := num.GeneratorsWithFactors(cycloOrdMod, primes, exps)
	cube := make([]int, len(cubeGen))
	if primes[0] == 2 && exps[0] == 1 {
		for i := range cube {
			pExp := num.Exp(primes[i+1], exps[i+1], nil)
			cube[i] = int(pExp - pExp/primes[i])
			cubeGen[i] = num.Inv(cubeGen[i], cycloOrdMod)
		}
	} else {
		for i := range cube {
			pExp := num.Exp(primes[i], exps[i], nil)
			cube[i] = int(pExp - pExp/primes[i])
			cubeGen[i] = num.Inv(cubeGen[i], cycloOrdMod)
		}
	}

	return &anyCyclotomicPacker{
		params: params,
		mod:    mod,

		packLen: packLen,
		ntt:     ntt,

		cube:    cube,
		cubeGen: cubeGen,

		cycloOrdMod: cycloOrdMod,
	}
}

// Params returns the ring parameters.
func (p *anyCyclotomicPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus used for the packing/unpacking.
func (p *anyCyclotomicPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *anyCyclotomicPacker) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *anyCyclotomicPacker) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *anyCyclotomicPacker) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) != p.packLen {
		panic("input(s) shape not consistent")
	}

	vec.MFormTo(vPack, v, p.mod)
	p.ntt.InverseTo(vPack, vPack)
}

// UnPack returns the unpacking of vPack.
func (p *anyCyclotomicPacker) UnPack(v []uint64) []uint64 {
	vPack := make([]uint64, p.packLen)
	p.UnPackTo(vPack, v)
	return vPack
}

// UnPackTo unpacks vPack to v.
func (p *anyCyclotomicPacker) UnPackTo(v, vPack []uint64) {
	if len(v) != p.packLen || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	p.ntt.ForwardTo(v, vPack)
	vec.InvMFormTo(v, v, p.mod)
}

// Cube returns the form of the hypercube structure.
func (p *anyCyclotomicPacker) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *anyCyclotomicPacker) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *anyCyclotomicPacker) RotIdxToAutIdx(idx []int) int {
	if len(idx) != len(p.cubeGen) {
		panic("input(s) shape not consistent")
	}
	autIdx := uint64(1)
	for i := range p.cubeGen {
		autIdx = num.Mul(autIdx, num.Exp(p.cubeGen[i], uint64(idx[i]), p.cycloOrdMod), p.cycloOrdMod)
	}
	return int(autIdx)
}

// trivialPacker is a packer for arbitrary cyclotomic rings,
// where the modulus is not NTT-friendly.
type trivialPacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int

	// cube is the form of the hypercube structure.
	cube []int
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64
}

// newTrivialPacker creates a new [trivialPacker].
func newTrivialPacker(params dft.RingParameters, mod *num.Modulus) *trivialPacker {
	return &trivialPacker{
		params: params,
		mod:    mod,

		packLen: 1,

		cube:    []int{1},
		cubeGen: []uint64{1},
	}
}

// Params returns the ring parameters.
func (p *trivialPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus used for the packing/unpacking.
func (p *trivialPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *trivialPacker) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *trivialPacker) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *trivialPacker) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) != p.packLen {
		panic("input(s) shape not consistent")
	}

	clear(vPack)
	vPack[0] = v[0]
}

// UnPack returns the unpacking of vPack.
func (p *trivialPacker) UnPack(vPack []uint64) []uint64 {
	v := make([]uint64, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

// UnPackTo unpacks vPack to v.
func (p *trivialPacker) UnPackTo(v, vPack []uint64) {
	if len(v) != p.packLen || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	v[0] = vPack[0]
}

// Cube returns the form of the hypercube structure.
func (p *trivialPacker) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *trivialPacker) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *trivialPacker) RotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("input(s) shape not consistent")
	}
	return 1
}
