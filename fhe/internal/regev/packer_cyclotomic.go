package regev

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type pow2Cyclotomic1Mod4Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// packIdx is the index mapping for the packing.
	packIdx []int
	// ntt is the transformer for NTT.
	ntt *dft.Transformer

	// cube is the form of the hypercube structure.
	cube []int
	// cubeMod is [num.Modulus] form of cube.
	// Entries are nil for dimensions of size one.
	cubeMod []*num.Modulus
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	pool *pool.Pool[*[]uint64]
}

// newPow2Cyclotomic1Mod4Packer creates a new [pow2Cyclotomic1Mod4Packer].
func newPow2Cyclotomic1Mod4Packer(params dft.RingParameters, mod *num.Modulus) *pow2Cyclotomic1Mod4Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloIdx() >> 1
	for i := range primes {
		logPackLen := bits.TrailingZeros64(primes[i]-1) - 1
		packLen = int(min(uint64(packLen), uint64(1)<<logPackLen))
	}

	packLen64 := uint64(packLen)
	packLenDoubleMod := num.NewModulus(packLen64 << 1)
	packIdx := make([]int, packLen)
	revShiftBits := 64 - uint64(num.Log2(packLen))
	for i := uint64(0); i < packLen64>>1; i++ {
		pow5 := num.Exp(5, packLen64>>1-i-1, packLenDoubleMod)
		packIdx[i] = int(bits.Reverse64(pow5>>1) >> revShiftBits)
		packIdx[i+packLen64>>1] = packLen - 1 - packIdx[i]
	}

	var cubeMod *num.Modulus
	if packLen>>1 > 1 {
		cubeMod = num.NewModulus(packLen >> 1)
	}

	return &pow2Cyclotomic1Mod4Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		packIdx: packIdx,
		ntt:     dft.NewTransformer(dft.NewCyclotomicParameters(packLen<<1), mod),

		cube:    []int{packLen >> 1, 2},
		cubeMod: []*num.Modulus{cubeMod, num.NewModulus(2)},
		cubeGen: []uint64{5, uint64(params.CycloIdx() - 1)},

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, packLen)
			return &v
		}),
	}
}

func (pck *pow2Cyclotomic1Mod4Packer) packTo(vPack, v []uint64) {
	if len(vPack) != pck.params.Rank() || len(v) > pck.packLen {
		panic("inconsistent inputs")
	}

	if len(v) == 0 {
		clear(vPack)
		return
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

	for i := 0; i < pck.packLen; i++ {
		vBuf[pck.packIdx[i]] = v[i%len(v)]
	}

	vec.MFormTo(vBuf, vBuf, pck.mod)
	pck.ntt.InverseTo(vBuf, vBuf)

	clear(vPack)
	skip := pck.params.Rank() / pck.packLen
	for i := 0; i < pck.params.Rank()/skip; i++ {
		vPack[i*skip] = vBuf[i]
	}
}

func (pck *pow2Cyclotomic1Mod4Packer) unpackTo(v, vPack []uint64) {
	if len(v) == 0 || pck.packLen%len(v) != 0 || len(vPack) != pck.params.Rank() {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)
	clear(vBuf)

	skip := pck.params.Rank() / pck.packLen
	for i := 0; i < pck.params.Rank()/skip; i++ {
		vBuf[i] = vPack[i*skip]
	}
	pck.ntt.ForwardTo(vBuf, vBuf)
	vec.InvMFormTo(vBuf, vBuf, pck.mod)

	for i := range v {
		v[i] = vBuf[pck.packIdx[i]]
	}
}

func (pck *pow2Cyclotomic1Mod4Packer) rotIdxToAutIdx(idx []int) int {
	if len(idx) != len(pck.cube) {
		panic("inconsistent inputs")
	}

	autIdx := 1
	if pck.cubeMod[0] != nil {
		autIdx = int(num.Exp(5, num.Reduce(idx[0], pck.cubeMod[0]), pck.params.CycloIdxMod()))
	}
	if idx[1]&1 == 1 {
		autIdx = pck.params.CycloIdx() - autIdx
	}
	return autIdx
}

type pow2Cyclotomic3Mod4Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// packIdx is the index mapping for the packing in bit-reversed NTT order.
	packIdx []int

	// cube is the form of the hypercube structure.
	cube []int
	// cubeMod is [num.Modulus] form of cube.
	// Entries are nil for dimensions of size one.
	cubeMod []*num.Modulus
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	// tw is the twiddle factor for NTT.
	tw []gnum.GaussianInt
	// twS is the Shoup form of tw.
	twS []gnum.GaussianInt
	// twInv is the twiddle factor for InvNTT.
	twInv []gnum.GaussianInt
	// twInvS is the Shoup form of twInv.
	twInvS []gnum.GaussianInt

	// rankInv is the modular inverse of the rank.
	rankInv uint64
	// rankInvS is the Shoup form of rankInv.
	rankInvS uint64

	pool *pool.Pool[*[]gnum.GaussianInt]
}

// newPow2Cyclotomic3Mod4Packer creates a new [pow2Cyclotomic3Mod4Packer].
func newPow2Cyclotomic3Mod4Packer(params dft.RingParameters, mod *num.Modulus) *pow2Cyclotomic3Mod4Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloIdx() >> 2
	logPackLen := bits.TrailingZeros64(primes[0]+1) - 1
	packLen = int(min(uint64(packLen), uint64(1)<<logPackLen))

	packLen64 := uint64(packLen)
	packLenQuadMod := num.NewModulus(packLen64 << 2)
	packIdx := make([]int, packLen<<1)
	revShiftBits := 64 - uint64(num.Log2(packLen<<1))
	for i := uint64(0); i < packLen64; i++ {
		idx0 := num.Exp(5, packLen64<<1-i, packLenQuadMod)
		idx1 := num.Mul(idx0, primes[0], packLenQuadMod)
		packIdx[i] = int(bits.Reverse64(idx0>>1) >> revShiftBits)
		packIdx[i+packLen64] = int(bits.Reverse64(idx1>>1) >> revShiftBits)
	}

	root := gaussIntNthRoot(packLen<<2, mod)
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

	twS := make([]gnum.GaussianInt, len(tw))
	twInvS := make([]gnum.GaussianInt, len(twInv))
	for i := range tw {
		twS[i] = gnum.SForm(tw[i], mod)
		twInvS[i] = gnum.SForm(twInv[i], mod)
	}

	var cubeMod *num.Modulus
	if packLen > 1 {
		cubeMod = num.NewModulus(packLen)
	}

	rankInv := num.Inv(uint64(packLen<<1), mod)
	rankInvS := num.SForm(rankInv, mod)

	return &pow2Cyclotomic3Mod4Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		packIdx: packIdx,

		cube:    []int{packLen},
		cubeMod: []*num.Modulus{cubeMod},
		cubeGen: []uint64{5},

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv:  rankInv,
		rankInvS: rankInvS,

		pool: pool.NewPool(func() *[]gnum.GaussianInt {
			v := make([]gnum.GaussianInt, packLen<<1)
			return &v
		}),
	}
}

func (pck *pow2Cyclotomic3Mod4Packer) packTo(vPack, v []uint64) {
	if len(vPack) != pck.params.Rank() || len(v) > pck.packLen {
		panic("inconsistent inputs")
	}

	if len(v) == 0 {
		clear(vPack)
		return
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

	for i := 0; i < pck.packLen; i++ {
		idx0, idx1 := pck.packIdx[i], pck.packIdx[i+pck.packLen]

		vBuf[idx0].Real = v[i%len(v)]
		vBuf[idx0].Imag = 0

		vBuf[idx1].Real = v[i%len(v)]
		vBuf[idx1].Imag = 0
	}

	invNTTGaussIntInPlacePow2(vBuf, pck.twInv, pck.twInvS, pck.mod.Value())

	clear(vPack)
	skip := pck.params.Rank() / (pck.packLen << 1)
	for i := range vBuf {
		vPack[i*skip] = num.SMul(vBuf[i].Real, pck.rankInv, pck.rankInvS, pck.mod)
	}
}

func (pck *pow2Cyclotomic3Mod4Packer) unpackTo(v, vPack []uint64) {
	if len(v) == 0 || pck.packLen%len(v) != 0 || len(vPack) != pck.params.Rank() {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)
	clear(vBuf)

	skip := pck.params.Rank() / (pck.packLen << 1)
	for i := range vBuf {
		vBuf[i].Real = vPack[i*skip]
		vBuf[i].Imag = 0
	}

	fwdNTTGaussIntInPlacePow2(vBuf, pck.tw, pck.twS, pck.mod.Value())

	for i := range v {
		v[i] = num.Reduce4Q(vBuf[pck.packIdx[i]].Real, pck.mod)
	}
}

func (pck *pow2Cyclotomic3Mod4Packer) rotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("inconsistent inputs")
	}
	if pck.cubeMod[0] == nil {
		return 1
	}
	return int(num.Exp(5, num.Reduce(idx[0], pck.cubeMod[0]), pck.params.CycloIdxMod()))
}

type anyCyclotomicPacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// ntt is the transformer for NTT.
	ntt *dft.Transformer

	// cube is the form of the hypercube structure.
	cube []int
	// cubeMod is [num.Modulus] form of cube.
	// Entries are nil for dimensions of size one.
	cubeMod []*num.Modulus
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64
}

// newAnyCyclotomicPacker creates a new [anyCyclotomicPacker].
func newAnyCyclotomicPacker(params dft.RingParameters, mod *num.Modulus) *anyCyclotomicPacker {
	primes, exps := num.Factor(uint64(params.CycloIdx()))
	cubeGen := num.Generators(uint64(params.CycloIdx()))
	cube := make([]int, len(cubeGen))
	if primes[0] == 2 && exps[0] == 1 {
		for i := range cube {
			pExp := num.Exp(primes[i+1], exps[i+1], nil)
			cube[i] = int(pExp - pExp/primes[i+1])
			cubeGen[i] = num.Inv(cubeGen[i], params.CycloIdxMod())
		}
	} else if primes[0] == 2 && exps[0] > 2 {
		cube[0] = 1 << (exps[0] - 2)
		cubeGen[0] = num.Inv(cubeGen[0], params.CycloIdxMod())

		cube[1] = 2
		cubeGen[1] = num.Inv(cubeGen[1], params.CycloIdxMod())

		for i := 2; i < len(cube); i++ {
			pExp := num.Exp(primes[i-1], exps[i-1], nil)
			cube[i] = int(pExp - pExp/primes[i-1])
			cubeGen[i] = num.Inv(cubeGen[i], params.CycloIdxMod())
		}
	} else {
		for i := range cube {
			pExp := num.Exp(primes[i], exps[i], nil)
			cube[i] = int(pExp - pExp/primes[i])
			cubeGen[i] = num.Inv(cubeGen[i], params.CycloIdxMod())
		}
	}

	cubeMod := make([]*num.Modulus, len(cube))
	for i := range cube {
		if cube[i] > 1 {
			cubeMod[i] = num.NewModulus(cube[i])
		}
	}

	return &anyCyclotomicPacker{
		params: params,
		mod:    mod,

		packLen: int(params.Rank()),
		ntt:     dft.NewTransformer(params, mod),

		cube:    cube,
		cubeMod: cubeMod,
		cubeGen: cubeGen,
	}
}

func (pck *anyCyclotomicPacker) packTo(vPack, v []uint64) {
	if len(vPack) != pck.params.Rank() || len(v) != pck.packLen {
		panic("inconsistent inputs")
	}

	vec.MFormTo(vPack, v, pck.mod)
	pck.ntt.InverseTo(vPack, vPack)
}

func (pck *anyCyclotomicPacker) unpackTo(v, vPack []uint64) {
	if len(v) != pck.packLen || len(vPack) != pck.params.Rank() {
		panic("inconsistent inputs")
	}

	pck.ntt.ForwardTo(v, vPack)
	vec.InvMFormTo(v, v, pck.mod)
}

func (pck *anyCyclotomicPacker) rotIdxToAutIdx(idx []int) int {
	if len(idx) != len(pck.cubeGen) {
		panic("inconsistent inputs")
	}

	autIdx := uint64(1)
	for i := range pck.cubeMod {
		if pck.cubeMod[i] == nil {
			continue
		}
		exp := num.Reduce(idx[i], pck.cubeMod[i])
		autIdx = num.Mul(autIdx, num.Exp(pck.cubeGen[i], exp, pck.params.CycloIdxMod()), pck.params.CycloIdxMod())
	}
	return int(autIdx)
}
