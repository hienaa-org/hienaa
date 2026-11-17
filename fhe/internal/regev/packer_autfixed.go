package regev

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type pow2AutFixed1Mod4Packer struct {
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

// newPow2AutFixed1Mod4Packer creates a new [pow2AutFixed1Mod4Packer].
func newPow2AutFixed1Mod4Packer(params dft.RingParameters, mod *num.Modulus) *pow2AutFixed1Mod4Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloIdx() >> 2
	for i := range primes {
		bits := bits.TrailingZeros64(primes[i]-1) - 2
		if packLen > (1 << bits) {
			packLen = 1 << bits
		}
	}

	packLen64 := uint64(packLen)
	packLenQuadMod := num.NewModulus(packLen << 2)
	packIdx := make([]int, packLen)
	revShiftBits := 64 - uint64(num.Log2(packLen)+1)
	for i := uint64(0); i < packLen64; i++ {
		invPow5 := num.Exp(5, packLen64-i, packLenQuadMod)
		idx := invPow5 >> 1
		if invPow5 > packLen64<<1 {
			idx = packLen64<<1 - 1 - idx
		}

		idxOut := bits.Reverse64(idx) >> revShiftBits
		if idxOut >= packLen64 {
			idxOut = packLen64<<1 - 1 - idxOut
		}
		packIdx[i] = int(idxOut)
	}

	var cubeMod *num.Modulus
	if packLen > 1 {
		cubeMod = num.NewModulus(packLen)
	}

	return &pow2AutFixed1Mod4Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		packIdx: packIdx,
		ntt:     dft.NewTransformer(dft.NewAutFixedParameters(packLen<<2, packLen), mod),

		cube:    []int{packLen},
		cubeMod: []*num.Modulus{cubeMod},
		cubeGen: []uint64{5},

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, packLen)
			return &v
		}),
	}
}

func (pck *pow2AutFixed1Mod4Packer) packTo(vPack, v []uint64) {
	if len(vPack) != pck.params.Rank() || len(v) > pck.packLen {
		panic("inconsistent inputs")
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

func (pck *pow2AutFixed1Mod4Packer) unpackTo(v, vPack []uint64) {
	if len(v)%pck.packLen != 0 || len(vPack) != pck.params.Rank() {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

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

func (pck *pow2AutFixed1Mod4Packer) rotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("inconsistent inputs")
	}
	if pck.cubeMod[0] == nil {
		return 1
	}
	return int(num.Exp(5, num.Reduce(idx[0], pck.cubeMod[0]), pck.params.CycloIdxMod()))
}

type pow2AutFixed3Mod4Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// nttRank is the rank of the NTT.
	nttRank int
	// packIdx is the index mapping for the packing.
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

	pool *pool.Pool[*[]gnum.GaussianInt]
}

// newPow2AutFixed3Mod4Packer creates a new [pow2AutFixed3Mod4Packer].
func newPow2AutFixed3Mod4Packer(params dft.RingParameters, mod *num.Modulus) *pow2AutFixed3Mod4Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloIdx() >> 2
	logPackLen := bits.TrailingZeros64(primes[0]+1) - 1
	if packLen > (1 << logPackLen) {
		packLen = 1 << logPackLen
	}

	nttRank := packLen
	if int(primes[0])%params.CycloIdx() != params.CycloIdx()-1 {
		packLen >>= 1
	}

	nttRank64QuadMod, packLen64 := num.NewModulus(uint64(nttRank)<<2), uint64(packLen)
	packIdx := make([]int, nttRank)
	revShiftBits := 64 - int(num.Log2(nttRank)+1)
	for i := uint64(0); i < packLen64; i++ {
		idx0 := int(num.Exp(5, packLen64<<1-i, nttRank64QuadMod))
		packIdx[i] = min(idx0, nttRank<<2-idx0)
		if nttRank != packLen {
			idx1 := int(num.Mul(uint64(idx0), primes[0], nttRank64QuadMod))
			idx1 = min(idx1, nttRank<<2-idx1)
			idx0 = packIdx[i]
			packIdx[i], packIdx[i+packLen64] = min(idx0, idx1), max(idx0, idx1)
		}
	}

	for i := 0; i < nttRank; i++ {
		idx := int(bits.Reverse64(uint64(packIdx[i]>>1)) >> revShiftBits)
		packIdx[i] = min(idx, nttRank<<1-1-idx)
	}

	twLarge := make([]gnum.GaussianInt, nttRank<<1)
	twInvLarge := make([]gnum.GaussianInt, nttRank<<1)
	twLarge[0], twLarge[1] = gnum.GaussianInt{Real: 1, Imag: 0}, gaussIntNthRoot(nttRank<<2, mod)
	twInvLarge[0], twInvLarge[1] = gnum.GaussianInt{Real: 1, Imag: 0}, gnum.Inv(twLarge[1], mod)
	for i := 2; i < 2*nttRank; i++ {
		twLarge[i] = gnum.Mul(twLarge[i-1], twLarge[1], mod)
		twInvLarge[i] = gnum.Mul(twInvLarge[i-1], twInvLarge[1], mod)
	}
	vec.RadixReverseInPlace(twLarge, 2)
	vec.RadixReverseInPlace(twInvLarge, 2)

	tw := make([]gnum.GaussianInt, nttRank)
	twInv := make([]gnum.GaussianInt, nttRank)

	tw[0] = twLarge[1]
	twInv[0] = twInvLarge[1]
	for m := 1; m <= nttRank/2; m <<= 1 {
		copy(tw[m:2*m], twLarge[2*m:3*m])
		copy(twInv[m:2*m], twInvLarge[2*m:3*m])
	}

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

	return &pow2AutFixed3Mod4Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		nttRank: nttRank,
		packIdx: packIdx,

		cube:    []int{packLen},
		cubeMod: []*num.Modulus{cubeMod},
		cubeGen: []uint64{5},

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: num.Inv(uint64(nttRank)<<1, mod),

		pool: pool.NewPool(func() *[]gnum.GaussianInt {
			v := make([]gnum.GaussianInt, nttRank)
			return &v
		}),
	}
}

func (pck *pow2AutFixed3Mod4Packer) packTo(vPack, v []uint64) {
	if len(vPack) != pck.params.Rank() || len(v) > pck.packLen {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

	for i := 0; i < pck.packLen; i++ {
		idx := pck.packIdx[i]
		vBuf[idx].Real = v[i%len(v)]
		vBuf[idx].Imag = 0

		if pck.packLen != pck.nttRank {
			idx = pck.packIdx[i+pck.packLen]
			vBuf[idx].Real = v[i%len(v)]
			vBuf[idx].Imag = 0
		}
	}

	invNTTGaussIntInPlacePow2(vBuf, pck.twInv, pck.twInvS, pck.mod.Value())

	clear(vPack)
	skip := pck.params.Rank() / pck.nttRank
	vPack[0] = vBuf[0].Real << 1
	twoQ := pck.mod.Value() << 1
	for i := 1; i < pck.nttRank; i++ {
		uReal := num.SMulLazy(vBuf[pck.nttRank-i].Real, pck.tw[0].Real, pck.twS[0].Real, pck.mod) + twoQ
		uReal -= num.SMulLazy(vBuf[pck.nttRank-i].Imag, pck.tw[0].Imag, pck.twS[0].Imag, pck.mod)
		if uReal >= twoQ {
			uReal -= twoQ
		}
		vPack[i*skip] = vBuf[i].Real + uReal
	}

	vec.MulTo(vPack, vPack, pck.rankInv, pck.mod)
}

func (pck *pow2AutFixed3Mod4Packer) unpackTo(v, vPack []uint64) {
	if len(v)%pck.packLen != 0 || len(vPack) != pck.params.Rank() {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

	skip := pck.params.Rank() / pck.nttRank
	vBuf[0].Real = vPack[0]
	vBuf[0].Imag = 0
	for i := 1; i < pck.nttRank; i++ {
		uReal := num.SMul(vPack[(pck.nttRank-i)*skip], pck.tw[0].Real, pck.twS[0].Real, pck.mod)
		uImag := num.SMul(vPack[(pck.nttRank-i)*skip], pck.tw[0].Imag, pck.twS[0].Imag, pck.mod)
		vBuf[i].Real = num.Sub(vPack[i*skip], uReal, pck.mod)
		vBuf[i].Imag = num.Neg(uImag, pck.mod)
	}

	fwdNTTGaussIntInPlacePow2(vBuf, pck.tw, pck.twS, pck.mod.Value())

	for i := 0; i < len(v); i++ {
		v[i] = num.Reduce4Q(vBuf[pck.packIdx[i]].Real, pck.mod)
	}
}

func (pck *pow2AutFixed3Mod4Packer) rotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("inconsistent inputs")
	}
	if pck.cubeMod[0] == nil {
		return 1
	}
	return int(num.Exp(5, num.Reduce(idx[0], pck.cubeMod[0]), pck.params.CycloIdxMod()))
}

type primeAutFixedPacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int

	// cube is the form of the hypercube structure.
	cube []int
	// cubeMod is [num.Modulus] form of cube.
	// Entries are nil for dimensions of size one.
	cubeMod []*num.Modulus
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	// resol is the resolution of unity.
	resol [][]uint64
	// resolInv is the inverse of resolution of unity.
	resolInv [][]uint64

	// ambRank is the rank of the ambient NTT.
	ambRank int
	// ambMod is the ambient modulus.
	ambMod []*num.Modulus
	// ambNTT is the ambient NTT.
	ambNTT []*dft.Transformer
	// embedder is the embedder for the packing.
	embedder *crt.VecEmbedder

	pool *pool.Pool[*[][]uint64]
}

// newPrimeAutFixedPacker creates a new [primeAutFixedPacker].
func newPrimeAutFixedPacker(params dft.RingParameters, mod *num.Modulus) *primeAutFixedPacker {
	// Compute the resolution of unity.
	primes, exps := num.Factor(mod.Value())
	if len(primes) != 1 {
		panic("mod must be a prime power")
	}
	prime, exp := primes[0], exps[0]
	resol := resolutionOfUnity(params, prime, exp)
	resolInv := make([]uint64, len(resol))
	ord := num.Order(prime, uint64(params.CycloIdx()))
	if ord&1 == 1 {
		copy(resolInv[len(resol)/2:], resol[:len(resol)/2])
		copy(resolInv[:len(resol)/2], resol[len(resol)/2:])
	} else {
		copy(resolInv, resol)
	}

	vec.MulTo(resolInv, resolInv, num.Reduce(params.CycloIdx(), mod), mod)
	vec.AddTo(resolInv, resolInv, num.Reduce(ord, mod), mod)

	// Reduce the resolution of unity to the packing length.
	packLen := num.GCD(params.Rank(), len(resol))
	for i := 1; i < len(resol)/packLen; i++ {
		vec.AddTo(resol[:packLen], resol[:packLen], resol[packLen*i:packLen*(i+1)], mod)
		vec.AddTo(resolInv[:packLen], resolInv[:packLen], resolInv[packLen*i:packLen*(i+1)], mod)
	}
	resol = resol[:packLen]
	resolInv = resolInv[:packLen]

	var ambRank int
	if num.IsProdPowerOf(packLen, []int{2, 3, 5}) {
		ambRank = packLen
	} else {
		ambRank = num.NextProdPower((2*packLen-1)+1, []int{2})
	}

	ambParams := dft.NewCyclicParameters(ambRank)
	maxBits := 2*num.Log2(mod.Value()) + num.Log2(ambRank)
	ambMod := dft.AmbientNTTPrimes(ambParams, maxBits)
	ambNTT := make([]*dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	resolMod := make([][]uint64, len(ambMod))
	resolInvMod := make([][]uint64, len(ambMod))
	for i := 0; i < len(ambMod); i++ {
		resolMod[i] = make([]uint64, ambRank)
		resolInvMod[i] = make([]uint64, ambRank)

		copy(resolMod[i][:packLen], resol)
		copy(resolInvMod[i][:packLen], resolInv)

		ambNTT[i].ForwardTo(resolMod[i], resolMod[i])
		ambNTT[i].ForwardTo(resolInvMod[i], resolInvMod[i])
	}

	cubeGen := uint64(1)
	if params.Rank() > 1 {
		cubeGen = num.Inv(num.Generators(uint64(params.CycloIdx()))[0], params.CycloIdxMod())
	}
	var cubeMod *num.Modulus
	if packLen > 1 {
		cubeMod = num.NewModulus(packLen)
	}

	return &primeAutFixedPacker{
		params: params,
		mod:    mod,

		packLen: packLen,

		cube:    []int{packLen},
		cubeMod: []*num.Modulus{cubeMod},
		cubeGen: []uint64{cubeGen},

		resol:    resolMod,
		resolInv: resolInvMod,

		ambRank:  ambRank,
		ambMod:   ambMod,
		ambNTT:   ambNTT,
		embedder: crt.NewVecEmbedder([]*num.Modulus{mod}, ambMod),

		pool: pool.NewPool(func() *[][]uint64 {
			v := make([][]uint64, len(ambMod))
			for i := range v {
				v[i] = make([]uint64, ambRank)
			}
			return &v
		}),
	}
}

func (pck *primeAutFixedPacker) packTo(vPack, v []uint64) {
	if len(vPack) != pck.params.Rank() || pck.packLen%len(v) != 0 {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

	for i := range pck.ambMod {
		for j := 0; j < pck.packLen/len(v); j++ {
			vBuf[i][j*len(v)] = v[0]
			copy(vBuf[i][j*len(v)+1:(j+1)*len(v)], v[1:])
			slices.Reverse(vBuf[i][j*len(v)+1 : (j+1)*len(v)])
		}
		clear(vBuf[i][pck.packLen:])

		pck.ambNTT[i].ForwardTo(vBuf[i], vBuf[i])
		vec.MMulLazyTo(vBuf[i], vBuf[i], pck.resol[i], pck.ambMod[i])
		pck.ambNTT[i].InverseTo(vBuf[i], vBuf[i])
	}

	pck.embedder.EmbedTo(vBuf[:1], vBuf)
	if pck.packLen != pck.ambRank {
		vec.AddTo(vBuf[0][:pck.packLen-1], vBuf[0][:pck.packLen-1], vBuf[0][pck.packLen:2*pck.packLen-1], pck.mod)
	}

	for i := 0; i < pck.params.Rank()/pck.packLen; i++ {
		copy(vPack[i*pck.packLen:(i+1)*pck.packLen], vBuf[0][:pck.packLen])
	}
}

func (pck *primeAutFixedPacker) unpackTo(v, vPack []uint64) {
	if len(v) != pck.packLen || len(vPack) != pck.params.Rank() {
		panic("inconsistent inputs")
	}

	vBufPtr := pck.pool.Get()
	vBuf := *vBufPtr
	defer pck.pool.Put(vBufPtr)

	for i := range pck.ambMod {
		vBuf[i][0] = vPack[0]
		copy(vBuf[i][1:pck.packLen], vPack[1:pck.packLen])
		slices.Reverse(vBuf[i][1:pck.packLen])
		clear(vBuf[i][pck.packLen:])

		pck.ambNTT[i].ForwardTo(vBuf[i], vBuf[i])
		vec.MMulLazyTo(vBuf[i], vBuf[i], pck.resolInv[i], pck.ambMod[i])
		pck.ambNTT[i].InverseTo(vBuf[i], vBuf[i])
	}
	pck.embedder.EmbedTo(vBuf[:1], vBuf)

	if pck.packLen != pck.ambRank {
		vec.AddTo(vBuf[0][:pck.packLen-1], vBuf[0][:pck.packLen-1], vBuf[0][pck.packLen:2*pck.packLen-1], pck.mod)
	}

	copy(v, vBuf[0][:len(v)])
}

func (pck *primeAutFixedPacker) rotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("inconsistent inputs")
	}
	if pck.cubeMod[0] == nil {
		return 1
	}
	return int(num.Exp(pck.cubeGen[0], num.Reduce(idx[0], pck.cubeMod[0]), pck.params.CycloIdxMod()))
}
