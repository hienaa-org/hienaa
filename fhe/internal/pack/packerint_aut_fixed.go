package pack

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/fhe/internal/gr"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// pow2AutFixedMod1IntPacker is a packer for the power-of-two autfixed ring,
// where the modulus is a multiple of prime powers that are 1 mod 4.
type pow2AutFixedMod1IntPacker struct {
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

	cycloIdxMod *num.Modulus

	pool *pool.Pool[*[]uint64]
}

// newPow2AutFixedMod1IntPacker creates a new [pow2AutFixedMod1IntPacker].
func newPow2AutFixedMod1IntPacker(params dft.RingParameters, mod *num.Modulus) *pow2AutFixedMod1IntPacker {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloIndex() >> 2
	for i := range primes {
		ithLogPackLen := bits.TrailingZeros64(primes[i]-1) - 2
		if packLen > (1 << ithLogPackLen) {
			packLen = 1 << ithLogPackLen
		}
	}

	nttParams := dft.NewAutFixedParameters(packLen<<2, packLen)
	ntt := dft.NewTransformer(nttParams, mod)

	return &pow2AutFixedMod1IntPacker{
		params: params,
		mod:    mod,

		packLen: packLen,
		ntt:     ntt,

		cube:    []int{packLen},
		cubeGen: []uint64{5},

		cycloIdxMod: num.NewModulus(params.CycloIndex()),

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, packLen)
			return &v
		}),
	}
}

// Params returns the ring parameters.
func (p *pow2AutFixedMod1IntPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus used for the packing/unpacking.
func (p *pow2AutFixedMod1IntPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *pow2AutFixedMod1IntPacker) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *pow2AutFixedMod1IntPacker) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *pow2AutFixedMod1IntPacker) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) > p.packLen {
		panic("input(s) shape not consistent")
	}

	vLen := len(v)

	vBufPtr := p.pool.Get()
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	vBufPow5Ptr := p.pool.Get()
	vBufPow5 := *vBufPow5Ptr
	defer p.pool.Put(vBufPow5Ptr)

	for i := 0; i < p.packLen/vLen; i++ {
		copy(vBuf[vLen*i:vLen*(i+1)], v)
	}

	invPow5 := 1
	inv5 := int(num.Inv(5, p.cycloIdxMod))
	mask := p.packLen<<2 - 1
	revShiftBits := 64 - int(num.Log2(uint64(p.packLen))+1)
	for i := 0; i < p.packLen; i++ {
		idx := invPow5 >> 1
		if invPow5 > p.packLen<<1 {
			idx = p.packLen<<1 - 1 - idx
		}

		idxOut := int(bits.Reverse64(uint64(idx)) >> revShiftBits)
		if idxOut >= p.packLen {
			idxOut = p.packLen<<1 - 1 - idxOut
		}

		vBufPow5[idxOut] = vBuf[i]
		invPow5 = (inv5 * invPow5) & mask
	}

	vec.MFormTo(vBufPow5, vBufPow5, p.mod)
	p.ntt.InverseTo(vBufPow5, vBufPow5)

	clear(vPack)
	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		vPack[i*skip] = vBufPow5[i]
	}
}

// UnPack returns the unpacking of vPack.
func (p *pow2AutFixedMod1IntPacker) UnPack(vPack []uint64) []uint64 {
	v := make([]uint64, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

// UnPackTo unpacks vPack to v.
func (p *pow2AutFixedMod1IntPacker) UnPackTo(v, vPack []uint64) {
	if len(v)%p.packLen != 0 || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	vLen := len(v)

	vBufPtr := p.pool.Get()
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	vBufPow5Ptr := p.pool.Get()
	vBufPow5 := *vBufPow5Ptr
	defer p.pool.Put(vBufPow5Ptr)

	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		vBuf[i] = vPack[i*skip]
	}
	p.ntt.ForwardTo(vBuf, vBuf)
	vec.InvMFormTo(vBuf, vBuf, p.mod)

	invPow5 := 1
	inv5 := int(num.Inv(5, p.cycloIdxMod))
	mask := p.packLen<<2 - 1
	revShiftBits := 64 - int(num.Log2(uint64(p.packLen))+1)
	for i := 0; i < p.packLen; i++ {
		idx := invPow5 >> 1
		if invPow5 > p.packLen<<1 {
			idx = p.packLen<<1 - 1 - idx
		}

		idxIn := int(bits.Reverse64(uint64(idx)) >> revShiftBits)
		if idxIn >= p.packLen {
			idxIn = p.packLen<<1 - 1 - idxIn
		}

		vBufPow5[i] = vBuf[idxIn]
		invPow5 = (inv5 * invPow5) & mask
	}

	copy(v, vBufPow5[:vLen])
}

// Cube returns the form of the hypercube structure.
func (p *pow2AutFixedMod1IntPacker) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *pow2AutFixedMod1IntPacker) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *pow2AutFixedMod1IntPacker) RotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("input(s) shape not consistent")
	}
	return int(num.Exp(5, uint64(idx[0]), nil)) & (p.params.CycloIndex() - 1)
}

// pow2AutFixedMod3IntPacker is a packer for the power-of-two autfixed ring,
// where the modulus is a multiple of prime powers that are 3 mod 4.
type pow2AutFixedMod3IntPacker struct {
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
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	// tw is the twiddle factor for NTT.
	tw []gnum.GaussianInt
	// twInv is the twiddle factor for InvNTT.
	twInv []gnum.GaussianInt

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	pool *pool.Pool[*[]gnum.GaussianInt]
}

// newPow2AutFixedMod3IntPacker creates a new [pow2AutFixedMod3IntPacker].
func newPow2AutFixedMod3IntPacker(params dft.RingParameters, mod *num.Modulus) *pow2AutFixedMod3IntPacker {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloIndex() >> 2
	logPackLen := bits.TrailingZeros64(primes[0]+1) - 1
	if packLen > (1 << logPackLen) {
		packLen = 1 << logPackLen
	}

	nttRank := packLen
	if int(primes[0])%params.CycloIndex() != params.CycloIndex()-1 {
		packLen >>= 1
	}

	twLarge := make([]gnum.GaussianInt, nttRank<<1)
	twInvLarge := make([]gnum.GaussianInt, nttRank<<1)
	twLarge[0], twLarge[1] = gnum.GaussianInt{Real: 1, Imag: 0}, gIntNthRoot(nttRank<<2, mod)
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

	mask := nttRank<<2 - 1
	packIdx := make([]int, nttRank)
	revShiftBits := 64 - int(num.Log2(uint64(nttRank))+1)
	if nttRank == packLen {
		for i := 0; i < packLen; i++ {
			packIdx[i] = int(num.Exp(5, uint64(packLen<<1-i), nil)) & mask
			if packIdx[i] > packLen<<1 {
				packIdx[i] = mask - packIdx[i] + 1
			}

			packIdx[i] >>= 1
			packIdx[i] = int(bits.Reverse64(uint64(packIdx[i])) >> revShiftBits)
			if packIdx[i] >= nttRank {
				packIdx[i] = nttRank<<1 - 1 - packIdx[i]
			}
		}
	} else {
		for i := 0; i < packLen; i++ {
			idx1 := int(num.Exp(5, uint64(packLen<<1-i), nil)) & mask
			idx2 := (idx1 * int(primes[0])) & mask
			idx3 := mask - idx1 + 1
			idx4 := mask - idx2 + 1

			if idx1 < packLen<<1 {
				packIdx[i] = idx1
			} else if idx2 < packLen<<1 {
				packIdx[i] = idx2
			} else if idx3 < packLen<<1 {
				packIdx[i] = idx3
			} else {
				packIdx[i] = idx4
			}

			if packIdx[i] < idx1 && idx1 < packLen<<2 {
				packIdx[i+packLen] = idx1
			}
			if packIdx[i] < idx2 && idx2 < packLen<<2 {
				packIdx[i+packLen] = idx2
			}
			if packIdx[i] < idx3 && idx3 < packLen<<2 {
				packIdx[i+packLen] = idx3
			}
			if packIdx[i] < idx4 && idx4 < packLen<<2 {
				packIdx[i+packLen] = idx4
			}

			packIdx[i] >>= 1
			packIdx[i+packLen] >>= 1

			packIdx[i] = int(bits.Reverse64(uint64(packIdx[i])) >> revShiftBits)
			if packIdx[i] >= nttRank {
				packIdx[i] = nttRank<<1 - 1 - packIdx[i]
			}

			packIdx[i+packLen] = int(bits.Reverse64(uint64(packIdx[i+packLen])) >> revShiftBits)
			if packIdx[i+packLen] >= nttRank {
				packIdx[i+packLen] = nttRank<<1 - 1 - packIdx[i+packLen]
			}
		}
	}

	return &pow2AutFixedMod3IntPacker{
		params: params,
		mod:    mod,

		packLen: packLen,
		nttRank: nttRank,
		packIdx: packIdx,

		cube:    []int{packLen},
		cubeGen: []uint64{5},

		tw:    tw,
		twInv: twInv,

		rankInv: num.Inv(uint64(nttRank<<1), mod),

		pool: pool.NewPool(func() *[]gnum.GaussianInt {
			v := make([]gnum.GaussianInt, nttRank)
			return &v
		}),
	}
}

// Params returns the ring parameters.
func (p *pow2AutFixedMod3IntPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus used for the packing/unpacking.
func (p *pow2AutFixedMod3IntPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *pow2AutFixedMod3IntPacker) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *pow2AutFixedMod3IntPacker) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *pow2AutFixedMod3IntPacker) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) > p.packLen {
		panic("input(s) shape not consistent")
	}

	vLen := len(v)

	vBufPtr := p.pool.Get()
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	for i := 0; i < p.packLen; i++ {
		idx := p.packIdx[i]
		vBuf[idx].Real = v[i&(vLen-1)]
		vBuf[idx].Imag = 0

		if p.packLen != p.nttRank {
			idx = p.packIdx[i+p.packLen]
			vBuf[idx].Real = v[i&(vLen-1)]
			vBuf[idx].Imag = 0
		}
	}

	invNTTGaloisRingInPlacePow2(vBuf, p.twInv, p.mod)

	var u gnum.GaussianInt
	skip := p.params.Rank() / p.nttRank
	vPack[0] = num.Add(vBuf[0].Real, vBuf[0].Real, p.mod)
	for i := 1; i < p.nttRank; i++ {
		u = gnum.Mul(vBuf[p.nttRank-i], p.tw[0], p.mod)
		vPack[i*skip] = num.Add(vBuf[i].Real, u.Real, p.mod)
	}

	vec.MulScalarTo(vPack, vPack, p.rankInv, p.mod)
}

// UnPack returns the unpacking of vPack.
func (p *pow2AutFixedMod3IntPacker) UnPack(vPack []uint64) []uint64 {
	v := make([]uint64, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

// UnPackTo unpacks vPack to v.
func (p *pow2AutFixedMod3IntPacker) UnPackTo(v, vPack []uint64) {
	if len(v)%p.packLen != 0 || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	vLen := len(v)

	vBufPtr := p.pool.Get()
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	skip := p.params.Rank() / p.nttRank
	vBuf[0].Real = vPack[0]
	vBuf[0].Imag = 0
	for i := 1; i < p.nttRank; i++ {
		real := num.Mul(p.tw[0].Real, vPack[(p.nttRank-i)*skip], p.mod)
		imag := num.Mul(p.tw[0].Imag, vPack[(p.nttRank-i)*skip], p.mod)
		vBuf[i].Real = num.Sub(vPack[i*skip], real, p.mod)
		vBuf[i].Imag = num.Neg(imag, p.mod)
	}

	nttGaloisRingInPlacePow2(vBuf, p.tw, p.mod)

	for i := 0; i < vLen; i++ {
		v[i] = vBuf[p.packIdx[i]].Real
	}
}

// Cube returns the form of the hypercube structure.
func (p *pow2AutFixedMod3IntPacker) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *pow2AutFixedMod3IntPacker) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *pow2AutFixedMod3IntPacker) RotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("input(s) shape not consistent")
	}
	return int(num.Exp(5, uint64(idx[0]), nil)) & (p.params.CycloIndex() - 1)
}

// primeAutFixedIntPacker is a packer for prime autfixed ring.
type primeAutFixedIntPacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int

	// cube is the form of the hypercube structure.
	cube []int
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	// resol is the resolution of unity.
	resol [][]uint64
	// invResol is the inverse resolution of unity.
	invResol [][]uint64

	// cycloIdxMod is the cyclotomic index modulus.
	cycloIdxMod *num.Modulus

	// ambRank is the rank of the ambient NTT.
	ambRank int
	// ambMod is the ambient modulus.
	ambMod []*num.Modulus
	// ambNTT is the ambient NTT.
	ambNTT []dft.Transformer
	// embedder is the embedder for the packing.
	embedder *crt.VecEmbedder

	pool *pool.Pool[*[][]uint64]
}

// newAutFixedPrimeIntPacker creates a new [primeAutFixedIntPacker].
func newAutFixedPrimeIntPacker(params dft.RingParameters, mod *num.Modulus) *primeAutFixedIntPacker {
	cycloIdx := params.CycloIndex()

	// Compute the resolution of unity.
	primes, exps := num.Factor(mod.Value())
	if len(primes) != 1 {
		panic("mod must be a prime power")
	}
	prime, exp := primes[0], exps[0]
	resolution := findResolutionOfUnity(cycloIdx, prime, exp)
	invResolution := make([]uint64, len(resolution))
	ord := num.Order(prime, num.NewModulus(cycloIdx))
	if ord&1 == 1 {
		copy(invResolution[len(resolution)/2:], resolution[:len(resolution)/2])
		copy(invResolution[:len(resolution)/2], resolution[len(resolution)/2:])
	} else {
		copy(invResolution, resolution)
	}

	vec.MulScalarTo(invResolution, invResolution, num.Reduce(uint64(cycloIdx), mod), mod)
	vec.AddScalarTo(invResolution, invResolution, num.Reduce(ord, mod), mod)

	// Reduce the resolution of unity to the packing length.
	packLen := num.GCD(params.Rank(), len(resolution))
	for i := 1; i < len(resolution)/packLen; i++ {
		vec.AddTo(resolution[:packLen], resolution[:packLen], resolution[packLen*i:packLen*(i+1)], mod)
		vec.AddTo(invResolution[:packLen], invResolution[:packLen], invResolution[packLen*i:packLen*(i+1)], mod)
	}
	resolution = resolution[:packLen]
	invResolution = invResolution[:packLen]

	var ambRank int
	if num.IsProdPowerOf(packLen, []int{2, 3, 5}) {
		ambRank = packLen
	} else {
		ambRank = num.NextProdPower(2*packLen-1, []int{2})
	}

	ambParams := dft.NewCyclicParameters(ambRank)
	maxBits := 2*num.Log2(mod.Value()) + num.Log2(ambRank)
	ambMod := dft.MustFindAmbientPrimes(ambParams, maxBits)
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	resol := make([][]uint64, len(ambMod))
	invResol := make([][]uint64, len(ambMod))
	for i := 0; i < len(ambMod); i++ {
		resol[i] = make([]uint64, ambRank)
		invResol[i] = make([]uint64, ambRank)

		copy(resol[i][:packLen], resolution)
		copy(invResol[i][:packLen], invResolution)

		ambNTT[i].ForwardTo(resol[i], resol[i])
		ambNTT[i].ForwardTo(invResol[i], invResol[i])
	}

	embedder := crt.NewVecEmbedder([]*num.Modulus{mod}, ambMod)

	cycloIdxMod := num.NewModulus(cycloIdx)

	return &primeAutFixedIntPacker{
		params: params,
		mod:    mod,

		packLen: packLen,

		cube:    []int{packLen},
		cubeGen: []uint64{num.Inv(num.Generators(cycloIdxMod)[0], cycloIdxMod)},

		resol:    resol,
		invResol: invResol,

		cycloIdxMod: cycloIdxMod,

		ambRank:  ambRank,
		ambMod:   ambMod,
		ambNTT:   ambNTT,
		embedder: embedder,

		pool: pool.NewPool(func() *[][]uint64 {
			v := make([][]uint64, len(ambMod))
			for i := range v {
				v[i] = make([]uint64, ambRank)
			}
			return &v
		}),
	}
}

// findResolutionOfUnity finds the resolution of unity for the given cyclotomic index, prime, and exponent.
// The algorithm is from https://eprint.iacr.org/2024/2032.
func findResolutionOfUnity(cycloIdx int, prime uint64, exp uint64) []uint64 {
	if !num.IsPrime(cycloIdx) {
		panic("cycloIdx must be a prime number")
	}
	if !num.IsPrime(prime) {
		panic("prime must be a prime number")
	}

	cycloIdxMod := num.NewModulus(cycloIdx)
	ord := num.Order(prime, cycloIdxMod)
	rank := int(num.Totient(uint64(cycloIdx)) / ord)
	r := gr.NewGaloisRing(uint64(num.Exp(prime, exp, nil)), int(ord))
	root := grNthRoot(r, cycloIdx)

	// Compute the first factor of the cyclotomic polynomial.
	factorPoly := make([]*gr.Element, int(ord)+1)
	tmpPoly := make([]*gr.Element, int(ord)+1)
	for i := 0; i <= int(ord); i++ {
		if i == 0 {
			factorPoly[i] = r.NewElementFromUint64(1)
		} else {
			factorPoly[i] = r.NewElement()
		}
		tmpPoly[i] = r.NewElement()
	}

	rootPow := r.NewElement()
	for i := 0; i < int(ord); i++ {
		r.ExpTo(rootPow, root, num.Exp(prime, uint64(i), cycloIdxMod))

		// tmpPoly = factor * rootPow
		for j := 0; j <= i; j++ {
			r.MulTo(tmpPoly[j], factorPoly[j], rootPow)
		}

		// factor = factor * X
		for j := i + 1; j > 0; j-- {
			factorPoly[j].CopyFrom(factorPoly[j-1])
		}
		factorPoly[0].Clear()

		// factor = factor - tmpPoly
		for j := 0; j <= i; j++ {
			r.SubTo(factorPoly[j], factorPoly[j], tmpPoly[j])
		}
	}

	factor := make([]uint64, int(ord)+1)
	factorInt := make([]int64, int(ord)+1)
	for i := 0; i <= int(ord); i++ {
		factor[i] = uint64(factorPoly[i].Coeffs()[0])
		factorInt[i] = int64(factorPoly[i].Coeffs()[0])
	}

	cycloPoly := dft.CyclotomicPolynomial(cycloIdx)
	cycloPolyMod := make([]uint64, len(cycloPoly))
	modulus := num.NewModulus(num.Exp(prime, exp, nil))
	for i := 0; i < len(cycloPoly); i++ {
		cycloPolyMod[i] = num.Reduce(cycloPoly[i], modulus)
	}

	// remPoly = (cycloPoly/factor) % factor.
	quoPoly, _ := quoRem(cycloPolyMod, factor, modulus)
	_, remPoly := quoRem(quoPoly, factor, modulus)

	// remEl = (1/remPoly) % factor.
	rFactor := gr.NewGaloisRingCustom(modulus.Value(), factorInt)
	remEl := rFactor.NewElement()
	for i := 0; i < len(remPoly); i++ {
		remEl.Coeffs()[i] = remPoly[i]
	}
	rFactor.InvTo(remEl, remEl)

	// Compute the resolution of unity.
	tf := crt.NewOperator(dft.NewCyclicParameters(num.NextProdPower(cycloIdx, []int{2})), []*num.Modulus{modulus})
	resolPoly := tf.NewPoly()
	remInvPoly := tf.NewPoly()

	// resolPoly = quoPoly * remEl.
	copy(resolPoly.Coeffs[0][:len(quoPoly)], quoPoly)
	copy(remInvPoly.Coeffs[0][:len(remEl.Coeffs())], remEl.Coeffs())

	tf.FwdNTTTo(resolPoly, resolPoly)
	tf.FwdNTTTo(remInvPoly, remInvPoly)
	tf.MulTo(resolPoly, resolPoly, remInvPoly)
	tf.InvNTTTo(resolPoly, resolPoly)

	// Extract the meaningful values.
	resolution := make([]uint64, rank)
	gen := num.Generators(cycloIdxMod)[0]
	idx := uint64(1)
	for i := 0; i < rank; i++ {
		resolution[i] = resolPoly.Coeffs[0][idx] + modulus.Value() - resolPoly.Coeffs[0][0]
		if resolution[i] >= modulus.Value() {
			resolution[i] -= modulus.Value()
		}
		idx = num.Mul(idx, gen, cycloIdxMod)
	}

	return resolution
}

// Params returns the ring parameters.
func (p *primeAutFixedIntPacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus used for the packing/unpacking.
func (p *primeAutFixedIntPacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the length of the packing/unpacking.
func (p *primeAutFixedIntPacker) PackLen() int {
	return p.packLen
}

// Pack returns the packing of v.
func (p *primeAutFixedIntPacker) Pack(v []uint64) []uint64 {
	vPack := make([]uint64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// PackTo packs v to vPack.
func (p *primeAutFixedIntPacker) PackTo(vPack, v []uint64) {
	if len(vPack) != p.params.Rank() || len(v) > p.packLen {
		panic("input(s) shape not consistent")
	}

	vBufPtr := p.pool.Get()
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	for i := range p.ambMod {
		for j := 0; j < p.packLen/len(v); j++ {
			vBuf[i][j*len(v)] = v[0]
			copy(vBuf[i][j*len(v)+1:(j+1)*len(v)], v[1:])
			slices.Reverse(vBuf[i][j*len(v)+1 : (j+1)*len(v)])
		}
		clear(vBuf[i][p.packLen:])

		p.ambNTT[i].ForwardTo(vBuf[i], vBuf[i])
		vec.MMulLazyTo(vBuf[i], vBuf[i], p.resol[i], p.ambMod[i])
		p.ambNTT[i].InverseTo(vBuf[i], vBuf[i])
	}

	p.embedder.EmbedTo(vBuf[:1], vBuf)
	if p.packLen != p.ambRank {
		vec.AddTo(vBuf[0][:p.packLen-1], vBuf[0][:p.packLen-1], vBuf[0][p.packLen:2*p.packLen-1], p.mod)
	}

	for i := 0; i < p.params.Rank()/p.packLen; i++ {
		copy(vPack[i*p.packLen:(i+1)*p.packLen], vBuf[0][:p.packLen])
	}
}

// UnPack returns the unpacking of vPack.
func (p *primeAutFixedIntPacker) UnPack(vPack []uint64) []uint64 {
	v := make([]uint64, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

// UnPackTo unpacks vPack to v.
func (p *primeAutFixedIntPacker) UnPackTo(v, vPack []uint64) {
	if len(v) != p.packLen || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	vBufPtr := p.pool.Get()
	vBuf := *vBufPtr
	defer p.pool.Put(vBufPtr)

	for i := range p.ambMod {
		vBuf[i][0] = vPack[0]
		copy(vBuf[i][1:p.packLen], vPack[1:p.packLen])
		slices.Reverse(vBuf[i][1:p.packLen])
		clear(vBuf[i][p.packLen:])

		p.ambNTT[i].ForwardTo(vBuf[i], vBuf[i])
		vec.MMulLazyTo(vBuf[i], vBuf[i], p.invResol[i], p.ambMod[i])
		p.ambNTT[i].InverseTo(vBuf[i], vBuf[i])
	}
	p.embedder.EmbedTo(vBuf[:1], vBuf)

	if p.packLen != p.ambRank {
		vec.AddTo(vBuf[0][:p.packLen-1], vBuf[0][:p.packLen-1], vBuf[0][p.packLen:2*p.packLen-1], p.mod)
	}

	copy(v, vBuf[0][:len(v)])
}

// Cube returns the form of the hypercube structure.
func (p *primeAutFixedIntPacker) Cube() []int {
	return p.cube
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *primeAutFixedIntPacker) CubeGen() []uint64 {
	return p.cubeGen
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *primeAutFixedIntPacker) RotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("input(s) shape not consistent")
	}
	return int(num.Exp(p.cubeGen[0], uint64(idx[0]), p.cycloIdxMod))
}
