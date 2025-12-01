package pack

import (
	"math"
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/fhe/internal/gr"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type autFixedPow2Mod1Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// ntt is the transformer for NTT.
	ntt dft.Transformer

	buf packerBuffer
}

// newAutFixedPow2Mod1Packer creates a new [autFixedPow2Mod1Packer].
func newAutFixedPow2Mod1Packer(params dft.RingParameters, mod *num.Modulus) *autFixedPow2Mod1Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloOrder() >> 2
	for i := range primes {
		ithLogPackLen := bits.TrailingZeros64(primes[i]-1) - 2
		if packLen > (1 << ithLogPackLen) {
			packLen = 1 << ithLogPackLen
		}
	}

	nttParams := dft.NewAutFixedParameters(packLen<<2, packLen)
	ntt := dft.NewTransformer(nttParams, mod)

	return &autFixedPow2Mod1Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		ntt:     ntt,

		buf: newPackerBuffer(2, packLen),
	}
}

// Params returns the parameters of the packer.
func (p *autFixedPow2Mod1Packer) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *autFixedPow2Mod1Packer) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *autFixedPow2Mod1Packer) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *autFixedPow2Mod1Packer) SafeCopy() PackerInt {
	return &autFixedPow2Mod1Packer{
		params:  p.params,
		mod:     p.mod,
		packLen: p.packLen,

		ntt: p.ntt.SafeCopy(),

		buf: newPackerBuffer(2, p.packLen),
	}
}

// Pack packs the input vector into a polynomial.
func (p *autFixedPow2Mod1Packer) Pack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.params.Rank())
	p.PackTo(vOut, vIn)
	return vOut
}

// PackTo packs the input vector into a polynomial.
func (p *autFixedPow2Mod1Packer) PackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vIn)

	for i := 0; i < p.packLen/vLen; i++ {
		copy(p.buf.coeffs[0][vLen*i:vLen*(i+1)], vIn)
	}

	pow5 := 1
	mask := p.packLen<<2 - 1
	revShiftBits := 64 - int(num.Log2(uint64(p.packLen))+1)
	for i := 0; i < p.packLen; i++ {
		idx := pow5 >> 1
		if pow5 > p.packLen<<1 {
			idx = p.packLen<<1 - 1 - idx
		}

		idxOut := int(bits.Reverse64(uint64(idx)) >> revShiftBits)
		if idxOut >= p.packLen {
			idxOut = p.packLen<<1 - 1 - idxOut
		}

		p.buf.coeffs[1][idxOut] = p.buf.coeffs[0][i]
		pow5 = (5 * pow5) & mask
	}

	vec.MFormTo(p.buf.coeffs[0], p.buf.coeffs[1], p.mod)
	p.ntt.InverseTo(p.buf.coeffs[0], p.buf.coeffs[0])

	clear(vOut)
	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		vOut[i*skip] = p.buf.coeffs[0][i]
	}
}

// UnPack unpacks the polynomial into a vector.
func (p *autFixedPow2Mod1Packer) UnPack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, vIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *autFixedPow2Mod1Packer) UnPackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vOut)

	skip := p.params.Rank() / p.packLen
	for i := 0; i < p.params.Rank()/skip; i++ {
		p.buf.coeffs[0][i] = vIn[i*skip]
	}
	p.ntt.ForwardTo(p.buf.coeffs[0], p.buf.coeffs[0])
	vec.InvMFormTo(p.buf.coeffs[0], p.buf.coeffs[0], p.mod)

	pow5 := 1
	mask := p.packLen<<2 - 1
	revShiftBits := 64 - int(num.Log2(uint64(p.packLen))+1)
	for i := 0; i < p.packLen; i++ {
		idx := pow5 >> 1
		if pow5 > p.packLen<<1 {
			idx = p.packLen<<1 - 1 - idx
		}

		idxIn := int(bits.Reverse64(uint64(idx)) >> revShiftBits)
		if idxIn >= p.packLen {
			idxIn = p.packLen<<1 - 1 - idxIn
		}

		p.buf.coeffs[1][i] = p.buf.coeffs[0][idxIn]
		pow5 = (5 * pow5) & mask
	}

	copy(vOut, p.buf.coeffs[1][:vLen])
}

type autFixedPow2Mod3Packer struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int
	// nttRank is the rank of the NTT.
	nttRank int
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

// newAutFixedPow2Mod3Packer creates a new [autFixedPow2Mod3Packer].
func newAutFixedPow2Mod3Packer(params dft.RingParameters, mod *num.Modulus) *autFixedPow2Mod3Packer {
	primes, _ := num.Factor(mod.Value())
	packLen := params.CycloOrder() >> 2
	LogPackLen := bits.TrailingZeros64(primes[0]+1) - 1
	if packLen > (1 << LogPackLen) {
		packLen = 1 << LogPackLen
	}

	nttRank := packLen
	if int(primes[0])%params.CycloOrder() != params.CycloOrder()-1 {
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
	bitReverseInPlace(twLarge)
	bitReverseInPlace(twInvLarge)

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
			packIdx[i] = int(num.Exp(5, uint64(i), nil)) & mask
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
			idx1 := int(num.Exp(5, uint64(i), nil)) & mask
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

	return &autFixedPow2Mod3Packer{
		params: params,
		mod:    mod,

		packLen: packLen,
		nttRank: nttRank,
		packIdx: packIdx,

		tw:    tw,
		twInv: twInv,

		rankInv: num.Inv(uint64(nttRank<<1), mod),

		buf: newPow2Mod3PackerBuffer(nttRank),
	}
}

// Params returns the parameters of the packer.
func (p *autFixedPow2Mod3Packer) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *autFixedPow2Mod3Packer) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *autFixedPow2Mod3Packer) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *autFixedPow2Mod3Packer) SafeCopy() PackerInt {
	return &autFixedPow2Mod3Packer{
		params:  p.params,
		mod:     p.mod,
		packLen: p.packLen,

		packIdx: p.packIdx,

		tw:    p.tw,
		twInv: p.twInv,

		rankInv: p.rankInv,

		buf: newPow2Mod3PackerBuffer(p.nttRank),
	}
}

// Pack packs the input vector into a polynomial.
func (p *autFixedPow2Mod3Packer) Pack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.params.Rank())
	p.PackTo(vOut, vIn)
	return vOut
}

// PackTo packs the input vector into a polynomial.
func (p *autFixedPow2Mod3Packer) PackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vIn)

	for i := 0; i < p.packLen; i++ {
		idx := p.packIdx[i]
		p.buf.coeffs[idx].Real = vIn[i&(vLen-1)]
		p.buf.coeffs[idx].Imag = 0

		if p.packLen != p.nttRank {
			idx = p.packIdx[i+p.packLen]
			p.buf.coeffs[idx].Real = vIn[i&(vLen-1)]
			p.buf.coeffs[idx].Imag = 0
		}
	}

	invNTTGaloisRingInPlacePow2(p.buf.coeffs, p.twInv, p.mod)

	var u gnum.GaussianInt
	skip := p.params.Rank() / p.nttRank
	vOut[0] = num.Add(p.buf.coeffs[0].Real, p.buf.coeffs[0].Real, p.mod)
	for i := 1; i < p.nttRank; i++ {
		u = gnum.Mul(p.buf.coeffs[p.nttRank-i], p.tw[0], p.mod)
		vOut[i*skip] = num.Add(p.buf.coeffs[i].Real, u.Real, p.mod)
	}

	vec.ScalarMulTo(vOut, vOut, p.rankInv, p.mod)
}

func (p *autFixedPow2Mod3Packer) UnPack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, vIn)
	return vOut
}

func (p *autFixedPow2Mod3Packer) UnPackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vOut)

	skip := p.params.Rank() / p.nttRank
	p.buf.coeffs[0].Real = vIn[0]
	p.buf.coeffs[0].Imag = 0
	for i := 1; i < p.nttRank; i++ {
		real := num.Mul(p.tw[0].Real, vIn[(p.nttRank-i)*skip], p.mod)
		imag := num.Mul(p.tw[0].Imag, vIn[(p.nttRank-i)*skip], p.mod)
		p.buf.coeffs[i].Real = num.Sub(vIn[i*skip], real, p.mod)
		p.buf.coeffs[i].Imag = num.Neg(imag, p.mod)
	}

	nttGaloisRingInPlacePow2(p.buf.coeffs, p.tw, p.mod)

	for i := 0; i < vLen; i++ {
		idx := p.packIdx[i]
		vOut[i] = p.buf.coeffs[idx].Real
	}
}

// autFixedPrimePacker is a packer for autfixed ring.
type autFixedPrimePacker struct {
	params dft.RingParameters
	mod    *num.Modulus

	// packLen is the packing length.
	packLen int

	// resol is the resolution of unity.
	resol [][]uint64
	// invResol is the inverse resolution of unity.
	invResol [][]uint64

	// ambRank is the rank of the ambient NTT.
	ambRank int
	// ambModLen is the length of the ambient modulus.
	ambModLen int
	// ambMod is the ambient modulus.
	ambMod []*num.Modulus
	// ambNTT is the ambient NTT.
	ambNTT []dft.Transformer
	// embedder is the embedder for the packing.
	embedder *crt.Embedder

	buf packerBuffer
}

// newAutFixedPrimePacker creates a new [autFixedPacker].
func newAutFixedPrimePacker(params dft.RingParameters, mod *num.Modulus) *autFixedPrimePacker {
	cycloOrd := params.CycloOrder()

	// Compute the resolution of unity.
	primes, exps := num.Factor(mod.Value())
	if len(primes) != 1 {
		panic("newAutFixedPacker: mod must be a prime power")
	}
	prime, exp := primes[0], exps[0]
	resolution := FindResolutionOfUnity(cycloOrd, prime, exp)
	resolLen := len(resolution)
	invResolution := make([]uint64, resolLen)
	ord := num.Order(prime, num.NewModulus(cycloOrd))
	if ord&1 == 1 {
		copy(invResolution[resolLen/2:], resolution[:resolLen/2])
		copy(invResolution[:resolLen/2], resolution[resolLen/2:])
	} else {
		copy(invResolution, resolution)
	}

	vec.ScalarMulTo(invResolution, invResolution, num.Reduce(uint64(cycloOrd), mod), mod)
	vec.ScalarAddTo(invResolution, invResolution, num.Reduce(ord, mod), mod)

	// Reduce the resolution of unity to the packing length.
	packLen := num.GCD(params.Rank(), len(resolution))
	for i := 1; i < resolLen/packLen; i++ {
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
	maxBits := 2*num.Log2(mod.Value()) + num.Log2(ambRank) + 1
	ambModLen := int(math.Ceil(maxBits / num.MaxModulusBits))
	ambMod := dft.MustFindPrevNTTPrimes(ambParams, num.MaxModulusBits, ambModLen)
	ambNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		ambNTT[i] = dft.NewTransformer(ambParams, ambMod[i])
	}

	resol := make([][]uint64, ambModLen)
	invResol := make([][]uint64, ambModLen)
	for i := 0; i < ambModLen; i++ {
		resol[i] = make([]uint64, ambRank)
		invResol[i] = make([]uint64, ambRank)

		copy(resol[i][:packLen], resolution)
		copy(invResol[i][:packLen], invResolution)

		ambNTT[i].ForwardTo(resol[i], resol[i])
		ambNTT[i].ForwardTo(invResol[i], invResol[i])
	}

	embedder := crt.NewEmbedder([]*num.Modulus{mod}, ambMod)

	return &autFixedPrimePacker{
		params: params,
		mod:    mod,

		packLen: packLen,

		resol:    resol,
		invResol: invResol,

		ambRank:   ambRank,
		ambModLen: ambModLen,
		ambMod:    ambMod,
		ambNTT:    ambNTT,
		embedder:  embedder,

		buf: newPackerBuffer(ambModLen, ambRank),
	}
}

// FindResolutionOfUnity finds the resolution of unity for the given cyclotomic order, prime, and exponent.
// The algorithm is from https://eprint.iacr.org/2024/2032.
func FindResolutionOfUnity(cycloOrd int, prime uint64, exp uint64) []uint64 {
	if !num.IsPrime(cycloOrd) {
		panic("findResolutionOfUnity: cycloOrd must be a prime number")
	}
	if !num.IsPrime(prime) {
		panic("findResolutionOfUnity: prime must be a prime number")
	}

	cycloOrdMod := num.NewModulus(cycloOrd)
	ord := num.Order(prime, cycloOrdMod)
	rank := int(num.Totient(uint64(cycloOrd)) / ord)
	r := gr.NewGaloisRing(uint64(num.Exp(prime, exp, nil)), int(ord))
	root := grNthRoot(r, cycloOrd)

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
		r.ExpTo(rootPow, root, num.Exp(prime, uint64(i), cycloOrdMod))

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

	cycloPoly := dft.CyclotomicPolynomial(cycloOrd)
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
	params := dft.NewCyclicParameters(num.NextProdPower(cycloOrd, []int{2}))
	tf := crt.NewPolyEvaluator(params, []*num.Modulus{modulus})
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
	gen := num.Generators(cycloOrdMod)[0]
	idx := uint64(1)
	for i := 0; i < rank; i++ {
		resolution[i] = resolPoly.Coeffs[0][idx] + modulus.Value() - resolPoly.Coeffs[0][0]
		if resolution[i] >= modulus.Value() {
			resolution[i] -= modulus.Value()
		}
		idx = num.Mul(idx, gen, cycloOrdMod)
	}

	return resolution
}

// Params returns the parameters of the packer.
func (p *autFixedPrimePacker) Params() dft.RingParameters {
	return p.params
}

// Modulus returns the modulus of the packer.
func (p *autFixedPrimePacker) Modulus() *num.Modulus {
	return p.mod
}

// PackLen returns the packing length of the packer.
func (p *autFixedPrimePacker) PackLen() int {
	return p.packLen
}

// SafeCopy returns a safe copy of the packer.
func (p *autFixedPrimePacker) SafeCopy() PackerInt {
	return &autFixedPrimePacker{
		params:  p.params,
		mod:     p.mod,
		packLen: p.packLen,

		resol:    p.resol,
		invResol: p.invResol,

		ambModLen: p.ambModLen,
		ambMod:    p.ambMod,
		ambNTT:    p.ambNTT,
		embedder:  p.embedder.SafeCopy(),

		buf: newPackerBuffer(p.ambRank, p.ambModLen),
	}
}

// Pack packs the input vector into a polynomial.
func (p *autFixedPrimePacker) Pack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.params.Rank())
	p.PackTo(vOut, vIn)
	return vOut
}

// PackTo packs the input vector into a polynomial.
func (p *autFixedPrimePacker) PackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vIn)

	for i := 0; i < p.ambModLen; i++ {
		for j := 0; j < p.packLen/vLen; j++ {
			p.buf.coeffs[i][vLen*j] = vIn[0]
			copy(p.buf.coeffs[i][vLen*j+1:vLen*(j+1)], vIn[1:])
			slices.Reverse(p.buf.coeffs[i][vLen*j+1 : vLen*(j+1)])
		}
		clear(p.buf.coeffs[i][p.packLen:])

		p.ambNTT[i].ForwardTo(p.buf.coeffs[i], p.buf.coeffs[i])
		vec.MMulLazyTo(p.buf.coeffs[i], p.buf.coeffs[i], p.resol[i], p.ambMod[i])
		p.ambNTT[i].InverseTo(p.buf.coeffs[i], p.buf.coeffs[i])
	}

	p.embedder.EmbedVecTo(p.buf.coeffs[:1], p.buf.coeffs[:p.ambModLen])
	if p.packLen != p.ambRank {
		vec.AddTo(p.buf.coeffs[0][:p.packLen-1], p.buf.coeffs[0][:p.packLen-1], p.buf.coeffs[0][p.packLen:2*p.packLen-1], p.mod)
	}

	for i := 0; i < p.params.Rank()/p.packLen; i++ {
		copy(vOut[i*p.packLen:(i+1)*p.packLen], p.buf.coeffs[0][:p.packLen])
	}
}

// UnPack unpacks the polynomial into a vector.
func (p *autFixedPrimePacker) UnPack(vIn []uint64) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, vIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *autFixedPrimePacker) UnPackTo(vOut []uint64, vIn []uint64) {
	vLen := len(vOut)

	for i := 0; i < p.ambModLen; i++ {
		p.buf.coeffs[i][0] = vIn[0]
		copy(p.buf.coeffs[i][1:p.packLen], vIn[1:p.packLen])
		slices.Reverse(p.buf.coeffs[i][1:p.packLen])
		clear(p.buf.coeffs[i][p.packLen:])

		p.ambNTT[i].ForwardTo(p.buf.coeffs[i], p.buf.coeffs[i])
		vec.MMulLazyTo(p.buf.coeffs[i], p.buf.coeffs[i], p.invResol[i], p.ambMod[i])
		p.ambNTT[i].InverseTo(p.buf.coeffs[i], p.buf.coeffs[i])
	}
	p.embedder.EmbedVecTo(p.buf.coeffs[:1], p.buf.coeffs[:p.ambModLen])

	if p.packLen != p.ambRank {
		vec.AddTo(p.buf.coeffs[0][:p.packLen-1], p.buf.coeffs[0][:p.packLen-1], p.buf.coeffs[0][p.packLen:2*p.packLen-1], p.mod)
	}

	copy(vOut, p.buf.coeffs[0][:vLen])
}
