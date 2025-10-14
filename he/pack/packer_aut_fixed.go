package pack

import (
	"math"
	"math/big"
	"slices"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/gr"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// autFixedPrimePacker is a packer for autfixed ring.
type autFixedPrimePacker struct {
	params  dft.RingParameters
	mod     *num.Modulus
	packLen int

	resol    [][]uint64
	invResol [][]uint64

	ambRank   int
	ambModLen int
	ambMod    []*num.Modulus
	ambNTT    []dft.Transformer
	embedder  *crt.Embedder

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
	packLen := int(num.GCD(uint64(params.Rank()), uint64(len(resolution))))
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
	ambMod := dft.FindPrevNTTPrimes(ambParams, num.MaxModulusBits, ambModLen)
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
		params:  params,
		mod:     mod,
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

	// Sampler.
	// Use the same seed nil to fix the randomness.
	rSrc := csprng.NewUniformSamplerWithSeed(nil)

	// genEl is the generator of the multiplicative group of the Galois ring.
	// root is cycloOrd-th root of unity of the multiplicative group.
	genEl := r.NewElement()
	root := r.NewElement()
	one := r.NewElementFromUint64(1)
	ordOverCycloOrd := new(big.Int).Div(r.Ord(), big.NewInt(int64(cycloOrd)))
	tmpEl := r.NewElement()
	genCoeffs := genEl.Coeffs()
	for {
		for i := 0; i < int(ord); i++ {
			genCoeffs[i] = rSrc.SampleN(r.Modulus())
		}
		r.ExpBigTo(tmpEl, genEl, r.Ord())
		if tmpEl.IsEqual(one) {
			r.ExpBigTo(root, genEl, ordOverCycloOrd)
			if !root.IsEqual(one) {
				break
			}
		}
	}

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

	tf.NTTTo(resolPoly, resolPoly)
	tf.NTTTo(remInvPoly, remInvPoly)
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
func (p *autFixedPrimePacker) Pack(vIn []uint64) *crt.Poly {
	pOut := crt.NewPoly(p.params.Rank(), 1)
	p.PackTo(pOut, vIn)
	return pOut
}

// PackTo packs the input vector into a polynomial.
func (p *autFixedPrimePacker) PackTo(pOut *crt.Poly, vIn []uint64) {
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
		copy(pOut.Coeffs[0][i*p.packLen:(i+1)*p.packLen], p.buf.coeffs[0][:p.packLen])
	}
}

// UnPack unpacks the polynomial into a vector.
func (p *autFixedPrimePacker) UnPack(pIn *crt.Poly) []uint64 {
	vOut := make([]uint64, p.packLen)
	p.UnPackTo(vOut, pIn)
	return vOut
}

// UnPackTo unpacks the polynomial into a vector.
func (p *autFixedPrimePacker) UnPackTo(vOut []uint64, pIn *crt.Poly) {
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

	for i := 0; i < p.ambModLen; i++ {
		p.buf.coeffs[i][0] = pIn.Coeffs[0][0]
		copy(p.buf.coeffs[i][1:p.packLen], pIn.Coeffs[0][1:p.packLen])
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
