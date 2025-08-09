package dft

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// autFixedPow2Transformer is a transformer for power-of-two conjugate invariant ring.
type autFixedPow2Transformer struct {
	params RingParameters
	mod    *num.Modulus

	// tw is the twiddle factor for NTT.
	tw []uint64
	// twInv is the twiddle factor for InvNTT.
	twInv []uint64
	// twS is the Shoup form of tw.
	twS []uint64
	// twInvS is the Shoup form of twInv.
	twInvS []uint64

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	buf transformerBuffer
}

// newAutFixedPow2Transformer creates a new [autFixedPow2Transformer].
func newAutFixedPow2Transformer(params RingParameters, mod *num.Modulus) *autFixedPow2Transformer {
	root := num.Generators(mod)

	twoRank := params.rank << 1
	twLarge := make([]uint64, twoRank)
	twInvLarge := make([]uint64, twoRank)
	twLarge[0], twLarge[1] = 1, num.NthRoot(params.cycloOrd, root, mod)
	twInvLarge[0], twInvLarge[1] = 1, num.Inv(twLarge[1], mod)
	for i := 2; i < twoRank; i++ {
		twLarge[i] = num.Mul(twLarge[i-1], twLarge[1], mod)
		twInvLarge[i] = num.Mul(twInvLarge[i-1], twInvLarge[1], mod)
	}
	vec.RadixReverseInPlace(twLarge, 2)
	vec.RadixReverseInPlace(twInvLarge, 2)

	tw := make([]uint64, params.rank)
	twInv := make([]uint64, params.rank)

	tw[0] = twLarge[1]
	twInv[0] = twInvLarge[1]
	for m := 1; m <= params.rank/2; m <<= 1 {
		copy(tw[m:2*m], twLarge[2*m:3*m])
		copy(twInv[m:2*m], twInvLarge[2*m:3*m])
	}

	twS := make([]uint64, params.rank)
	twInvS := make([]uint64, params.rank)
	for i := 0; i < params.rank; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	return &autFixedPow2Transformer{
		params: params,
		mod:    mod,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: num.InvMForm(num.Inv(uint64(twoRank), mod), mod),

		buf: newTransformerBuffer(params.rank),
	}
}

func (ntt *autFixedPow2Transformer) ForwardInPlace(coeffs []uint64) {
	copy(ntt.buf.coeffs, coeffs)
	slices.Reverse(coeffs[1:])
	vec.ScalarMulSubLazyTo(ntt.buf.coeffs[1:], coeffs[1:], ntt.tw[0], ntt.mod)

	dftops.NTTInPlacePow2(ntt.buf.coeffs, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(coeffs, ntt.buf.coeffs, ntt.mod)
}

func (ntt *autFixedPow2Transformer) InverseInPlace(coeffs []uint64) {
	dftops.INTTInPlacePow2(coeffs, ntt.twInv, ntt.twInvS, ntt.mod.Value())

	copy(ntt.buf.coeffs, coeffs)
	slices.Reverse(coeffs[1:])
	ntt.buf.coeffs[0] = coeffs[0] + coeffs[0]
	vec.ScalarMulAddLazyTo(ntt.buf.coeffs[1:], coeffs[1:], ntt.tw[0], ntt.mod)

	vec.ScalarMulTo(coeffs, ntt.buf.coeffs, ntt.rankInv, ntt.mod)
}

func (ntt *autFixedPow2Transformer) Params() RingParameters {
	return ntt.params
}

func (ntt *autFixedPow2Transformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *autFixedPow2Transformer) SafeCopy() Transformer {
	return &autFixedPow2Transformer{
		params: ntt.params,
		mod:    ntt.mod,

		tw:     ntt.tw,
		twS:    ntt.twS,
		twInv:  ntt.twInv,
		twInvS: ntt.twInvS,

		rankInv: ntt.rankInv,

		buf: newTransformerBuffer(ntt.params.rank),
	}
}

// autFixedPrimeTransformer is a transformer for prime order decomposition ring.
type autFixedPrimeTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT *dftops.CyclicPow2Transformer

	// isPow2 is true of the rank is power-of-two.
	isPow2 bool

	// root is the sum of powers of primitive root.
	// Pre-transformed for a fast convolution.
	root []uint64
	// rootInv is the sum of powers of inverse primitive root.
	// Pre-transformed for a fast convolution.
	rootInv []uint64

	// fold is cyclotomic order divided by rank.
	fold uint64
	// ambRankInvM is the modular inverse of the rank of the ambient NTT in Montgomery form.
	ambRankInvM uint64
	// cycloOrdInv is the modular inverse of the cyclotomic order.
	cycloOrdInv uint64

	buf transformerBuffer
}

// newAutFixedPrimeTransformer creates a new [autFixedPrimeTransformer].
func newAutFixedPrimeTransformer(params RingParameters, mod *num.Modulus) *autFixedPrimeTransformer {
	cycloOrd := uint64(params.cycloOrd)
	rank := uint64(params.rank)
	fold := int((cycloOrd - 1) / rank)

	isPow2 := num.IsPowerOfTwo(rank)

	cycloOrdMod := num.NewModulus(cycloOrd)
	cycloRoot := num.Generators(cycloOrdMod)[0]
	modRoot := num.NthRoot(int(cycloOrd), num.Generators(mod), mod)

	var ambRank int
	if isPow2 {
		ambRank = int(rank)
	} else {
		ambRank = int(num.NextProdPower(2*rank-1, []uint64{2}))
	}
	ambNTT := dftops.NewCyclicPow2Transformer(ambRank, mod)
	ambRankInv := num.MForm(num.Inv(uint64(ambRank), mod), mod)

	modRootPowSum := make([]uint64, ambRank)
	modRootPowInvSum := make([]uint64, ambRank)

	cycloRootPowRank := num.Exp(cycloRoot, rank, cycloOrdMod)
	modRootPow := modRoot
	modRootPowInv := num.Inv(modRoot, mod)

	var modRootPowNext, modRootPowInvNext uint64
	for i := 0; i < fold; i++ {
		modRootPow = num.Exp(modRootPow, cycloRootPowRank, mod)
		modRootPowInv = num.Exp(modRootPowInv, cycloRootPowRank, mod)

		modRootPowNext = modRootPow
		modRootPowInvNext = modRootPowInv

		for j := 0; j < int(rank); j++ {
			modRootPowSum[j] = num.Add(modRootPowSum[j], modRootPowNext, mod)
			modRootPowInvSum[j] = num.Add(modRootPowInvSum[j], modRootPowInvNext, mod)

			modRootPowNext = num.Exp(modRootPowNext, cycloRoot, mod)
			modRootPowInvNext = num.Exp(modRootPowInvNext, cycloRoot, mod)
		}
	}

	ambNTT.ForwardInPlace(modRootPowSum)
	ambNTT.ForwardInPlace(modRootPowInvSum)

	return &autFixedPrimeTransformer{
		params: params,
		mod:    mod,

		ambNTT: ambNTT,

		isPow2: isPow2,

		root:    modRootPowSum,
		rootInv: modRootPowInvSum,

		fold:        uint64(fold),
		ambRankInvM: ambRankInv,
		cycloOrdInv: num.Inv(cycloOrd, mod),

		buf: newTransformerBuffer(ambRank),
	}
}

func (ntt *autFixedPrimeTransformer) ForwardInPlace(coeffs []uint64) {
	clear(ntt.buf.coeffs)

	ntt.buf.coeffs[0] = coeffs[0]
	for i := 1; i < ntt.params.rank; i++ {
		ntt.buf.coeffs[i] = coeffs[ntt.params.rank-i]
	}

	dftops.NTTInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.Tw, ntt.ambNTT.TwS, ntt.mod.Value())
	vec.MFormTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.mod)

	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.root, ntt.mod)

	dftops.INTTInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.TwInv, ntt.ambNTT.TwInvS, ntt.mod.Value())
	vec.ScalarMMulTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.ambRankInvM, ntt.mod)

	if ntt.isPow2 {
		copy(coeffs, ntt.buf.coeffs)
	} else {
		vec.AddTo(coeffs, ntt.buf.coeffs[:ntt.params.rank], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank], ntt.mod)
	}
}

func (ntt *autFixedPrimeTransformer) InverseInPlace(coeffs []uint64) {
	clear(ntt.buf.coeffs)

	sumFold := coeffs[0]
	ntt.buf.coeffs[0] = coeffs[0]
	for i := 1; i < ntt.params.rank; i++ {
		ntt.buf.coeffs[i] = coeffs[ntt.params.rank-i]
		sumFold = num.Add(sumFold, ntt.buf.coeffs[i], ntt.mod)
	}

	dftops.NTTInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.Tw, ntt.ambNTT.TwS, ntt.mod.Value())

	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.rootInv, ntt.mod)

	dftops.INTTInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.TwInv, ntt.ambNTT.TwInvS, ntt.mod.Value())
	vec.ScalarMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.ambNTT.RankInv, ntt.mod)

	if !ntt.isPow2 {
		vec.AddTo(ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank-1], ntt.mod)
	}

	sumFold = num.MMul(sumFold, ntt.fold, ntt.mod)
	vec.ScalarSubTo(coeffs, ntt.buf.coeffs[:ntt.params.rank], sumFold, ntt.mod)
	vec.ScalarMulTo(coeffs, coeffs, ntt.cycloOrdInv, ntt.mod)
}

func (ntt *autFixedPrimeTransformer) Params() RingParameters {
	return ntt.params
}

func (ntt *autFixedPrimeTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *autFixedPrimeTransformer) SafeCopy() Transformer {
	return &autFixedPrimeTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		ambNTT: ntt.ambNTT,

		isPow2: ntt.isPow2,

		root:    ntt.root,
		rootInv: ntt.rootInv,

		fold:        ntt.fold,
		ambRankInvM: ntt.ambRankInvM,
		cycloOrdInv: ntt.cycloOrdInv,

		buf: newTransformerBuffer(ntt.ambNTT.Rank),
	}
}
