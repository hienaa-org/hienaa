package dft

import (
	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type autFixedPrimeTransformer struct {
	params RingParameters
	mod    *num.Modulus

	// root is the sum of powers of primitive root.
	// Pre-transformed for a fast convolution.
	root []uint64
	// rootInv is the sum of powers of inverse primitive root.
	// Pre-transformed for a fast convolution.
	rootInv []uint64

	// ambNTT is the ambient NTT for the fast convolution between roots and input vector.
	ambNTT *dftops.CyclicPow2Transformer

	// fold is cyclotomic order divided by rank.
	fold uint64
	// ambRankInv is the modular inverse of the rank of the ambient NTT in Montgomery form.
	ambRankInv uint64
	// cycloOrdInv is the modular inverse of the cyclotomic order.
	cycloOrdInv uint64

	buf transformerBuffer
}

func newAutFixedPrimeTransformer(ringParams RingParameters, mod *num.Modulus) *autFixedPrimeTransformer {
	cycloOrd := uint64(ringParams.cycloOrd)
	rank := uint64(ringParams.rank)
	fold := int((cycloOrd - 1) / rank)

	cycloOrdMod := num.NewModulus(cycloOrd)
	cycloRoot := num.Generators(cycloOrdMod)[0]
	modRoot := num.NthRoot(int(cycloOrd), num.Generators(mod), mod)

	var convLen int
	if num.IsPowerOfTwo(rank) {
		convLen = int(rank)
	} else {
		convLen = int(num.NextProdPower(2*rank-1, []uint64{2}))
	}

	ambNTT := dftops.NewCyclicPow2Transformer(convLen, mod)
	ambRankInv := num.MForm(num.Inv(uint64(convLen), mod), mod)

	modRootPowSum := make([]uint64, convLen)
	modRootPowInvSum := make([]uint64, convLen)

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

	buf := newTransformerBuffer(convLen)

	return &autFixedPrimeTransformer{
		params: ringParams,
		mod:    mod,

		root:    modRootPowSum,
		rootInv: modRootPowInvSum,

		ambNTT: ambNTT,

		fold:        uint64(fold),
		ambRankInv:  ambRankInv,
		cycloOrdInv: num.Inv(cycloOrd, mod),

		buf: buf,
	}
}

func (ntt *autFixedPrimeTransformer) Params() RingParameters {
	return ntt.params
}

func (ntt *autFixedPrimeTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *autFixedPrimeTransformer) ForwardInPlace(coeffs []uint64) {
	clear(ntt.buf.coeffs)

	ntt.buf.coeffs[0] = coeffs[0]
	for i := 1; i < ntt.params.rank; i++ {
		ntt.buf.coeffs[i] = coeffs[ntt.params.rank-i]
	}

	ntt.ambNTT.ForwardInPlace(ntt.buf.coeffs)
	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.root, ntt.mod)
	dftops.INTTInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.TwInv, ntt.ambNTT.TwInvS, ntt.mod.Value())
	vec.ScalarMMulTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.ambRankInv, ntt.mod)

	if num.IsPowerOfTwo(uint64(ntt.params.rank)) {
		copy(coeffs, ntt.buf.coeffs)
	} else {
		vec.AddTo(coeffs, ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank-1], ntt.mod)
		coeffs[ntt.params.rank-1] = ntt.buf.coeffs[ntt.params.rank-1]
	}
}

func (ntt *autFixedPrimeTransformer) InverseInPlace(coeffs []uint64) {
	clear(ntt.buf.coeffs)

	vec.InvMFormTo(coeffs, coeffs, ntt.mod)

	sumFold := coeffs[0]
	ntt.buf.coeffs[0] = coeffs[0]
	for i := 1; i < ntt.params.rank; i++ {
		ntt.buf.coeffs[i] = coeffs[ntt.params.rank-i]
		sumFold = num.Add(sumFold, ntt.buf.coeffs[i], ntt.mod)
	}

	ntt.ambNTT.ForwardInPlace(ntt.buf.coeffs)
	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.rootInv, ntt.mod)
	ntt.ambNTT.InverseInPlace(ntt.buf.coeffs)

	if !num.IsPowerOfTwo(uint64(ntt.params.rank)) {
		vec.AddTo(ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank-1], ntt.mod)
	}

	sumFold = num.Mul(sumFold, ntt.fold, ntt.mod)
	// vec.SubScalarTo(coeffs, ntt.buf.coeffs[:ntt.params.rank], sumFold, ntt.mod)
	// vec.ScalarMulTo(coeffs, coeffs, ntt.cycloOrdInv, ntt.mod)
	for i := 0; i < ntt.params.rank; i++ {
		coeffs[i] = num.Sub(ntt.buf.coeffs[i], sumFold, ntt.mod)
	}
	vec.ScalarMulTo(coeffs, coeffs, ntt.cycloOrdInv, ntt.mod)
}

func (ntt *autFixedPrimeTransformer) SafeCopy() Transformer {
	return &autFixedPrimeTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		root:    ntt.root,
		rootInv: ntt.rootInv,

		ambNTT: ntt.ambNTT,

		cycloOrdInv: ntt.cycloOrdInv,

		buf: ntt.buf,
	}
}
