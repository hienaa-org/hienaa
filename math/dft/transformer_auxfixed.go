package dft

import (
	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

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
		params: ringParams,
		mod:    mod,

		ambNTT: ambNTT,

		isPow2: isPow2,

		root:    modRootPowSum,
		rootInv: modRootPowInvSum,

		fold:        uint64(fold),
		ambRankInv:  ambRankInv,
		cycloOrdInv: num.Inv(cycloOrd, mod),

		buf: newTransformerBuffer(ambRank),
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
	ntt.ambNTT.InverseInPlace(ntt.buf.coeffs)
	vec.ScalarMulTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.ambRankInv, ntt.mod)

	if ntt.isPow2 {
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

	if !ntt.isPow2 {
		vec.AddTo(ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank-1], ntt.mod)
	}

	sumFold = num.Mul(sumFold, ntt.fold, ntt.mod)
	vec.ScalarSubTo(coeffs, ntt.buf.coeffs[:ntt.params.rank], sumFold, ntt.mod)
	vec.ScalarMulTo(coeffs, coeffs, ntt.cycloOrdInv, ntt.mod)
}

func (ntt *autFixedPrimeTransformer) SafeCopy() Transformer {
	return &autFixedPrimeTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		ambNTT: ntt.ambNTT,

		isPow2: ntt.isPow2,

		root:    ntt.root,
		rootInv: ntt.rootInv,

		cycloOrdInv: ntt.cycloOrdInv,

		buf: newTransformerBuffer(ntt.ambNTT.Rank),
	}
}
