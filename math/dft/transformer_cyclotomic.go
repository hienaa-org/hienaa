package dft

import (
	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclotomicPow2Transformer is a transformer for power-of-two cyclotomic NTT.
type cyclotomicPow2Transformer struct {
	params RingParameters
	mod    *num.Modulus

	// tw is the twiddle factor for NTT.
	tw []uint64
	// twS is the Shoup form of tw.
	twS []uint64
	// twInv is the twiddle factor for InvNTT.
	twInv []uint64
	// twInvS is the Shoup form of twInv.
	twInvS []uint64

	// rankInv is the modular inverse of the rank.
	rankInv uint64
}

// newCyclotomicPow2Transformer creates a new [pow2CyclotomicTransformer].
func newCyclotomicPow2Transformer(ringParams RingParameters, mod *num.Modulus) *cyclotomicPow2Transformer {
	root := num.Generators(mod)

	tw := make([]uint64, ringParams.rank)
	twInv := make([]uint64, ringParams.rank)
	tw[0], tw[1] = 1, num.NthRoot(ringParams.cycloOrd, root, mod)
	twInv[0], twInv[1] = 1, num.Inv(tw[1], mod)
	for i := 2; i < ringParams.rank; i++ {
		tw[i] = num.Mul(tw[i-1], tw[1], mod)
		twInv[i] = num.Mul(twInv[i-1], twInv[1], mod)
	}
	vec.RadixReverseInPlace(tw, 2)
	vec.RadixReverseInPlace(twInv, 2)

	twS := make([]uint64, ringParams.rank)
	twInvS := make([]uint64, ringParams.rank)
	for i := 0; i < ringParams.rank; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	rankInv := num.InvMForm(num.Inv(uint64(ringParams.rank), mod), mod)

	return &cyclotomicPow2Transformer{
		params: ringParams,
		mod:    mod,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: rankInv,
	}
}

func (ntt *cyclotomicPow2Transformer) ForwardInPlace(coeffs []uint64) {
	dftops.NTTInPlacePow2(coeffs, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(coeffs, coeffs, ntt.mod)
}

func (ntt *cyclotomicPow2Transformer) InverseInPlace(coeffs []uint64) {
	dftops.INTTInPlacePow2(coeffs, ntt.twInv, ntt.twInvS, ntt.mod.Value())
	vec.ScalarMulTo(coeffs, coeffs, ntt.rankInv, ntt.mod)
}

func (ntt *cyclotomicPow2Transformer) Params() RingParameters {
	return ntt.params
}

func (ntt *cyclotomicPow2Transformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *cyclotomicPow2Transformer) SafeCopy() Transformer {
	return ntt
}

// cyclotomicAnyTransformer is a transformer for power-of-two cyclotomic NTT.
type cyclotomicAnyTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT  Transformer
	reducer *dftops.CyclotomicReducerNTTModulus

	// idx is the CRT mapping index.
	idx []uint64

	buf transformerBuffer
}

// newCyclotomicAnyTransformer creates a new [cyclotomicAnyTransformer].
func newCyclotomicAnyTransformer(ringParams RingParameters, mod *num.Modulus) *cyclotomicAnyTransformer {
	cycloOrd := uint64(ringParams.cycloOrd)
	cycloOrdMod := num.NewModulus(cycloOrd)
	primes, exps := num.Factor(cycloOrd)

	dims := make([]uint64, len(primes))
	for i := range dims {
		pExp := uint64(1)
		for j := 0; j < int(exps[i]); j++ {
			pExp *= primes[i]
		}
		dims[i] = pExp - pExp/primes[i]
	}

	root := num.GeneratorsWithFactors(cycloOrdMod, primes, exps)
	idx := make([]uint64, ringParams.rank)

	for i := 0; i < ringParams.rank; i++ {
		idxIn := uint64(i)
		idxDigits := make([]uint64, len(dims))
		for j := 0; j < len(dims); j++ {
			idxDigits[j] = idxIn % dims[j]
			idxIn /= dims[j]
		}
		idxOut := uint64(1)
		for j := 0; j < len(dims); j++ {
			idxOut = num.Mul(idxOut, num.Exp(root[j], idxDigits[j], cycloOrdMod), cycloOrdMod)
		}
		idx[uint64(i)] = idxOut
	}

	return &cyclotomicAnyTransformer{
		params: ringParams,
		mod:    mod,

		ambNTT:  NewTransformer(NewCyclicParameters(ringParams.cycloOrd), mod),
		reducer: dftops.NewCyclotomicReducerNTTModulus(ringParams.cycloOrd, ringParams.rank, mod),

		idx: idx,

		buf: newTransformerBuffer(ringParams.cycloOrd),
	}
}

func (ntt *cyclotomicAnyTransformer) ForwardInPlace(coeffs []uint64) {
	copy(ntt.buf.coeffs, coeffs)
	clear(ntt.buf.coeffs[ntt.params.rank:])

	ntt.ambNTT.ForwardInPlace(ntt.buf.coeffs)

	for i := 0; i < ntt.params.rank; i++ {
		coeffs[i] = ntt.buf.coeffs[ntt.idx[i]]
	}
}

func (ntt *cyclotomicAnyTransformer) InverseInPlace(coeffs []uint64) {
	clear(ntt.buf.coeffs)
	for i := 0; i < ntt.params.rank; i++ {
		ntt.buf.coeffs[ntt.idx[i]] = coeffs[i]
	}

	ntt.ambNTT.InverseInPlace(ntt.buf.coeffs)

	ntt.reducer.ReduceTo(ntt.buf.coeffs, ntt.buf.coeffs)
	copy(coeffs, ntt.buf.coeffs)
}

func (ntt *cyclotomicAnyTransformer) Params() RingParameters {
	return ntt.params
}

func (ntt *cyclotomicAnyTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *cyclotomicAnyTransformer) SafeCopy() Transformer {
	return &cyclotomicAnyTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		ambNTT:  ntt.ambNTT.SafeCopy(),
		reducer: ntt.reducer.SafeCopy(),

		idx: ntt.idx,

		buf: newTransformerBuffer(ntt.params.cycloOrd),
	}
}
