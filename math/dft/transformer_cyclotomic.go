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
func newCyclotomicPow2Transformer(params RingParameters, mod *num.Modulus) *cyclotomicPow2Transformer {
	root := num.Generators(mod)

	tw := make([]uint64, params.rank)
	twInv := make([]uint64, params.rank)
	tw[0], tw[1] = 1, num.NthRoot(params.cycloOrd, root, mod)
	twInv[0], twInv[1] = 1, num.Inv(tw[1], mod)
	for i := 2; i < params.rank; i++ {
		tw[i] = num.Mul(tw[i-1], tw[1], mod)
		twInv[i] = num.Mul(twInv[i-1], twInv[1], mod)
	}
	vec.RadixReverseInPlace(tw, 2)
	vec.RadixReverseInPlace(twInv, 2)

	twS := make([]uint64, params.rank)
	twInvS := make([]uint64, params.rank)
	for i := 0; i < params.rank; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	rankInv := num.InvMForm(num.Inv(uint64(params.rank), mod), mod)

	return &cyclotomicPow2Transformer{
		params: params,
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
	idx []int

	buf transformerBuffer
}

// newCyclotomicAnyTransformer creates a new [cyclotomicAnyTransformer].
func newCyclotomicAnyTransformer(params RingParameters, mod *num.Modulus) *cyclotomicAnyTransformer {
	cycloOrdMod := num.NewModulus(params.cycloOrd)
	primes, exps := num.Factor(params.cycloOrd)

	dims := make([]int, len(primes))
	for i := range dims {
		pExp := 1
		for j := 0; j < int(exps[i]); j++ {
			pExp *= primes[i]
		}
		dims[i] = pExp - pExp/primes[i]
	}
	if primes[0] == 2 {
		if exps[0] == 1 {
			dims = dims[1:]
		} else if exps[0] > 2 {
			dims = append([]int{0}, dims...)
			dims[0], dims[1] = dims[1]/2, 2
		}
	}

	root := num.GeneratorsWithFactors(cycloOrdMod, vec.Cast[uint64](primes), vec.Cast[uint64](exps))
	idx := make([]int, params.rank)

	idxDigits := make([]int, len(dims))
	for i := 0; i < params.rank; i++ {
		idxIn := i
		for j := 0; j < len(dims); j++ {
			idxDigits[j] = idxIn % dims[j]
			idxIn /= dims[j]
		}
		idxOut := 1
		for j := 0; j < len(dims); j++ {
			idxOut = int(num.Mul(uint64(idxOut), num.Exp(root[j], uint64(idxDigits[j]), cycloOrdMod), cycloOrdMod))
		}
		idx[i] = idxOut
	}

	if num.IsProdPowerOf(params.cycloOrd, cyclicNTTFactors) {
		cycloOrdFactors := make([]int, len(cyclicNTTFactors))
		cycloOrdFactorsMod := make([]*num.Modulus, len(cyclicNTTFactors))
		exps := make([]int, len(cyclicNTTFactors))

		t := params.cycloOrd
		for i, f := range cyclicNTTFactors {
			cycloOrdFactors[i] = 1
			for t%f == 0 {
				t /= f
				cycloOrdFactors[i] *= f
				exps[i]++
			}

			if cycloOrdFactors[i] != 1 {
				cycloOrdFactorsMod[i] = num.NewModulus(cycloOrdFactors[i])
			}
		}

		for i := 0; i < params.rank; i++ {
			idxOut := idx[i]
			idx[i] = 0
			for j := range cycloOrdFactors {
				var t, r int
				if cycloOrdFactors[j] != 1 {
					t = int(num.Mul(uint64(idxOut), num.Inv(uint64(params.cycloOrd/cycloOrdFactors[j]), cycloOrdFactorsMod[j]), cycloOrdFactorsMod[j]))
				}

				for d := 0; d < int(exps[j]); d++ {
					r = r*cyclicNTTFactors[j] + (t % cyclicNTTFactors[j])
					t /= cyclicNTTFactors[j]
				}

				for k := 0; k < j; k++ {
					r *= cycloOrdFactors[k]
				}
				idx[i] += r
			}
		}
	}

	return &cyclotomicAnyTransformer{
		params: params,
		mod:    mod,

		ambNTT:  NewTransformer(NewCyclicParameters(params.cycloOrd), mod),
		reducer: dftops.NewCyclotomicReducerNTTModulus(params.cycloOrd, params.rank, mod),

		idx: idx,

		buf: newTransformerBuffer(params.cycloOrd),
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
