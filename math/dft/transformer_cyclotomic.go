package dft

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// pow2CyclotomicTransformer is a transformer for power-of-two cyclotomic ring.
type pow2CyclotomicTransformer struct {
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
	// rankInvS is the Shoup form of rankInv.
	rankInvS uint64
}

// newPow2CyclotomicTransformer creates a new [pow2CyclotomicTransformer].
func newPow2CyclotomicTransformer(params RingParameters, mod *num.Modulus) *pow2CyclotomicTransformer {
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

	rankInv := num.InvMForm(num.Inv(uint64(params.rank), mod), mod)

	return &pow2CyclotomicTransformer{
		params: params,
		mod:    mod,

		tw:     tw,
		twS:    vec.SForm(tw, mod),
		twInv:  twInv,
		twInvS: vec.SForm(twInv, mod),

		rankInv:  rankInv,
		rankInvS: num.SForm(rankInv, mod),
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *pow2CyclotomicTransformer) ForwardTo(vNTT, v []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	copy(vNTT, v)
	nttInPlacePow2(vNTT, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(vNTT, vNTT, ntt.mod)
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *pow2CyclotomicTransformer) InverseTo(v, vNTT []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	copy(v, vNTT)
	inttInPlacePow2(v, ntt.twInv, ntt.twInvS, ntt.mod.Value())
	vec.SMulScalarTo(v, v, ntt.rankInv, ntt.rankInvS, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *pow2CyclotomicTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *pow2CyclotomicTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

// anyCyclotomicTransformer is a transformer for arbitrary order cyclotomic ring.
type anyCyclotomicTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT  Transformer
	reducer *cyclotomicReducer

	// idx is the CRT mapping index.
	idx []int

	pool *sync.Pool
}

// newAnyCyclotomicTransformer creates a new [anyCyclotomicTransformer].
func newAnyCyclotomicTransformer(params RingParameters, mod *num.Modulus) *anyCyclotomicTransformer {
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

	return &anyCyclotomicTransformer{
		params: params,
		mod:    mod,

		ambNTT:  NewTransformer(NewCyclicParameters(params.cycloOrd), mod),
		reducer: newCyclotomicReducer(params, mod),

		idx: idx,

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, params.cycloOrd)
				return &v
			},
		},
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *anyCyclotomicTransformer) ForwardTo(vNTT, v []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	vBufPtr := ntt.pool.Get().(*[]uint64)
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	copy(vBuf, v)
	clear(vBuf[ntt.params.rank:])

	ntt.ambNTT.ForwardTo(vBuf, vBuf)

	for i := 0; i < ntt.params.rank; i++ {
		vNTT[i] = vBuf[ntt.idx[i]]
	}
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *anyCyclotomicTransformer) InverseTo(v, vNTT []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	vBufPtr := ntt.pool.Get().(*[]uint64)
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	clear(vBuf)
	for i := 0; i < ntt.params.rank; i++ {
		vBuf[ntt.idx[i]] = vNTT[i]
	}

	ntt.ambNTT.InverseTo(vBuf, vBuf)

	ntt.reducer.reduceTo(v, vBuf)
}

// Params returns the ring parameters.
func (ntt *anyCyclotomicTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *anyCyclotomicTransformer) Modulus() *num.Modulus {
	return ntt.mod
}
