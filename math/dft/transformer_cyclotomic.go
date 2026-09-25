package dft

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// pow2CyclotomicTransformer is a transformer for [typePow2Cyclotomic].
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
	root := num.Generators(mod.Value())

	tw := make([]uint64, params.rank)
	twInv := make([]uint64, params.rank)
	tw[0], twInv[0] = 1, 1
	if params.rank > 1 {
		tw[1] = num.NthRoot(uint64(params.cycloIdx), root, mod.Value())
		twInv[1] = num.Inv(tw[1], mod)
		for i := 2; i < params.rank; i++ {
			tw[i] = num.Mul(tw[i-1], tw[1], mod)
			twInv[i] = num.Mul(twInv[i-1], twInv[1], mod)
		}
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

func (ntt *pow2CyclotomicTransformer) forwardTo(vNTT, v []uint64) {
	copy(vNTT, v)
	fwdNTTInPlacePow2(vNTT, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(vNTT, vNTT, ntt.mod)
}

func (ntt *pow2CyclotomicTransformer) inverseTo(v, vNTT []uint64) {
	copy(v, vNTT)
	invNTTInPlacePow2(v, ntt.twInv, ntt.twInvS, ntt.mod.Value())
	vec.SMulTo(v, v, ntt.rankInv, ntt.rankInvS, ntt.mod)
}

// anyCyclotomicTransformer is a transformer for [typeAnyCyclotomic].
type anyCyclotomicTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT  *Transformer
	reducer *cyclotomicReducer

	// idx is the CRT mapping index.
	idx []int

	pool *pool.Pool[*[]uint64]
}

// newAnyCyclotomicTransformer creates a new [anyCyclotomicTransformer].
func newAnyCyclotomicTransformer(params RingParameters, mod *num.Modulus) *anyCyclotomicTransformer {
	primes, exps := num.Factor(params.cycloIdx)

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

	root := num.Generators(uint64(params.cycloIdx))
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
			idxOut = int(num.Mul(uint64(idxOut), num.Exp(root[j], uint64(idxDigits[j]), uint64(params.cycloIdx)), uint64(params.cycloIdx)))
		}
		idx[i] = idxOut
	}

	if num.IsProdPowerOf(params.cycloIdx, cyclicNTTFactors) {
		cycloIdxFactors := make([]int, len(cyclicNTTFactors))
		exps := make([]int, len(cyclicNTTFactors))

		t := params.cycloIdx
		for i, f := range cyclicNTTFactors {
			cycloIdxFactors[i] = 1
			for t%f == 0 {
				t /= f
				cycloIdxFactors[i] *= f
				exps[i]++
			}
		}

		for i := 0; i < params.rank; i++ {
			idxOut := idx[i]
			idx[i] = 0
			for j := range cycloIdxFactors {
				var t, r int
				if cycloIdxFactors[j] != 1 {
					t = int(num.Mul(uint64(idxOut), num.Inv(uint64(params.cycloIdx/cycloIdxFactors[j]), uint64(cycloIdxFactors[j])), uint64(cycloIdxFactors[j])))
				}

				for k := 0; k < exps[j]; k++ {
					r = r*cyclicNTTFactors[j] + (t % cyclicNTTFactors[j])
					t /= cyclicNTTFactors[j]
				}

				for k := 0; k < j; k++ {
					r *= cycloIdxFactors[k]
				}
				idx[i] += r
			}
		}
	}

	return &anyCyclotomicTransformer{
		params: params,
		mod:    mod,

		ambNTT:  NewTransformer(NewCyclicParameters(params.cycloIdx), mod),
		reducer: newCyclotomicReducer(params, mod),

		idx: idx,

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.cycloIdx)
			return &v
		}),
	}
}

func (ntt *anyCyclotomicTransformer) forwardTo(vNTT, v []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	copy(vBuf, v)
	clear(vBuf[ntt.params.rank:])

	ntt.ambNTT.ForwardTo(vBuf, vBuf)

	for i := 0; i < ntt.params.rank; i++ {
		vNTT[i] = vBuf[ntt.idx[i]]
	}
}

func (ntt *anyCyclotomicTransformer) inverseTo(v, vNTT []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	clear(vBuf)
	for i := 0; i < ntt.params.rank; i++ {
		vBuf[ntt.idx[i]] = vNTT[i]
	}

	ntt.ambNTT.InverseTo(vBuf, vBuf)
	ntt.reducer.reduceTo(v, vBuf)
}
