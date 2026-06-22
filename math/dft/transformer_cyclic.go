package dft

import (
	"slices"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

var (
	// cyclicNTTFactors are the factors of the rank for native cyclic NTT.
	cyclicNTTFactors = []int{2, 3, 5}
)

// cyclicTwiddleFactor computes the (inverse) twiddle factor for Cyclic (Inv)NTT.
func cyclicTwiddleFactor(rank, radix int, root []uint64, mod *num.Modulus) (tw, twInv []uint64) {
	if rank == 1 {
		return []uint64{1}, []uint64{1}
	}

	z := num.NthRoot(rank, root, mod)
	zInv := num.Inv(z, mod)

	tw = make([]uint64, rank)
	twInv = make([]uint64, rank)
	tw[0] = 1
	twInv[0] = 1
	for i := 1; i < rank/radix; i++ {
		tw[i] = num.Mul(tw[i-1], z, mod)
		twInv[i] = num.Mul(twInv[i-1], zInv, mod)
	}

	vec.RadixReverseInPlace(tw[:rank/radix], radix)
	vec.RadixReverseInPlace(twInv[:rank/radix], radix)

	t := rank / radix
	for i := 0; i < t; i++ {
		tw[i+t] = tw[i]
		twInv[i+t] = twInv[i]
		for j := 2 * t; j < radix*t; j += t {
			tw[i+j] = num.Mul(tw[i+j-t], tw[i], mod)
			twInv[i+j] = num.Mul(twInv[i+j-t], twInv[i], mod)
		}
	}

	for t > 1 {
		t /= radix
		for i := 0; i < t; i++ {
			tw[t+i] = tw[radix*t+i]
			twInv[t+i] = twInv[radix*t+i]
		}
		for j := 2; j < radix; j++ {
			for i := 0; i < t; i++ {
				tw[j*t+i] = num.Mul(tw[t*(j-1)+i], tw[t+i], mod)
				twInv[j*t+i] = num.Mul(twInv[t*(j-1)+i], twInv[t+i], mod)
			}
		}
	}

	return tw, twInv
}

// pow235CyclicTransformer is a transformer for cyclic ring with ranks multiple of 2, 3 and 5.
type pow235CyclicTransformer struct {
	params RingParameters
	mod    *num.Modulus

	// rankFactors are the factors of the rank.
	rankFactors []int

	// tw is the twiddle factor for NTT.
	// Ordered as [cyclicNTTFactors][Rank].
	tw [][]uint64
	// twM is the MulForm of tw.
	// Ordered as [cyclicNTTFactors][Rank].
	twM []vec.MulForm
	// twInv is the twiddle factor for InvNTT.
	// Ordered as [cyclicNTTFactors][Rank].
	twInv [][]uint64
	// twInvM is the MulForm of twInv.
	// Ordered as [cyclicNTTFactors][Rank].
	twInvM []vec.MulForm

	// root is the powers of primitive roots.
	// Ordered as [cyclicNTTFactors][Radix].
	root [][]uint64
	// rootM is the MulForm of root.
	// Ordered as [cyclicNTTFactors][Radix].
	rootM []vec.MulForm

	// rankInv is the modular inverse of the rank.
	rankInv uint64
	// rankInvM is the MulForm of rankInv.
	rankInvM num.MulForm

	// idx is the CRT mapping index.
	// Empty if the mapping is not needed, or in other words, rank is a prime power.
	idx []int

	pool *pool.Pool[*[]uint64]
}

// newCyclicPow235Transformer creates a new [cyclicNativeTransformer].
func newCyclicPow235Transformer(params RingParameters, mod *num.Modulus) *pow235CyclicTransformer {
	rankFactors := make([]int, len(cyclicNTTFactors))
	rankTmp := params.rank
	for i, f := range cyclicNTTFactors {
		rankFactors[i] = 1
		for rankTmp%f == 0 {
			rankTmp /= f
			rankFactors[i] *= f
		}
	}

	root := num.Generators(mod)

	tw := make([][]uint64, len(cyclicNTTFactors))
	twM := make([]vec.MulForm, len(cyclicNTTFactors))
	twInv := make([][]uint64, len(cyclicNTTFactors))
	twInvM := make([]vec.MulForm, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		tw[i], twInv[i] = cyclicTwiddleFactor(rankFactors[i], r, root, mod)
		twM[i] = vec.ToMulForm(tw[i], mod)
		twInvM[i] = vec.ToMulForm(twInv[i], mod)
	}

	rootExp := make([][]uint64, len(cyclicNTTFactors))
	rootExpM := make([]vec.MulForm, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		if rankFactors[i] == 1 {
			continue
		}
		rootExp[i] = make([]uint64, r)
		rootExp[i][0] = 1
		rootExp[i][1] = num.NthRoot(int(r), root, mod)
		for j := 2; j < int(r); j++ {
			rootExp[i][j] = num.Mul(rootExp[i][j-1], rootExp[i][1], mod)
		}
		rootExpM[i] = vec.ToMulForm(rootExp[i], mod)
	}

	rankInv := num.Inv(uint64(params.rank), mod)

	var idx []int
	if !slices.Contains(rankFactors, params.rank) {
		idx = make([]int, params.rank)
		for i2 := 0; i2 < rankFactors[0]; i2++ {
			for i3 := 0; i3 < rankFactors[1]; i3++ {
				for i5 := 0; i5 < rankFactors[2]; i5++ {
					idxIn := i2 + i3*rankFactors[0] + i5*rankFactors[0]*rankFactors[1]
					idxOut := i2*(params.rank/rankFactors[0]) + i3*(params.rank/rankFactors[1]) + i5*(params.rank/rankFactors[2])
					idx[idxIn%params.rank] = idxOut % params.rank
				}
			}
		}
	}

	return &pow235CyclicTransformer{
		params:      params,
		mod:         mod,
		rankFactors: rankFactors,

		tw:     tw,
		twM:    twM,
		twInv:  twInv,
		twInvM: twInvM,

		root:  rootExp,
		rootM: rootExpM,

		rankInv:  rankInv,
		rankInvM: num.ToMulForm(rankInv, mod),

		idx: idx,

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.rank)
			return &v
		}),
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *pow235CyclicTransformer) ForwardTo(vNTT, v []uint64) {
	if len(ntt.idx) > 0 {
		vBufPtr := ntt.pool.Get()
		vBuf := *vBufPtr
		defer ntt.pool.Put(vBufPtr)

		for i := 0; i < ntt.params.rank; i++ {
			vBuf[i] = v[ntt.idx[i]]
		}
		copy(vNTT, vBuf)
	} else {
		copy(vNTT, v)
	}

	if ntt.rankFactors[0] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] {
			fwdNTTInPlacePow2(vNTT[i:i+ntt.rankFactors[0]], ntt.tw[0], ntt.twM[0].SForm, ntt.mod.Value())
		}
	}

	if ntt.rankFactors[1] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] * ntt.rankFactors[1] {
			fwdNTTInPlacePow3(ntt.rankFactors[0], vNTT[i:i+ntt.rankFactors[0]*ntt.rankFactors[1]], ntt.tw[1], ntt.twM[1].SForm, ntt.root[1], ntt.rootM[1].SForm, ntt.mod.Value())
		}
	}

	if ntt.rankFactors[2] > 1 {
		fwdNTTInPlacePow5(ntt.rankFactors[0]*ntt.rankFactors[1], vNTT, ntt.tw[2], ntt.twM[2].SForm, ntt.root[2], ntt.rootM[2].SForm, ntt.mod.Value())
	}
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *pow235CyclicTransformer) InverseTo(v, vNTT []uint64) {
	copy(v, vNTT)

	if ntt.rankFactors[0] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] {
			invNTTInPlacePow2(v[i:i+ntt.rankFactors[0]], ntt.twInv[0], ntt.twInvM[0].SForm, ntt.mod.Value())
		}
	}

	if ntt.rankFactors[1] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] * ntt.rankFactors[1] {
			invNTTInPlacePow3(ntt.rankFactors[0], v[i:i+ntt.rankFactors[0]*ntt.rankFactors[1]], ntt.twInv[1], ntt.twInvM[1].SForm, ntt.root[1], ntt.rootM[1].SForm, ntt.mod.Value())
		}
	}

	if ntt.rankFactors[2] > 1 {
		invNTTInPlacePow5(ntt.rankFactors[0]*ntt.rankFactors[1], v, ntt.twInv[2], ntt.twInvM[2].SForm, ntt.root[2], ntt.rootM[2].SForm, ntt.mod.Value())
	}

	if len(ntt.idx) > 0 {
		vBufPtr := ntt.pool.Get()
		vBuf := *vBufPtr
		defer ntt.pool.Put(vBufPtr)

		for i := 0; i < ntt.params.rank; i++ {
			vBuf[ntt.idx[i]] = v[i]
		}
		copy(v, vBuf)
	}

	vec.FMulScalarTo(v, v, ntt.rankInv, ntt.rankInvM, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *pow235CyclicTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *pow235CyclicTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

// anyCyclicTransformer is a transformer for aribtrary rank cyclic ring.
// Internally, it uses Bluestein NTT.
type anyCyclicTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT *pow235CyclicTransformer

	// z is the factor for Z-transform.
	z []uint64
	// zM is the MulForm of z.
	zM vec.MulForm
	// zInv is the inverse factor for Z-transform.
	zInv []uint64
	// zInvM is the MulForm of zInv.
	zInvM vec.MulForm

	// chirp is the chirp factor for Bluestein NTT in Montgomery form.
	chirp []uint64
	// chirpM is the MulForm of chirpM.
	chirpM vec.MulForm
	// chirpInv is the inverse chirp factor for Bluestein NTT.
	chirpInv []uint64
	// chirpInvM is the MulForm of chirpInv.
	chirpInvM vec.MulForm

	pool *pool.Pool[*[]uint64]
}

// newAnyCyclicTransformer creates a new [anyCyclicTransformer].
func newAnyCyclicTransformer(params RingParameters, mod *num.Modulus) *anyCyclicTransformer {
	ambRank := num.NextProdPower(2*params.rank-1, []int{2})

	root := num.Generators(mod)
	zz := num.NthRoot(2*params.rank, root, mod)
	zzInv := num.Inv(zz, mod)

	z := make([]uint64, params.rank)
	zInv := make([]uint64, params.rank)
	for i := 0; i < params.rank; i++ {
		idx := (i * i) % (2 * params.rank)
		z[i] = num.Exp(zz, uint64(idx), mod)
		zInv[i] = num.Exp(zzInv, uint64(idx), mod)
	}

	ambNTT := newCyclicPow235Transformer(NewCyclicParameters(ambRank), mod)

	chirp := make([]uint64, ambRank)
	copy(chirp, zInv)
	copy(chirp[ambRank-params.rank+1:], zInv[1:])
	slices.Reverse(chirp[ambRank-params.rank+1:])

	vec.MulScalarTo(chirp, chirp, num.Inv(uint64(ambRank), mod), mod)
	fwdNTTInPlacePow2(chirp, ambNTT.tw[0], ambNTT.twM[0].SForm, mod.Value())
	vec.ReduceTo(chirp, chirp, mod)

	chirpInv := make([]uint64, ambRank)
	copy(chirpInv, z)
	copy(chirpInv[ambRank-params.rank+1:], z[1:])
	slices.Reverse(chirpInv[ambRank-params.rank+1:])

	vec.MulScalarTo(chirpInv, chirpInv, num.Inv(uint64(ambRank*params.rank), mod), mod)
	fwdNTTInPlacePow2(chirpInv, ambNTT.tw[0], ambNTT.twM[0].SForm, mod.Value())
	vec.ReduceTo(chirpInv, chirpInv, mod)

	return &anyCyclicTransformer{
		params: params,
		mod:    mod,

		ambNTT: ambNTT,

		z:     z,
		zM:    vec.ToMulForm(z, mod),
		zInv:  zInv,
		zInvM: vec.ToMulForm(zInv, mod),

		chirp:     chirp,
		chirpM:    vec.ToMulForm(chirp, mod),
		chirpInv:  chirpInv,
		chirpInvM: vec.ToMulForm(chirpInv, mod),

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, ambRank)
			return &v
		}),
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *anyCyclicTransformer) ForwardTo(vNTT, v []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	vec.FMulTo(vBuf[:ntt.params.rank], v, ntt.z, ntt.zM, ntt.mod)
	clear(vBuf[ntt.params.rank:])

	fwdNTTInPlacePow2(vBuf, ntt.ambNTT.tw[0], ntt.ambNTT.twM[0].SForm, ntt.mod.Value())

	vec.FMulTo(vBuf, vBuf, ntt.chirp, ntt.chirpM, ntt.mod)

	invNTTInPlacePow2(vBuf, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvM[0].SForm, ntt.mod.Value())

	vec.FMulTo(vNTT, vBuf[:ntt.params.rank], ntt.z, ntt.zM, ntt.mod)
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *anyCyclicTransformer) InverseTo(v, vNTT []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	vec.FMulTo(vBuf[:ntt.params.rank], vNTT, ntt.zInv, ntt.zInvM, ntt.mod)
	clear(vBuf[ntt.params.rank:])

	fwdNTTInPlacePow2(vBuf, ntt.ambNTT.tw[0], ntt.ambNTT.twM[0].SForm, ntt.mod.Value())

	vec.FMulTo(vBuf, vBuf, ntt.chirpInv, ntt.chirpInvM, ntt.mod)

	invNTTInPlacePow2(vBuf, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvM[0].SForm, ntt.mod.Value())

	vec.FMulTo(v, vBuf[:ntt.params.rank], ntt.zInv, ntt.zInvM, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *anyCyclicTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *anyCyclicTransformer) Modulus() *num.Modulus {
	return ntt.mod
}
