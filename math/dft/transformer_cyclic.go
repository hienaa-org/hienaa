package dft

import (
	"slices"

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
	// twS is the Shoup form of tw.
	// Ordered as [cyclicNTTFactors][Rank].
	twS [][]uint64
	// twInv is the twiddle factor for InvNTT.
	// Ordered as [cyclicNTTFactors][Rank].
	twInv [][]uint64
	// twInvS is the Shoup form of twInv.
	// Ordered as [cyclicNTTFactors][Rank].
	twInvS [][]uint64

	// root is the powers of primitive roots.
	// Ordered as [cyclicNTTFactors][Radix].
	root [][]uint64
	// rootS is the Shoup form of root.
	// Ordered as [cyclicNTTFactors][Radix].
	rootS [][]uint64

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	// idx is the CRT mapping index.
	// Empty if the mapping is not needed, or in other words, rank is a prime power.
	idx []int

	buf transformerBuffer
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
	twS := make([][]uint64, len(cyclicNTTFactors))
	twInv := make([][]uint64, len(cyclicNTTFactors))
	twInvS := make([][]uint64, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		tw[i], twInv[i] = cyclicTwiddleFactor(rankFactors[i], r, root, mod)
		twS[i] = vec.SForm(tw[i], mod)
		twInvS[i] = vec.SForm(twInv[i], mod)
	}

	rootExp := make([][]uint64, len(cyclicNTTFactors))
	rootExpS := make([][]uint64, len(cyclicNTTFactors))
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
		rootExpS[i] = vec.SForm(rootExp[i], mod)
	}

	rankInv := num.InvMForm(num.Inv(uint64(params.rank), mod), mod)

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

	var buf transformerBuffer
	if len(idx) > 0 {
		buf = newTransformerBuffer(params.rank)
	}

	return &pow235CyclicTransformer{
		params:      params,
		mod:         mod,
		rankFactors: rankFactors,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		root:  rootExp,
		rootS: rootExpS,

		rankInv: rankInv,

		idx: idx,

		buf: buf,
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *pow235CyclicTransformer) ForwardTo(vNTT, v []uint64) {
	if len(ntt.idx) > 0 {
		for i := 0; i < ntt.params.rank; i++ {
			ntt.buf.coeffs[i] = v[ntt.idx[i]]
		}
		copy(vNTT, ntt.buf.coeffs)
	} else {
		copy(vNTT, v)
	}

	if ntt.rankFactors[0] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] {
			nttInPlacePow2(vNTT[i:i+ntt.rankFactors[0]], ntt.tw[0], ntt.twS[0], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[1] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] * ntt.rankFactors[1] {
			nttInPlacePow3(ntt.rankFactors[0], vNTT[i:i+ntt.rankFactors[0]*ntt.rankFactors[1]], ntt.tw[1], ntt.twS[1], ntt.root[1], ntt.rootS[1], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[2] > 1 {
		nttInPlacePow5(ntt.rankFactors[0]*ntt.rankFactors[1], vNTT, ntt.tw[2], ntt.twS[2], ntt.root[2], ntt.rootS[2], ntt.mod.Value())
	}

	vec.MFormTo(vNTT, vNTT, ntt.mod)
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *pow235CyclicTransformer) InverseTo(v, vNTT []uint64) {
	copy(v, vNTT)

	if ntt.rankFactors[0] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] {
			inttInPlacePow2(v[i:i+ntt.rankFactors[0]], ntt.twInv[0], ntt.twInvS[0], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[1] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] * ntt.rankFactors[1] {
			inttInPlacePow3(ntt.rankFactors[0], v[i:i+ntt.rankFactors[0]*ntt.rankFactors[1]], ntt.twInv[1], ntt.twInvS[1], ntt.root[1], ntt.rootS[1], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[2] > 1 {
		inttInPlacePow5(ntt.rankFactors[0]*ntt.rankFactors[1], v, ntt.twInv[2], ntt.twInvS[2], ntt.root[2], ntt.rootS[2], ntt.mod.Value())
	}

	if len(ntt.idx) > 0 {
		for i := 0; i < ntt.params.rank; i++ {
			ntt.buf.coeffs[ntt.idx[i]] = v[i]
		}
		copy(v, ntt.buf.coeffs)
	}

	vec.ScalarMulTo(v, v, ntt.rankInv, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *pow235CyclicTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *pow235CyclicTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

// SafeCopy returns a thread-safe copy.
func (ntt *pow235CyclicTransformer) SafeCopy() Transformer {
	if len(ntt.idx) == 0 {
		return ntt
	}

	return &pow235CyclicTransformer{
		params:      ntt.params,
		mod:         ntt.mod,
		rankFactors: ntt.rankFactors,

		tw:     ntt.tw,
		twS:    ntt.twS,
		twInv:  ntt.twInv,
		twInvS: ntt.twInvS,

		root:  ntt.root,
		rootS: ntt.rootS,

		rankInv: ntt.rankInv,

		idx: ntt.idx,

		buf: newTransformerBuffer(ntt.params.rank),
	}
}

// anyCyclicTransformer is a transformer for aribtrary rank cyclic ring.
// Internally, it uses Bluestein NTT.
type anyCyclicTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT *pow235CyclicTransformer

	// z is the factor for Z-transform.
	z []uint64
	// zS is the Shoup form of z.
	zS []uint64
	// zInv is the inverse factor for Z-transform.
	zInv []uint64
	// zInvS is the Shoup form of zInv.
	zInvS []uint64

	// chirpM is the chirp factor for Bluestein NTT in Montgomery form.
	chirpM []uint64
	// chirpMS is the Shoup form of chirpM.
	chirpMS []uint64
	// chirpInv is the inverse chirp factor for Bluestein NTT.
	chirpInv []uint64

	buf transformerBuffer
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

	chirpM := make([]uint64, ambRank)
	copy(chirpM, zInv)
	copy(chirpM[ambRank-params.rank+1:], zInv[1:])
	slices.Reverse(chirpM[ambRank-params.rank+1:])

	vec.ScalarMulTo(chirpM, chirpM, num.Inv(uint64(ambRank), mod), mod)
	nttInPlacePow2(chirpM, ambNTT.tw[0], ambNTT.twS[0], mod.Value())
	vec.MFormTo(chirpM, chirpM, mod)

	chirpInv := make([]uint64, ambRank)
	copy(chirpInv, z)
	copy(chirpInv[ambRank-params.rank+1:], z[1:])
	slices.Reverse(chirpInv[ambRank-params.rank+1:])

	vec.ScalarMulTo(chirpInv, chirpInv, num.Inv(uint64(ambRank*params.rank), mod), mod)
	nttInPlacePow2(chirpInv, ambNTT.tw[0], ambNTT.twS[0], mod.Value())
	vec.ReduceTo(chirpInv, chirpInv, mod)

	return &anyCyclicTransformer{
		params: params,
		mod:    mod,

		ambNTT: ambNTT,

		z:     z,
		zS:    vec.SForm(z, mod),
		zInv:  zInv,
		zInvS: vec.SForm(zInv, mod),

		chirpM:   chirpM,
		chirpMS:  vec.SForm(chirpM, mod),
		chirpInv: chirpInv,

		buf: newTransformerBuffer(ambRank),
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *anyCyclicTransformer) ForwardTo(vNTT, v []uint64) {
	vec.SMulLazyTo(ntt.buf.coeffs[:ntt.params.rank], v, ntt.z, ntt.zS, ntt.mod)
	clear(ntt.buf.coeffs[ntt.params.rank:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())

	vec.SMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.chirpM, ntt.chirpMS, ntt.mod)

	inttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())

	vec.SMulTo(vNTT, ntt.buf.coeffs[:ntt.params.rank], ntt.z, ntt.zS, ntt.mod)
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *anyCyclicTransformer) InverseTo(v, vNTT []uint64) {
	vec.SMulLazyTo(ntt.buf.coeffs[:ntt.params.rank], vNTT, ntt.zInv, ntt.zInvS, ntt.mod)
	clear(ntt.buf.coeffs[ntt.params.rank:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())

	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.chirpInv, ntt.mod)

	inttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())

	vec.SMulTo(v, ntt.buf.coeffs[:ntt.params.rank], ntt.zInv, ntt.zInvS, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *anyCyclicTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *anyCyclicTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

// SafeCopy returns a thread-safe copy.
func (ntt *anyCyclicTransformer) SafeCopy() Transformer {
	return &anyCyclicTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		ambNTT: ntt.ambNTT.SafeCopy().(*pow235CyclicTransformer),

		z:     ntt.z,
		zS:    ntt.zS,
		zInv:  ntt.zInv,
		zInvS: ntt.zInvS,

		chirpM:   ntt.chirpM,
		chirpMS:  ntt.chirpMS,
		chirpInv: ntt.chirpInv,

		buf: newTransformerBuffer(ntt.ambNTT.params.rank),
	}
}
