package dft

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

var (
	// cyclicNTTFactors are the factors of the rank for native cyclic NTT.
	cyclicNTTFactors = []uint64{2, 3, 5}
)

// cyclicTransformerBuffer is a buffer for [cyclicNativeTransformer].
type cyclicTransformerBuffer struct {
	// coeffs is the input.
	coeffs []uint64
}

func newCyclicTransformerBuffer(rank int) cyclicTransformerBuffer {
	return cyclicTransformerBuffer{
		coeffs: make([]uint64, rank),
	}
}

// cyclicPow235Transformer is a transformer for ranks multiple of [cyclicNTTFactors].
type cyclicPow235Transformer struct {
	params RingParameters
	mod    *num.Modulus

	// rankFactors are the factors of the rank.
	rankFactors []int

	// tw is the twiddle factor for NTT.
	// Ordered as [cyclicNTTFactors][Degree].
	tw [][]uint64
	// twS is the Shoup form of tw.
	// Ordered as [cyclicNTTFactors][Degree].
	twS [][]uint64
	// twInv is the twiddle factor for InvNTT.
	// Ordered as [cyclicNTTFactors][Degree].
	twInv [][]uint64
	// twInvS is the Shoup form of twInv.
	// Ordered as [cyclicNTTFactors][Degree].
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

	buf cyclicTransformerBuffer
}

// newCyclicPow235Transformer creates a new [cyclicNativeTransformer].
func newCyclicPow235Transformer(params RingParameters, mod *num.Modulus) *cyclicPow235Transformer {
	rankFactors := cyclicFactorDegree(params.rank, cyclicNTTFactors)

	root := num.PrimitiveRoot(mod)

	tw := make([][]uint64, len(cyclicNTTFactors))
	twS := make([][]uint64, len(cyclicNTTFactors))
	twInv := make([][]uint64, len(cyclicNTTFactors))
	twInvS := make([][]uint64, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		tw[i], twInv[i] = cyclicTwiddleFactor(int(rankFactors[i]), int(r), root, mod)
		twS[i] = vec.SForm(tw[i], mod)
		twInvS[i] = vec.SForm(twInv[i], mod)
	}

	rootPow := make([][]uint64, len(cyclicNTTFactors))
	rootPowS := make([][]uint64, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		rootPow[i] = make([]uint64, r)
		rootPow[i][0] = 1
		rootPow[i][1] = num.NthRoot(int(r), root, mod)
		for j := 2; j < int(r); j++ {
			rootPow[i][j] = num.Mul(rootPow[i][j-1], rootPow[i][1], mod)
		}
		rootPowS[i] = vec.SForm(rootPow[i], mod)
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

	return &cyclicPow235Transformer{
		params:      params,
		mod:         mod,
		rankFactors: rankFactors,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		root:  rootPow,
		rootS: rootPowS,

		rankInv: rankInv,

		idx: idx,

		buf: newCyclicTransformerBuffer(params.rank),
	}
}

func (ntt *cyclicPow235Transformer) ForwardInPlace(coeffs []uint64) {
	if len(ntt.idx) > 0 {
		for i := 0; i < ntt.params.rank; i++ {
			ntt.buf.coeffs[i] = coeffs[ntt.idx[i]]
		}
		copy(coeffs, ntt.buf.coeffs)
	}

	if ntt.rankFactors[0] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] {
			nttInPlacePow2(coeffs[i:i+ntt.rankFactors[0]], ntt.tw[0], ntt.twS[0], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[1] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] * ntt.rankFactors[1] {
			nttInPlacePow3(ntt.rankFactors[0], coeffs[i:i+ntt.rankFactors[0]*ntt.rankFactors[1]], ntt.tw[1], ntt.twS[1], ntt.root[1], ntt.rootS[1], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[2] > 1 {
		nttInPlacePow5(ntt.rankFactors[0]*ntt.rankFactors[1], coeffs, ntt.tw[2], ntt.twS[2], ntt.root[2], ntt.rootS[2], ntt.mod.Value())
	}

	vec.MFormTo(coeffs, coeffs, ntt.mod)
}

func (ntt *cyclicPow235Transformer) InverseInPlace(coeffs []uint64) {
	if ntt.rankFactors[0] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] {
			inttInPlacePow2(coeffs[i:i+ntt.rankFactors[0]], ntt.twInv[0], ntt.twInvS[0], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[1] > 1 {
		for i := 0; i < ntt.params.rank; i += ntt.rankFactors[0] * ntt.rankFactors[1] {
			inttInPlacePow3(ntt.rankFactors[0], coeffs[i:i+ntt.rankFactors[0]*ntt.rankFactors[1]], ntt.twInv[1], ntt.twInvS[1], ntt.root[1], ntt.rootS[1], ntt.mod.Value())
		}
	}

	if ntt.rankFactors[2] > 1 {
		inttInPlacePow5(ntt.rankFactors[0]*ntt.rankFactors[1], coeffs, ntt.twInv[2], ntt.twInvS[2], ntt.root[2], ntt.rootS[2], ntt.mod.Value())
	}

	if len(ntt.idx) > 0 {
		for i := 0; i < ntt.params.rank; i++ {
			ntt.buf.coeffs[ntt.idx[i]] = coeffs[i]
		}
		copy(coeffs, ntt.buf.coeffs)
	}

	vec.ScalarMulTo(coeffs, coeffs, ntt.rankInv, ntt.mod)
}

func (ntt *cyclicPow235Transformer) SafeCopy() Transformer {
	return &cyclicPow235Transformer{
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

		buf: newCyclicTransformerBuffer(ntt.params.rank),
	}
}

// cyclicBluesteinTransformer is a transformer for arbitrary ranks.
// Internally, it uses the Bluestein NTT.
type cyclicBluesteinTransformer struct {
	params RingParameters
	mod    *num.Modulus

	*cyclicPow235Transformer

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
}

// newCyclicBluesteinTransformer creates a new [cyclicBluesteinTransformer] for the given ringParams and modulus.
func newCyclicBluesteinTransformer(params RingParameters, mod *num.Modulus) *cyclicBluesteinTransformer {
	embedRank := int(num.NextProdPower(uint64(2*params.rank-1), []uint64{2}))
	embedParams := NewCyclicParameters(embedRank)

	root := num.PrimitiveRoot(mod)
	zz := num.NthRoot(2*params.rank, root, mod)
	zzInv := num.Inv(zz, mod)

	z := make([]uint64, params.rank)
	zInv := make([]uint64, params.rank)
	for i := 0; i < params.rank; i++ {
		idx := (i * i) % (2 * params.rank)
		z[i] = num.Exp(zz, uint64(idx), mod)
		zInv[i] = num.Exp(zzInv, uint64(idx), mod)
	}

	embedNTT := newCyclicPow235Transformer(embedParams, mod)

	chirpM := make([]uint64, embedRank)
	copy(chirpM, zInv)
	copy(chirpM[embedRank-params.rank+1:], zInv[1:])
	slices.Reverse(chirpM[embedRank-params.rank+1:])

	vec.ScalarMulTo(chirpM, chirpM, num.Inv(uint64(embedRank), mod), mod)
	nttInPlacePow2(chirpM, embedNTT.tw[0], embedNTT.twS[0], mod.Value())
	vec.MFormTo(chirpM, chirpM, mod)

	chirpInv := make([]uint64, embedRank)
	copy(chirpInv, z)
	copy(chirpInv[embedRank-params.rank+1:], z[1:])
	slices.Reverse(chirpInv[embedRank-params.rank+1:])

	vec.ScalarMulTo(chirpInv, chirpInv, num.Inv(uint64(embedRank*params.rank), mod), mod)
	nttInPlacePow2(chirpInv, embedNTT.tw[0], embedNTT.twS[0], mod.Value())
	vec.ReduceTo(chirpInv, chirpInv, mod)

	return &cyclicBluesteinTransformer{
		params: params,
		mod:    mod,

		cyclicPow235Transformer: embedNTT,

		z:     z,
		zS:    vec.SForm(z, mod),
		zInv:  zInv,
		zInvS: vec.SForm(zInv, mod),

		chirpM:   chirpM,
		chirpMS:  vec.SForm(chirpM, mod),
		chirpInv: chirpInv,
	}
}

func (ntt *cyclicBluesteinTransformer) ForwardInPlace(coeffs []uint64) {
	vec.SMulTo(ntt.buf.coeffs[:ntt.params.rank], coeffs, ntt.z, ntt.zS, ntt.mod)
	clear(ntt.buf.coeffs[ntt.params.rank:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.tw[0], ntt.twS[0], ntt.mod.Value())

	vec.SMulTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.chirpM, ntt.chirpMS, ntt.mod)

	inttInPlacePow2(ntt.buf.coeffs, ntt.twInv[0], ntt.twInvS[0], ntt.mod.Value())

	vec.SMulTo(coeffs, ntt.buf.coeffs[:ntt.params.rank], ntt.z, ntt.zS, ntt.mod)
}

func (ntt *cyclicBluesteinTransformer) InverseInPlace(coeffs []uint64) {
	vec.SMulTo(ntt.buf.coeffs[:ntt.params.rank], coeffs, ntt.zInv, ntt.zInvS, ntt.mod)
	clear(ntt.buf.coeffs[ntt.params.rank:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.tw[0], ntt.twS[0], ntt.mod.Value())

	vec.MMulTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.chirpInv, ntt.mod)

	inttInPlacePow2(ntt.buf.coeffs, ntt.twInv[0], ntt.twInvS[0], ntt.mod.Value())

	vec.SMulTo(coeffs, ntt.buf.coeffs[:ntt.params.rank], ntt.zInv, ntt.zInvS, ntt.mod)
}

func (ntt *cyclicBluesteinTransformer) SafeCopy() Transformer {
	embedNTTCopy := ntt.cyclicPow235Transformer.SafeCopy().(*cyclicPow235Transformer)
	return &cyclicBluesteinTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		cyclicPow235Transformer: embedNTTCopy,

		z:     ntt.z,
		zS:    ntt.zS,
		zInv:  ntt.zInv,
		zInvS: ntt.zInvS,

		chirpM:   ntt.chirpM,
		chirpMS:  ntt.chirpMS,
		chirpInv: ntt.chirpInv,
	}
}
