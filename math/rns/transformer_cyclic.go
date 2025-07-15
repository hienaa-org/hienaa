package rns

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

var (
	// cyclicNTTFactors are the factors of the degree for native cyclic NTT.
	cyclicNTTFactors = []uint64{2, 3, 5}
)

// cyclicTransformerBuffer is a buffer for [cyclicNativeTransformer].
type cyclicTransformerBuffer struct {
	// coeffs is the input.
	coeffs []uint64
}

func newCyclicTransformerBuffer(degree int) cyclicTransformerBuffer {
	return cyclicTransformerBuffer{
		coeffs: make([]uint64, degree),
	}
}

// cyclicPow235Transformer is a transformer for degrees multiple of [cyclicNTTFactors].
type cyclicPow235Transformer struct {
	params  RingParameters
	modulus *mod.Modulus

	// degFactors are the factors of the degree.
	degFactors []int

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

	// degInv is the modular inverse of the degree.
	degInv uint64

	// idx is the CRT mapping index.
	// Empty if the mapping is not needed, or in other words, degree is a prime power.
	idx []int

	buf cyclicTransformerBuffer
}

// newCyclicPow235Transformer creates a new [cyclicNativeTransformer] for the given ringParams and modulus.
func newCyclicPow235Transformer(params RingParameters, modulus *mod.Modulus) *cyclicPow235Transformer {
	degFactors := cyclicFactorDegree(params.degree, cyclicNTTFactors)

	root := num.PrimitiveRoot(modulus)

	tw := make([][]uint64, len(cyclicNTTFactors))
	twS := make([][]uint64, len(cyclicNTTFactors))
	twInv := make([][]uint64, len(cyclicNTTFactors))
	twInvS := make([][]uint64, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		tw[i], twInv[i] = cyclicTwiddleFactor(int(degFactors[i]), int(r), root, modulus)
		twS[i] = mod.SFormVec(tw[i], modulus)
		twInvS[i] = mod.SFormVec(twInv[i], modulus)
	}

	rootPow := make([][]uint64, len(cyclicNTTFactors))
	rootPowS := make([][]uint64, len(cyclicNTTFactors))
	for i, r := range cyclicNTTFactors {
		rootPow[i] = make([]uint64, r)
		rootPow[i][0] = 1
		rootPow[i][1] = num.NthRoot(int(r), root, modulus)
		for j := 2; j < int(r); j++ {
			rootPow[i][j] = mod.Mul(rootPow[i][j-1], rootPow[i][1], modulus)
		}
		rootPowS[i] = mod.SFormVec(rootPow[i], modulus)
	}

	degInv := mod.InvMForm(mod.Inv(uint64(params.degree), modulus), modulus)

	var idx []int
	if !slices.Contains(degFactors, params.degree) {
		idx = make([]int, params.degree)
		for i2 := 0; i2 < degFactors[0]; i2++ {
			for i3 := 0; i3 < degFactors[1]; i3++ {
				for i5 := 0; i5 < degFactors[2]; i5++ {
					idxIn := i2 + i3*degFactors[0] + i5*degFactors[0]*degFactors[1]
					idxOut := i2*(params.degree/degFactors[0]) + i3*(params.degree/degFactors[1]) + i5*(params.degree/degFactors[2])
					idx[idxIn%params.degree] = idxOut % params.degree
				}
			}
		}
	}

	return &cyclicPow235Transformer{
		params:     params,
		modulus:    modulus,
		degFactors: degFactors,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		root:  rootPow,
		rootS: rootPowS,

		degInv: degInv,

		idx: idx,

		buf: newCyclicTransformerBuffer(params.degree),
	}
}

func (ntt *cyclicPow235Transformer) nttInPlace(coeffs []uint64) {
	if len(ntt.idx) > 0 {
		for i := 0; i < ntt.params.degree; i++ {
			ntt.buf.coeffs[i] = coeffs[ntt.idx[i]]
		}
		copy(coeffs, ntt.buf.coeffs)
	}

	if ntt.degFactors[0] > 1 {
		for i := 0; i < ntt.params.degree; i += ntt.degFactors[0] {
			nttInPlacePow2(coeffs[i:i+ntt.degFactors[0]], ntt.tw[0], ntt.twS[0], ntt.modulus.Value())
		}
	}

	if ntt.degFactors[1] > 1 {
		for i := 0; i < ntt.params.degree; i += ntt.degFactors[0] * ntt.degFactors[1] {
			nttInPlacePow3(ntt.degFactors[0], coeffs[i:i+ntt.degFactors[0]*ntt.degFactors[1]], ntt.tw[1], ntt.twS[1], ntt.root[1], ntt.rootS[1], ntt.modulus.Value())
		}
	}

	if ntt.degFactors[2] > 1 {
		nttInPlacePow5(ntt.degFactors[0]*ntt.degFactors[1], coeffs, ntt.tw[2], ntt.twS[2], ntt.root[2], ntt.rootS[2], ntt.modulus.Value())
	}

	mod.MFormVecTo(coeffs, ntt.modulus, coeffs)
}

func (ntt *cyclicPow235Transformer) invNTTInPlace(coeffs []uint64) {
	if ntt.degFactors[0] > 1 {
		for i := 0; i < ntt.params.degree; i += ntt.degFactors[0] {
			inttInPlacePow2(coeffs[i:i+ntt.degFactors[0]], ntt.twInv[0], ntt.twInvS[0], ntt.modulus.Value())
		}
	}

	if ntt.degFactors[1] > 1 {
		for i := 0; i < ntt.params.degree; i += ntt.degFactors[0] * ntt.degFactors[1] {
			inttInPlacePow3(ntt.degFactors[0], coeffs[i:i+ntt.degFactors[0]*ntt.degFactors[1]], ntt.twInv[1], ntt.twInvS[1], ntt.root[1], ntt.rootS[1], ntt.modulus.Value())
		}
	}

	if ntt.degFactors[2] > 1 {
		inttInPlacePow5(ntt.degFactors[0]*ntt.degFactors[1], coeffs, ntt.twInv[2], ntt.twInvS[2], ntt.root[2], ntt.rootS[2], ntt.modulus.Value())
	}

	if len(ntt.idx) > 0 {
		for i := 0; i < ntt.params.degree; i++ {
			ntt.buf.coeffs[ntt.idx[i]] = coeffs[i]
		}
		copy(coeffs, ntt.buf.coeffs)
	}

	mod.ScalarMulVecTo(coeffs, ntt.degInv, ntt.modulus, coeffs)
}

func (ntt *cyclicPow235Transformer) shallowCopy() singleTransformer {
	return &cyclicPow235Transformer{
		params:     ntt.params,
		modulus:    ntt.modulus,
		degFactors: ntt.degFactors,

		tw:     ntt.tw,
		twS:    ntt.twS,
		twInv:  ntt.twInv,
		twInvS: ntt.twInvS,

		root:  ntt.root,
		rootS: ntt.rootS,

		degInv: ntt.degInv,

		idx: ntt.idx,

		buf: newCyclicTransformerBuffer(ntt.params.degree),
	}
}

// cyclicBluesteinTransformer is a transformer for arbitrary degrees.
// Internally, it uses the Bluestein NTT.
type cyclicBluesteinTransformer struct {
	params  RingParameters
	modulus *mod.Modulus

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
func newCyclicBluesteinTransformer(params RingParameters, modulus *mod.Modulus) *cyclicBluesteinTransformer {
	embedDeg := int(num.NextProdPower(uint64(2*params.degree-1), []uint64{2}))
	embedParams := NewCyclicParameters(embedDeg)

	root := num.PrimitiveRoot(modulus)
	zz := num.NthRoot(2*params.degree, root, modulus)
	zzInv := mod.Inv(zz, modulus)

	z := make([]uint64, params.degree)
	zInv := make([]uint64, params.degree)
	for i := 0; i < params.degree; i++ {
		idx := (i * i) % (2 * params.degree)
		z[i] = mod.Exp(zz, uint64(idx), modulus)
		zInv[i] = mod.Exp(zzInv, uint64(idx), modulus)
	}

	embedNTT := newCyclicPow235Transformer(embedParams, modulus)

	chirpM := make([]uint64, embedDeg)
	copy(chirpM, zInv)
	copy(chirpM[embedDeg-params.degree+1:], zInv[1:])
	slices.Reverse(chirpM[embedDeg-params.degree+1:])

	mod.ScalarMulVecTo(chirpM, mod.Inv(uint64(embedDeg), modulus), modulus, chirpM)
	nttInPlacePow2(chirpM, embedNTT.tw[0], embedNTT.twS[0], modulus.Value())
	mod.MFormVecTo(chirpM, modulus, chirpM)

	chirpInv := make([]uint64, embedDeg)
	copy(chirpInv, z)
	copy(chirpInv[embedDeg-params.degree+1:], z[1:])
	slices.Reverse(chirpInv[embedDeg-params.degree+1:])

	mod.ScalarMulVecTo(chirpInv, mod.Inv(uint64(embedDeg*params.degree), modulus), modulus, chirpInv)
	nttInPlacePow2(chirpInv, embedNTT.tw[0], embedNTT.twS[0], modulus.Value())
	mod.ReduceVecTo(chirpInv, modulus, chirpInv)

	return &cyclicBluesteinTransformer{
		params:  params,
		modulus: modulus,

		cyclicPow235Transformer: embedNTT,

		z:     z,
		zS:    mod.SFormVec(z, modulus),
		zInv:  zInv,
		zInvS: mod.SFormVec(zInv, modulus),

		chirpM:   chirpM,
		chirpMS:  mod.SFormVec(chirpM, modulus),
		chirpInv: chirpInv,
	}
}

func (ntt *cyclicBluesteinTransformer) nttInPlace(coeffs []uint64) {
	mod.SMulVecTo(coeffs, ntt.z, ntt.zS, ntt.modulus, ntt.buf.coeffs[:ntt.params.degree])
	clear(ntt.buf.coeffs[ntt.params.degree:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.tw[0], ntt.twS[0], ntt.modulus.Value())

	mod.SMulVecTo(ntt.buf.coeffs, ntt.chirpM, ntt.chirpMS, ntt.modulus, ntt.buf.coeffs)

	inttInPlacePow2(ntt.buf.coeffs, ntt.twInv[0], ntt.twInvS[0], ntt.modulus.Value())

	mod.SMulVecTo(ntt.buf.coeffs[:ntt.params.degree], ntt.z, ntt.zS, ntt.modulus, coeffs)
}

func (ntt *cyclicBluesteinTransformer) invNTTInPlace(coeffs []uint64) {
	mod.SMulVecTo(coeffs, ntt.zInv, ntt.zInvS, ntt.modulus, ntt.buf.coeffs[:ntt.params.degree])
	clear(ntt.buf.coeffs[ntt.params.degree:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.tw[0], ntt.twS[0], ntt.modulus.Value())

	mod.MMulVecTo(ntt.buf.coeffs, ntt.chirpInv, ntt.modulus, ntt.buf.coeffs)

	inttInPlacePow2(ntt.buf.coeffs, ntt.twInv[0], ntt.twInvS[0], ntt.modulus.Value())

	mod.SMulVecTo(ntt.buf.coeffs[:ntt.params.degree], ntt.zInv, ntt.zInvS, ntt.modulus, coeffs)
}

func (ntt *cyclicBluesteinTransformer) shallowCopy() singleTransformer {
	embedNTTCopy := ntt.cyclicPow235Transformer.shallowCopy().(*cyclicPow235Transformer)
	return &cyclicBluesteinTransformer{
		params:  ntt.params,
		modulus: ntt.modulus,

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
