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
