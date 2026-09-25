package dft

import "github.com/hienaa-org/hienaa/math/num"

// dftType is a type for DFT algorithm.
type dftType uint64

const (
	// typePow2Cyclotomic is a power-of-two cyclotomic ring.
	typePow2Cyclotomic = iota
	// typeAnyCyclotomic is an arbitrary index cyclotomic ring.
	typeAnyCyclotomic
	// typePow235Cyclic is a cyclic ring with rank multiples of 2, 3, 5.
	typePow235Cyclic
	// typeAnyCyclic is an arbitrary rank cyclic ring.
	typeAnyCyclic
	// typePow2AutFixed is a power-of-two conjugate invariant ring.
	typePow2AutFixed
	// typePrimeAutFixed is a prime order fixed automorphism subring.
	typePrimeAutFixed
	// typeNone is a invalid type.
	typeNone
)

// dftTypeOf returns the [dftType] of [RingParameters].
func dftTypeOf(params RingParameters) dftType {
	switch params.ringType {
	case TypeCyclotomic:
		if num.IsPowerOfTwo(params.cycloIdx) {
			return typePow2Cyclotomic
		}
		return typeAnyCyclotomic

	case TypeCyclic:
		if num.IsProdPowerOf(params.rank, cyclicNTTFactors) {
			return typePow235Cyclic
		}
		return typeAnyCyclic

	case TypeAutFixed:
		if num.IsPowerOfTwo(params.cycloIdx) {
			return typePow2AutFixed
		}
		return typePrimeAutFixed
	}

	return typeNone
}

// Transformer computes NTT/InvNTT transform.
//
// After the transform, the coefficients are in Montgomery form for efficient convolution.
type Transformer struct {
	params  RingParameters
	dftType dftType
	mod     *num.Modulus

	*pow2CyclotomicTransformer
	*anyCyclotomicTransformer
	*pow235CyclicTransformer
	*anyCyclicTransformer
	*pow2AutFixedTransformer
	*primeAutFixedTransformer
}

// NewTransformer creates a new [Transformer].
//
// Panics when the ring parameters or modulus are unsupported.
func NewTransformer(params RingParameters, mod *num.Modulus) *Transformer {
	if !IsNTTFriendly(params, mod) {
		panic("unsupported ring type or parameters")
	}

	ntt := &Transformer{
		params:  params,
		dftType: dftTypeOf(params),
		mod:     mod,
	}

	switch ntt.dftType {
	case typePow2Cyclotomic:
		ntt.pow2CyclotomicTransformer = newPow2CyclotomicTransformer(params, mod)
	case typeAnyCyclotomic:
		ntt.anyCyclotomicTransformer = newAnyCyclotomicTransformer(params, mod)
	case typePow235Cyclic:
		ntt.pow235CyclicTransformer = newPow235CyclicTransformer(params, mod)
	case typeAnyCyclic:
		ntt.anyCyclicTransformer = newAnyCyclicTransformer(params, mod)
	case typePow2AutFixed:
		ntt.pow2AutFixedTransformer = newPow2AutFixedTransformer(params, mod)
	case typePrimeAutFixed:
		ntt.primeAutFixedTransformer = newPrimeAutFixedTransformer(params, mod)
	}

	return ntt
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *Transformer) ForwardTo(vNTT, v []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	switch ntt.dftType {
	case typePow2Cyclotomic:
		ntt.pow2CyclotomicTransformer.forwardTo(vNTT, v)
	case typeAnyCyclotomic:
		ntt.anyCyclotomicTransformer.forwardTo(vNTT, v)
	case typePow235Cyclic:
		ntt.pow235CyclicTransformer.forwardTo(vNTT, v)
	case typeAnyCyclic:
		ntt.anyCyclicTransformer.forwardTo(vNTT, v)
	case typePow2AutFixed:
		ntt.pow2AutFixedTransformer.forwardTo(vNTT, v)
	case typePrimeAutFixed:
		ntt.primeAutFixedTransformer.forwardTo(vNTT, v)
	}
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *Transformer) InverseTo(v, vNTT []uint64) {
	checkLength(ntt.params.rank, len(v), len(vNTT))

	switch ntt.dftType {
	case typePow2Cyclotomic:
		ntt.pow2CyclotomicTransformer.inverseTo(v, vNTT)
	case typeAnyCyclotomic:
		ntt.anyCyclotomicTransformer.inverseTo(v, vNTT)
	case typePow235Cyclic:
		ntt.pow235CyclicTransformer.inverseTo(v, vNTT)
	case typeAnyCyclic:
		ntt.anyCyclicTransformer.inverseTo(v, vNTT)
	case typePow2AutFixed:
		ntt.pow2AutFixedTransformer.inverseTo(v, vNTT)
	case typePrimeAutFixed:
		ntt.primeAutFixedTransformer.inverseTo(v, vNTT)
	}
}

// Params returns the [RingParameters] of transformer.
func (ntt *Transformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the [*num.Modulus] of transformer.
func (ntt *Transformer) Modulus() *num.Modulus {
	return ntt.mod
}
