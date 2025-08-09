package dft

import (
	"github.com/hienaa-org/hienaa/math/num"
)

// Transformer is an interface for NTT/InvNTT transforms.
//
// After the transform, the coefficients are in Montgomery form for efficient convolution.
type Transformer interface {
	// Params returns the ring parameters.
	Params() RingParameters
	// Modulus returns the modulus used for the transform.
	Modulus() *num.Modulus
	// ForwardInPlace transforms the uint64 vector to NTT form.
	ForwardInPlace(coeffs []uint64)
	// InverseInPlace transforms the uint64 vector to Standard form.
	InverseInPlace(coeffs []uint64)
	// SafeCopy returns a thread-safe copy.
	SafeCopy() Transformer
}

// NewTransformer creates a new [Transformer].
func NewTransformer(params RingParameters, mod *num.Modulus) Transformer {
	if !IsNTTFriendly(params, mod) {
		panic("NewTransformer: unsupported ring parameters or modulus")
	}

	switch params.ringType {
	case Cyclotomic:
		switch {
		case num.IsPowerOfTwo(uint64(params.cycloOrd)):
			return newCyclotomicPow2Transformer(params, mod)
		default:
			return newCyclotomicAnyTransformer(params, mod)
		}
	case Cyclic:
		switch {
		case num.IsProdPowerOf(uint64(params.rank), cyclicNTTFactors):
			return newCyclicPow235Transformer(params, mod)
		default:
			return newCyclicBluesteinTransformer(params, mod)
		}
	case AutFixed:
		switch {
		case num.IsPrime(uint64(params.cycloOrd)):
			return newAutFixedPrimeTransformer(params, mod)
		case num.IsPowerOfTwo(uint64(params.cycloOrd)):
			return newAutFixedPow2Transformer(params, mod)
		}
	}

	panic("NewTransformer: unsupported ring type or parameters")
}

// transformerBuffer is a buffer for [Transformer].
type transformerBuffer struct {
	// coeffs is the input.
	coeffs []uint64
}

// newTransformerBuffer creates a new [transformerBuffer].
func newTransformerBuffer(rank int) transformerBuffer {
	return transformerBuffer{
		coeffs: make([]uint64, rank),
	}
}
