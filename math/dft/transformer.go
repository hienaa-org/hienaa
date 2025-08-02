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
func NewTransformer(ringParams RingParameters, mod *num.Modulus) Transformer {
	switch ringParams.ringType {
	case Cyclotomic:
		switch {
		case num.IsPowerOfTwo(uint64(ringParams.cycloOrd)):
			return newCyclotomicPow2Transformer(ringParams, mod)
		default:
			return newCyclotomicAnyTransformer(ringParams, mod)
		}
	case Cyclic:
		switch {
		case num.IsProdPowerOf(uint64(ringParams.rank), cyclicNTTFactors):
			return newCyclicPow235Transformer(ringParams, mod)
		default:
			return newCyclicBluesteinTransformer(ringParams, mod)
		}
	case AutFixed:
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
