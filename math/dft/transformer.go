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
	// ForwardTo transforms the uint64 vector to NTT form.
	ForwardTo(vNTT, v []uint64)
	// InverseTo transforms the uint64 vector to Standard form.
	InverseTo(v, vNTT []uint64)
	// SafeCopy returns a thread-safe copy.
	SafeCopy() Transformer
}

// NewTransformer creates a new [Transformer].
//
// Panics when the ring parameters or modulus are unsupported.
func NewTransformer(params RingParameters, mod *num.Modulus) Transformer {
	if !IsNTTFriendly(params, mod) {
		panic("unsupported ring type or parameters")
	}

	switch params.ringType {
	case Cyclotomic:
		switch {
		case num.IsPowerOfTwo(params.cycloOrd):
			return newPow2CyclotomicTransformer(params, mod)
		default:
			return newAnyCyclotomicTransformer(params, mod)
		}
	case Cyclic:
		switch {
		case num.IsProdPowerOf(params.rank, cyclicNTTFactors):
			return newCyclicPow235Transformer(params, mod)
		default:
			return newAnyCyclicTransformer(params, mod)
		}
	case AutFixed:
		switch {
		case num.IsPowerOfTwo(params.cycloOrd):
			return newPow2AutFixedTransformer(params, mod)
		case num.IsPrime(params.cycloOrd):
			return newPrimeAutFixedTransformer(params, mod)
		}
	}

	panic("unsupported ring type or parameters")
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
