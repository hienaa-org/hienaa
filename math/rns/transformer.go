package rns

import (
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

// Transformer computes NTT/InvNTT transform.
type Transformer struct {
	modulus      []*mod.Modulus
	transformers []singleTransformer
}

// NewTransformer creates a new [Transformer] for the given modulus.
func NewTransformer(ringParams RingParameters, modulus []*mod.Modulus) *Transformer {
	transformers := make([]singleTransformer, len(modulus))
	for i, q := range modulus {
		transformers[i] = newSingleTransformer(ringParams, q)
	}

	return &Transformer{
		modulus:      modulus,
		transformers: transformers,
	}
}

// singleTransformer is an interface for NTT/InvNTT transforms for a single modulus.
type singleTransformer interface {
	// nttInPlace transforms the uint64 vector to NTT form.
	nttInPlace([]uint64)
	// invNTTInPlace transforms the uint64 vector to Standard form.
	invNTTInPlace([]uint64)
}

// newSingleTransformer creates a new [singleTransformer] for the given ring parameters and modulus.
func newSingleTransformer(ringParams RingParameters, modulus *mod.Modulus) singleTransformer {
	if isNTTFriendly(ringParams, modulus) {
		return newTrivialTransformer(ringParams, modulus)
	}

	switch ringParams.ringType {
	case Cyclotomic:
		switch {
		case num.IsPowerOfTwo(uint64(ringParams.cycloDegree)):
			return newCyclotomicPow2Transformer(ringParams, modulus)
		default:

		}
	case Cyclic:

	case AutFixed:
	}

	panic("newSingleTransformer: unsupported ring type or parameters")
}

// trivialTransformer is a no-op transformer that does nothing.
type trivialTransformer struct {
	params  RingParameters
	modulus *mod.Modulus
}

func newTrivialTransformer(ringParams RingParameters, modulus *mod.Modulus) *trivialTransformer {
	return &trivialTransformer{
		params:  ringParams,
		modulus: modulus,
	}
}

func (t *trivialTransformer) nttInPlace(coeffs []uint64) {
	mod.MFormVecTo(coeffs, t.modulus, coeffs)
}

func (t *trivialTransformer) invNTTInPlace(coeffs []uint64) {
	mod.InvMFormVecTo(coeffs, t.modulus, coeffs)
}
