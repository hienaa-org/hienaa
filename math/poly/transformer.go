package poly

import (
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

// Transformer is an interface for NTT/InvNTT transforms.
type Transformer interface {
	// RingParameters returns the [RingParameters] used by the transformer.
	RingParameters() RingParameters
	// NTTInPlace transforms the uint64 vector to NTT form.
	NTTInPlace([]uint64)
	// InvNTTInPlace transforms the uint64 vector to Standard form.
	InvNTTInPlace([]uint64)
	// Type is the type of the transformer.
	Type() TransformType
}

// TransformType is an algorithm type for NTT/InvNTT.
type TransformType uint64

const (
	// PowerOfTwo is an TransformType for N = 2^a.
	// This is the most common and fastest case.
	PowerOfTwo TransformType = iota
	// ProductOfPrimes is an TransformType for N = 2^a * 3^b * 5^c.
	// Uses Good-Thomas NTT.
	ProductOfPrimes
	// General is an TransformType for general N.
	// Uses Bluestein NTT.
	General
)

// NewTransformer creates a new [Transformer] for the given ringParams and modulus.
func NewTransformer(ringParams RingParameters, modulus *mod.Modulus) Transformer {
	switch ringParams.ringType {
	case Cyclotomic:
		switch {
		case num.IsPowerOfTwo(uint64(ringParams.cycloDegree)):
			return newCyclotomicPow2Transformer(ringParams, modulus)
		}
	}

	panic("NewTransformer: unsupported ring parameters")
}

// TrivialTransformer is a no-op transformer that does nothing.
type TrivialTransformer struct {
	Params        RingParameters
	TransformType TransformType
}

func (t *TrivialTransformer) RingParameters() RingParameters {
	return t.Params
}

func (t *TrivialTransformer) NTTInPlace(coeffs []uint64) {}

func (t *TrivialTransformer) InvNTTInPlace(coeffs []uint64) {}

func (t *TrivialTransformer) Type() TransformType {
	return t.TransformType
}
