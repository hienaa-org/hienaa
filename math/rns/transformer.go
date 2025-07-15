package rns

import (
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

// Transformer computes NTT/InvNTT transform.
type Transformer struct {
	Params       RingParameters
	Modulus      []*mod.Modulus
	transformers []singleTransformer
}

// NewTransformer creates a new [Transformer] for the given modulus.
func NewTransformer(ringParams RingParameters, modulus []*mod.Modulus) *Transformer {
	transformers := make([]singleTransformer, len(modulus))
	for i, q := range modulus {
		transformers[i] = newSingleTransformer(ringParams, q)
	}

	return &Transformer{
		Params:       ringParams,
		Modulus:      modulus,
		transformers: transformers,
	}
}

// NTT returns NTT(p).
func (ntt *Transformer) NTT(p *Poly) *Poly {
	pOut := NewNTTPoly(ntt.Params.degree, len(ntt.Modulus))
	ntt.NTTTo(p, pOut)
	return pOut
}

// NTTTo computes pOut = NTT(p).
func (ntt *Transformer) NTTTo(p, pOut *Poly) {
	if p.IsNTT {
		panic("NTTTo: input polynomial is in NTT form")
	}

	for i := range ntt.transformers {
		copy(pOut.Coeffs[i], p.Coeffs[i])
		ntt.transformers[i].nttInPlace(pOut.Coeffs[i])
	}
	pOut.IsNTT = true
}

// InvNTT returns InvNTT(p).
func (ntt *Transformer) InvNTT(p *Poly) *Poly {
	pOut := NewPoly(ntt.Params.degree, len(ntt.Modulus))
	ntt.InvNTTTo(p, pOut)
	return pOut
}

// InvNTTTo computes pOut = InvNTT(p).
func (ntt *Transformer) InvNTTTo(p, pOut *Poly) {
	if !p.IsNTT {
		panic("InvNTTTo: input polynomial is in Standard form")
	}

	for i := range ntt.transformers {
		copy(pOut.Coeffs[i], p.Coeffs[i])
		ntt.transformers[i].invNTTInPlace(pOut.Coeffs[i])
	}
	pOut.IsNTT = false
}

// singleTransformer is an interface for NTT/InvNTT transforms for a single modulus.
type singleTransformer interface {
	// nttInPlace transforms the uint64 vector to NTT form.
	nttInPlace([]uint64)
	// invNTTInPlace transforms the uint64 vector to Standard form.
	invNTTInPlace([]uint64)
	// shallowCopy creates a thread-safe copy of the singleTransformer.
	shallowCopy() singleTransformer
}

// newSingleTransformer creates a new [singleTransformer] for the given ring parameters and modulus.
func newSingleTransformer(ringParams RingParameters, modulus *mod.Modulus) singleTransformer {
	if !isNTTFriendly(ringParams, modulus) {
		return newTrivialTransformer(modulus)
	}

	switch ringParams.ringType {
	case Cyclotomic:
		switch {
		case num.IsPowerOfTwo(uint64(ringParams.cycloDegree)):
			return newCyclotomicPow2Transformer(ringParams, modulus)
		}
	case Cyclic:
		switch {
		case num.IsProdPowerOf(uint64(ringParams.degree), cyclicNTTFactors):
			return newCyclicPow235Transformer(ringParams, modulus)
		default:
			return newCyclicBluesteinTransformer(ringParams, modulus)
		}
	case AutFixed:
	}

	panic("newSingleTransformer: unsupported ring type or parameters")
}
