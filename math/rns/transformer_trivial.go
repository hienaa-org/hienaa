package rns

import "github.com/hienaa-org/hienaa/math/mod"

// trivialTransformer is a no-op transformer that does nothing.
type trivialTransformer struct {
	modulus *mod.Modulus
}

func newTrivialTransformer(modulus *mod.Modulus) *trivialTransformer {
	return &trivialTransformer{
		modulus: modulus,
	}
}

func (t *trivialTransformer) nttInPlace(coeffs []uint64) {
	if t.modulus.Inv() != 0 {
		mod.MFormVecTo(coeffs, coeffs, t.modulus)
	}
}

func (t *trivialTransformer) invNTTInPlace(coeffs []uint64) {
	if t.modulus.Inv() != 0 {
		mod.InvMFormVecTo(coeffs, coeffs, t.modulus)
	}
}

func (t *trivialTransformer) safeCopy() singleTransformer {
	return t
}
