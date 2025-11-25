package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
)

// SecretKey is a RLWE secret key.
type SecretKey struct {
	Value *crt.Poly
}

// Ciphertext is a RLWE ciphertext.
type Ciphertext struct {
	Body *crt.Poly
	Mask *crt.Poly
	IsPQ bool
}

// Tensor is a vector of polynomials.
type Tensor struct {
	Value []*crt.Poly
	IsPQ  bool
}

// GadgetEncryption is a gadget encryption.
type GadgetEncryption struct {
	Value  []*Ciphertext
	params GadgetParameters
}
