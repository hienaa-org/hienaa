package ckks

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
)

// SecretKey is a CKKS secret key.
type SecretKey struct {
	Value *crt.Poly
}

// Plaintext is a CKKS plaintext embedded in ciphertext modulus.
type Plaintext struct {
	Value *crt.Poly
}

// Ciphertext is a CKKS ciphertext.
type Ciphertext struct {
	Value *rlwe.Ciphertext
	noise *big.Int
}
