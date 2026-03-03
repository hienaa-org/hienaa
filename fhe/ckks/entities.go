package ckks

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
)

// SecretKey is a CKKS secret key.
type SecretKey struct {
	Value *crt.Element
}

// Plaintext is a CKKS plaintext embedded in ciphertext modulus.
type Plaintext struct {
	Value *crt.Element
}

// Ciphertext is a CKKS ciphertext.
type Ciphertext struct {
	Value *rlwe.Ciphertext
	scFac *big.Float
}
