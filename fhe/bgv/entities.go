package bgv

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
)

// SecretKey is a BGV secret key.
type SecretKey struct {
	Value *crt.Poly
}

// Plaintext is a BGV plaintext embedded in ciphertext modulus.
type Plaintext struct {
	Value *crt.Poly
}

// Ciphertext is a BGV ciphertext.
type Ciphertext struct {
	Value *rlwe.Ciphertext
	noise *big.Int
}
