package fhe

import (
	"github.com/hienaa-org/hienaa/fhe/bfv"
	"github.com/hienaa-org/hienaa/fhe/bgv"
	"github.com/hienaa-org/hienaa/fhe/ckks"
)

// Plaintext is a generic type for plaintext.
type Plaintext interface {
	*bfv.Plaintext | *bgv.Plaintext | *ckks.Plaintext
}

// Ciphertext is a generic type for ciphertext.
type Ciphertext interface {
	*bfv.Ciphertext | *bgv.Ciphertext | *ckks.Ciphertext
}
