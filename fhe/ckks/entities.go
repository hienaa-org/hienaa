package ckks

import (
	"slices"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
)

// Plaintext is a CKKS plaintext.
type Plaintext []float64

// NewScalar creates a new [Plaintext] for a scalar.
func NewScalar() Plaintext {
	return make([]float64, 1)
}

// NewScalarFrom creates a new [Plaintext] for a scalar from x.
func NewScalarFrom(x float64) Plaintext {
	res := NewScalar()
	res[0] = x
	return res
}

// NewPoly creates a new [Plaintext] for a polynomial.
func NewPoly(rank int) Plaintext {
	return make([]float64, rank)
}

// NewPolyFrom creates a new [Plaintext] for a polynomial from x.
func NewPolyFrom(x []float64) Plaintext {
	res := NewPoly(len(x))
	copy(res, x)
	return res
}

// Rank returns the rank.
func (e Plaintext) Rank() int {
	return len(e)
}

// Clear clears value.
func (e Plaintext) Clear() {
	clear(e)
}

// Copy returns a copy of e.
func (e Plaintext) Copy() Plaintext {
	res := make([]float64, len(e))
	copy(res, e)
	return res
}

// CopyFrom copies the coefficients from eIn to e.
func (e Plaintext) CopyFrom(eIn Plaintext) {
	copy(e, eIn)
}

// IsEqual checks if two values are equal.
func (e Plaintext) IsEqual(e0 Plaintext) bool {
	return slices.Equal(e, e0)
}

// Ciphertext is a CKKS ciphertext.
type Ciphertext struct {
	Value *rlwe.Ciphertext
	scFac float64
}

// NewCiphertext creates a new [Ciphertext].
func NewCiphertext(params rlwe.Parameters, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Value: rlwe.NewCiphertext(params, false, isNTT),
		scFac: 0,
	}
}

// NewCiphertextCustom creates a new [Ciphertext] with custom parameters.
func NewCiphertextCustom(rank, modLen int, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Value: rlwe.NewCiphertextCustom(rank, modLen, 0, isNTT),
		scFac: 0,
	}
}

// ScFac returns the scaling factor.
func (c *Ciphertext) ScFac() float64 {
	return c.scFac
}

// Rank returns the rank.
func (c *Ciphertext) Rank() int {
	return c.Value.Rank()
}

// BaseModLen returns the base modulus length.
func (c *Ciphertext) ModLen() int {
	return c.Value.BaseModLen()
}

// IsNTT returns whether value is in NTT form.
func (c *Ciphertext) IsNTT() bool {
	return c.Value.IsNTT()
}

// Clear clears value.
func (c *Ciphertext) Clear() {
	c.Value.Clear()
	c.scFac = 0
}

// WithModLen returns a shallow copy with the given modulus lengths.
func (c *Ciphertext) WithModLen(baseLen int) *Ciphertext {
	return &Ciphertext{
		Value: c.Value.WithModLen(baseLen, 0),
		scFac: c.scFac,
	}
}

// Copy returns a copy of c.
func (c *Ciphertext) Copy() *Ciphertext {
	return &Ciphertext{
		Value: c.Value.Copy(),
		scFac: c.scFac,
	}
}

// CopyFrom copies the coefficients from cIn to c.
func (c *Ciphertext) CopyFrom(cIn *Ciphertext) {
	c.Value.CopyFrom(cIn.Value)
	c.scFac = cIn.scFac
}

// IsEqual checks if two values are equal.
func (c *Ciphertext) IsEqual(c0 *Ciphertext) bool {
	return c.Value.IsEqual(c0.Value) && c.scFac == c0.scFac
}

// IsConsistent checks if two values have the same shape.
func (c *Ciphertext) IsConsistent(c0 *Ciphertext) bool {
	return c.Value.IsConsistent(c0.Value) && c.scFac == c0.scFac
}
