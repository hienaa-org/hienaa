package bfv

import (
	"slices"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Plaintext is a BFV plaintext embedded in ciphertext modulus.
type Plaintext []uint64

// NewScalar creates a new [Plaintext] for a scalar.
func NewScalar() Plaintext {
	return make([]uint64, 1)
}

// NewScalarFrom creates a new [Plaintext] for a scalar from x.
func NewScalarFrom(x uint64, msgMod *num.Modulus) Plaintext {
	res := NewScalar()
	res[0] = num.Reduce(x, msgMod)
	return res
}

// NewPoly creates a new [Plaintext] for a polynomial.
func NewPoly(rank int) Plaintext {
	return make([]uint64, rank)
}

// NewPolyFrom creates a new [Plaintext] for a polynomial from x.
func NewPolyFrom(x []uint64, msgMod *num.Modulus) Plaintext {
	res := NewPoly(len(x))
	vec.ReduceTo(res, x, msgMod)
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
	res := make([]uint64, len(e))
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

// Ciphertext is a BFV ciphertext.
type Ciphertext struct {
	Value *rlwe.Ciphertext
	noise float64
}

// NewCiphertext creates a new [Ciphertext].
func NewCiphertext(params rlwe.Parameters, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Value: rlwe.NewCiphertext(params, false, isNTT),
		noise: 0,
	}
}

// NewCiphertextCustom creates a new [Ciphertext] with custom parameters.
func NewCiphertextCustom(rank, modLen int, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Value: rlwe.NewCiphertextCustom(rank, modLen, 0, isNTT),
		noise: 0,
	}
}

// Noise returns the noise of c.
func (c *Ciphertext) Noise() float64 {
	return c.noise
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
	c.noise = 0
}

// Resize resizes the ciphertext to the given modulus length.
func (c *Ciphertext) Resize(modLen int) {
	c.Value.Resize(modLen, 0)
}

// WithModLen returns a shallow copy with the given modulus lengths.
func (c *Ciphertext) WithModLen(baseLen int) *Ciphertext {
	return &Ciphertext{
		Value: c.Value.WithModLen(baseLen, 0),
		noise: c.noise,
	}
}

// Copy returns a copy of c.
func (c *Ciphertext) Copy() *Ciphertext {
	return &Ciphertext{
		Value: c.Value.Copy(),
		noise: c.noise,
	}
}

// CopyFrom copies the coefficients from cIn to c.
func (c *Ciphertext) CopyFrom(cIn *Ciphertext) {
	c.Value.CopyFrom(cIn.Value)
	c.noise = cIn.noise
}

// IsEqual checks if two values are equal.
func (c *Ciphertext) IsEqual(c0 *Ciphertext) bool {
	return c.Value.IsEqual(c0.Value) && c.noise == c0.noise
}

// IsConsistent checks if two values have the same shape.
func (c *Ciphertext) IsConsistent(c0 *Ciphertext) bool {
	return c.Value.IsConsistent(c0.Value) && c.noise == c0.noise
}
