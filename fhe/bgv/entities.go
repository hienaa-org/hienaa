package bgv

import (
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Plaintext is a BGV plaintext embedded in ciphertext modulus.
type Plaintext crt.Element

// NewScalar creates a new [Plaintext] for a scalar.
func NewScalar() *Plaintext {
	return (*Plaintext)(crt.NewScalar(1))
}

// NewScalarFrom creates a new [Plaintext] for a scalar from x.
func NewScalarFrom(x uint64, msgMod *num.Modulus) *Plaintext {
	res := NewScalar()
	res.Coeffs[0][0] = num.Reduce(x, msgMod)
	return res
}

// NewPoly creates a new [Plaintext] for a polynomial.
func NewPoly(rank int) *Plaintext {
	return (*Plaintext)(crt.NewPolyCustom(rank, 1, false))
}

// NewPolyFrom creates a new [Plaintext] for a polynomial from x.
func NewPolyFrom(x []uint64, msgMod *num.Modulus) *Plaintext {
	res := NewPoly(len(x))
	vec.ReduceTo(res.Coeffs[0], x, msgMod)
	return res
}

// Rank returns the rank.
func (e *Plaintext) Rank() int {
	return (*crt.Element)(e).Rank()
}

// Clear clears value.
func (e *Plaintext) Clear() {
	(*crt.Element)(e).Clear()
}

// Copy returns a copy of e.
func (e *Plaintext) Copy() *Plaintext {
	return (*Plaintext)((*crt.Element)(e).Copy())
}

// CopyFrom copies the coefficients from eIn to e.
func (e *Plaintext) CopyFrom(eIn *Plaintext) {
	(*crt.Element)(e).CopyFrom((*crt.Element)(eIn))
}

// IsEqual checks if two values are equal.
func (e *Plaintext) IsEqual(e0 *Plaintext) bool {
	return (*crt.Element)(e).IsEqual((*crt.Element)(e0))
}

// Ciphertext is a BGV ciphertext.
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
