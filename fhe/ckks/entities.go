package ckks

import (
	"slices"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
)

// ValueType indicates whether the plaintext is real or integer.
type ValueType int

const (
	TypeReal ValueType = iota
	TypeInt
)

// Plaintext is a CKKS plaintext.
type Plaintext struct {
	Value     []float64
	valueType ValueType
}

// NewScalar creates a new [Plaintext] for a scalar.
func NewScalar(vType ValueType) *Plaintext {
	return &Plaintext{
		Value:     make([]float64, 1),
		valueType: vType,
	}
}

// NewScalarFrom creates a new [Plaintext] for a scalar from x.
func NewScalarFrom(x float64, vType ValueType) *Plaintext {
	res := NewScalar(vType)
	res.Value[0] = x
	return res
}

// NewPoly creates a new [Plaintext] for a polynomial.
func NewPoly(rank int, vType ValueType) *Plaintext {
	return &Plaintext{
		Value:     make([]float64, rank),
		valueType: vType,
	}
}

// NewPolyFrom creates a new [Plaintext] for a polynomial from x.
func NewPolyFrom(x []float64, vType ValueType) *Plaintext {
	res := NewPoly(len(x), vType)
	copy(res.Value, x)
	return res
}

// Rank returns the rank.
func (e Plaintext) Rank() int {
	return len(e.Value)
}

// Clear clears value.
func (e Plaintext) Clear() {
	clear(e.Value)
}

// Copy returns a copy of e.
func (e Plaintext) Copy() *Plaintext {
	res := make([]float64, len(e.Value))
	copy(res, e.Value)
	return &Plaintext{
		Value:     res,
		valueType: e.valueType,
	}
}

// CopyFrom copies the coefficients from eIn to e.
func (e Plaintext) CopyFrom(eIn *Plaintext) {
	copy(e.Value, eIn.Value)
}

// IsEqual checks if two values are equal.
func (e Plaintext) IsEqual(e0 *Plaintext) bool {
	return slices.Equal(e.Value, e0.Value) && e.valueType == e0.valueType
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

// ScalingFactor returns the scaling factor.
func (c *Ciphertext) ScalingFactor() float64 {
	return c.scFac
}

// SetScalingFactor sets the scaling factor.
func (c *Ciphertext) SetScalingFactor(scFac float64) {
	c.scFac = scFac
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

// Resize resizes the ciphertext to the given modulus length.
func (c *Ciphertext) Resize(modLen int) {
	c.Value.Resize(modLen, 0)
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
