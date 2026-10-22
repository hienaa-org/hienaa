package crt

import (
	"math/big"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Scalar is a scalar in CRT representation.
type Scalar struct {
	// Value is a scalar reduced by each modulus.
	Value []uint64
}

// NewScalar creates a new [Scalar].
func NewScalar(modLen int) *Scalar {
	return &Scalar{
		Value: make([]uint64, modLen),
	}
}

// NewScalarFrom creates a new [Scalar] from x.
func NewScalarFrom[T *big.Int | num.Integer](x T, mod []*num.Modulus) *Scalar {
	r := NewScalar(len(mod))

	switch x := any(x).(type) {
	case *big.Int:
		q, t := new(big.Int), new(big.Int)
		for i := range r.Value {
			q.SetUint64(mod[i].Value())
			r.Value[i] = t.Mod(x, q).Uint64()
		}
	case uint:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case uint8:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case uint16:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case uint32:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case uint64:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case int:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case int8:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case int16:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case int32:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	case int64:
		for i := range r.Value {
			r.Value[i] = num.Reduce(x, mod[i])
		}
	}

	return r
}

// ModLen returns the length of moduli of c.
func (c *Scalar) ModLen() int {
	return len(c.Value)
}

// Clear clears c.
func (c *Scalar) Clear() {
	clear(c.Value)
}

// Slice slices the modulus of [Scalar] and returns the new [Scalar].
func (c *Scalar) Slice(lo, hi int) *Scalar {
	return &Scalar{
		Value: c.Value[lo:hi],
	}
}

// Copy returns a copy of c.
func (c *Scalar) Copy() *Scalar {
	cOut := NewScalar(c.ModLen())
	copy(cOut.Value, c.Value)
	return cOut
}

// CopyFrom copies the coefficients from cIn to c.
//
// Panics when c and cIn are not consistent.
func (c *Scalar) CopyFrom(cIn *Scalar) {
	if !c.IsConsistent(cIn) {
		panic("inconsistent inputs")
	}

	copy(c.Value, cIn.Value)
}

// IsEqual checks if c is equal to cIn.
func (c *Scalar) IsEqual(cIn *Scalar) bool {
	return slices.Equal(c.Value, cIn.Value)
}

// IsConsistent checks if c has the same shape as cIn.
func (c *Scalar) IsConsistent(cIn *Scalar) bool {
	return c.ModLen() == cIn.ModLen()
}

// Poly is a polynomial in CRT representation.
// It can have Standard or NTT form.
type Poly struct {
	// Coeffs are the coefficients of the polynomial.
	// Ordered as [ModLen][Rank].
	Coeffs [][]uint64

	// Form represents the [dft.Form] of the polynomial.
	Form dft.Form
}

// NewPoly creates a new [Poly].
func NewPoly(rank, modLen int, form dft.Form) *Poly {
	coeffs := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		coeffs[i] = make([]uint64, rank)
	}

	return &Poly{
		Coeffs: coeffs,
		Form:   form,
	}
}

// SetCoeff sets the i-th coefficient to c.
func (p *Poly) SetCoeff(i int, c *Scalar) {
	if len(p.Coeffs) != len(c.Value) {
		panic("inconsistent inputs")
	}

	for j := range p.Coeffs {
		p.Coeffs[j][i] = c.Value[j]
	}
}

// Rank returns the length of the coefficients of p.
func (p *Poly) Rank() int {
	if len(p.Coeffs) == 0 {
		return 0
	}

	rank := len(p.Coeffs[0])
	for i := 1; i < len(p.Coeffs); i++ {
		if len(p.Coeffs[i]) != rank {
			panic("inconsistent rank")
		}
	}

	return rank
}

// ModLen returns the length of moduli of p.
func (p *Poly) ModLen() int {
	return len(p.Coeffs)
}

// Clear clears p.
func (p *Poly) Clear() {
	for i := range p.Coeffs {
		clear(p.Coeffs[i])
	}
}

// Slice slices the modulus of [Poly] and returns the new [Poly].
func (p *Poly) Slice(lo, hi int) *Poly {
	return &Poly{
		Coeffs: p.Coeffs[lo:hi],
		Form:   p.Form,
	}
}

// Copy returns a copy of p.
func (p *Poly) Copy() *Poly {
	pOut := NewPoly(p.Rank(), p.ModLen(), p.Form)
	for i := range p.Coeffs {
		copy(pOut.Coeffs[i], p.Coeffs[i])
	}
	return pOut
}

// CopyFrom copies the coefficients from pIn to p.
//
// Panics when p and pIn are not consistent.
func (p *Poly) CopyFrom(pIn *Poly) {
	if !p.IsConsistent(pIn) {
		panic("inconsistent inputs")
	}

	for i := range p.Coeffs {
		copy(p.Coeffs[i], pIn.Coeffs[i])
	}
	p.Form = pIn.Form
}

// IsEqual checks if p is equal to pIn.
func (p *Poly) IsEqual(pIn *Poly) bool {
	if !p.IsConsistent(pIn) {
		return false
	}

	if p.Form != pIn.Form {
		return false
	}

	for i := range p.Coeffs {
		if !slices.Equal(p.Coeffs[i], pIn.Coeffs[i]) {
			return false
		}
	}

	return true
}

// IsConsistent checks if p has the same shape as pIn.
func (p *Poly) IsConsistent(pIn *Poly) bool {
	if len(p.Coeffs) != len(pIn.Coeffs) {
		return false
	}

	for i := range p.Coeffs {
		if len(p.Coeffs[i]) != len(pIn.Coeffs[i]) {
			return false
		}
	}

	return true
}
