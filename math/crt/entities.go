package crt

import (
	"math/big"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
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

// WithModIdx returns a copy of c with the given modulus indices.
//
// Panics when idx is out of range.
func (c *Scalar) WithModIdx(idx ...int) *Scalar {
	for i := range idx {
		if !(0 <= idx[i] && idx[i] < c.ModLen()) {
			panic("index out of range")
		}
	}

	return &Scalar{
		Value: vec.Gather(c.Value, idx...),
	}
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

// CopyFrom copies the coefficients from c0 to c.
//
// Panics when c and c0 are not consistent.
func (c *Scalar) CopyFrom(c0 *Scalar) {
	if !c.IsConsistent(c0) {
		panic("inconsistent input(s)")
	}

	copy(c.Value, c0.Value)
}

// IsEqual checks if c is equal to c0.
func (c *Scalar) IsEqual(c0 *Scalar) bool {
	return slices.Equal(c.Value, c0.Value)
}

// IsConsistent checks if c has the same shape as c0.
func (c *Scalar) IsConsistent(c0 *Scalar) bool {
	return c.ModLen() == c0.ModLen()
}

// Poly is a polynomial in CRT representation.
// It can have Standard or NTT form.
// All coefficients in NTT form are also in Montgomery form.
type Poly struct {
	// Coeffs are the coefficients of the polynomial.
	// Ordered as [ModLen][Rank].
	Coeffs [][]uint64

	// Form represents the [dft.Form] of the polynomial.
	Form dft.Form
}

// NewPoly creates a new [Poly] in Coefficient form.
func NewPoly(rank, modLen int) *Poly {
	return NewPolyCustom(rank, modLen, dft.FormCoeff)
}

// NewNTTPoly creates a new [Poly] in NTT form.
func NewNTTPoly(rank, modLen int) *Poly {
	return NewPolyCustom(rank, modLen, dft.FormNTT)
}

// NewPolyCustom creates a new [Poly].
func NewPolyCustom(rank, modLen int, form dft.Form) *Poly {
	coeffs := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		coeffs[i] = make([]uint64, rank)
	}

	return &Poly{
		Coeffs: coeffs,
		Form:   form,
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

// WithModIdx returns a shallow copy of p with the given modulus indices.
//
// Panics when idx is out of range.
func (p *Poly) WithModIdx(idx ...int) *Poly {
	for i := range idx {
		if !(0 <= idx[i] && idx[i] < p.ModLen()) {
			panic("index out of range")
		}
	}

	coeffs := vec.Gather(p.Coeffs, idx...)

	return &Poly{
		Coeffs: coeffs,
		Form:   p.Form,
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
	pOut := NewPolyCustom(p.Rank(), p.ModLen(), p.Form)
	for i := range p.Coeffs {
		copy(pOut.Coeffs[i], p.Coeffs[i])
	}
	return pOut
}

// CopyFrom copies the coefficients from p0 to p.
//
// Panics when p and p0 are not consistent.
func (p *Poly) CopyFrom(p0 *Poly) {
	if !p.IsConsistent(p0) {
		panic("inconsistent input(s)")
	}

	for i := range p.Coeffs {
		copy(p.Coeffs[i], p0.Coeffs[i])
	}
	p.Form = p0.Form
}

// IsEqual checks if p is equal to p0.
func (p *Poly) IsEqual(p0 *Poly) bool {
	if !p.IsConsistent(p0) {
		return false
	}

	if p.Form != p0.Form {
		return false
	}

	for i := range p.Coeffs {
		if !slices.Equal(p.Coeffs[i], p0.Coeffs[i]) {
			return false
		}
	}

	return true
}

// IsConsistent checks if p has the same shape as p0.
func (p *Poly) IsConsistent(p0 *Poly) bool {
	if len(p.Coeffs) != len(p0.Coeffs) {
		return false
	}

	for i := range p.Coeffs {
		if len(p.Coeffs[i]) != len(p0.Coeffs[i]) {
			return false
		}
	}

	return true
}
