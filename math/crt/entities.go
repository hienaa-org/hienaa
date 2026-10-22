package crt

import (
	"math/big"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// ElementType is the type of [Element].
type ElementType bool

const (
	// TypeScalar is a scalar.
	TypeScalar ElementType = false
	// TypePoly is a polynomial.
	TypePoly ElementType = true
)

// Element represents an element with CRT representation.
// Can be a scalar or a polynomial.
//
// A polynomial can have Standard or NTT form.
type Element struct {
	// Coeffs are the coefficients of the element.
	// Ordered as [ModLen][Rank].
	Coeffs [][]uint64

	// Form represents the [dft.Form] of the polynomial.
	// Ignored when element type is [TypeScalar].
	Form dft.Form

	// Type is the [ElementType].
	Type ElementType
}

// NewScalar creates a new scalar [Element].
func NewScalar(modLen int) *Element {
	coeffs := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		coeffs[i] = make([]uint64, 1)
	}

	return &Element{
		Coeffs: coeffs,
		Type:   TypeScalar,
	}
}

// NewScalarFrom creates a new scalar [Element] from x.
func NewScalarFrom[T *big.Int | num.Integer](x T, mod []*num.Modulus) *Element {
	r := NewScalar(len(mod))

	switch x := any(x).(type) {
	case *big.Int:
		q, t := new(big.Int), new(big.Int)
		for i := range r.Coeffs {
			q.SetUint64(mod[i].Value())
			r.Coeffs[i][0] = t.Mod(x, q).Uint64()
		}
	case uint:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case uint8:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case uint16:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case uint32:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case uint64:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case int:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case int8:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case int16:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case int32:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	case int64:
		for i := range r.Coeffs {
			r.Coeffs[i][0] = num.Reduce(x, mod[i])
		}
	}

	return r
}

// NewPoly creates a new polynomial [Element].
func NewPoly(rank, modLen int, form dft.Form) *Element {
	coeffs := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		coeffs[i] = make([]uint64, rank)
	}

	return &Element{
		Coeffs: coeffs,
		Form:   form,
		Type:   TypePoly,
	}
}

// SetCoeff sets the i-th coefficient to c.
func (e *Element) SetCoeff(i int, c *Element) {
	if !(e.Type == TypePoly && c.Type == TypeScalar) {
		panic("inconsistent types")
	} else if len(e.Coeffs) != len(c.Coeffs) {
		panic("inconsistent inputs")
	}

	for j := range e.Coeffs {
		e.Coeffs[j][i] = c.Coeffs[j][0]
	}
}

// Rank returns the length of the coefficients of e.
// If e is a scalar, it returns 1.
func (e *Element) Rank() int {
	if len(e.Coeffs) == 0 {
		return 0
	}

	rank := len(e.Coeffs[0])
	for i := 1; i < len(e.Coeffs); i++ {
		if len(e.Coeffs[i]) != rank {
			panic("inconsistent rank")
		}
	}

	return rank
}

// ModLen returns the length of moduli of e.
func (e *Element) ModLen() int {
	return len(e.Coeffs)
}

// Clear clears e.
func (e *Element) Clear() {
	for i := range e.Coeffs {
		clear(e.Coeffs[i])
	}
}

// Slice slices the modulus of e and returns the new [Element].
func (e *Element) Slice(lo, hi int) *Element {
	return &Element{
		Coeffs: e.Coeffs[lo:hi],
		Form:   e.Form,
		Type:   e.Type,
	}
}

// Copy returns a copy of e.
func (e *Element) Copy() *Element {
	var eOut *Element
	if e.Type == TypeScalar {
		eOut = NewScalar(e.ModLen())
	} else {
		eOut = NewPoly(e.Rank(), e.ModLen(), e.Form)
	}

	for i := range e.Coeffs {
		copy(eOut.Coeffs[i], e.Coeffs[i])
	}
	return eOut
}

// CopyFrom copies the coefficients from eIn to e.
//
// Panics when e and eIn are not consistent.
func (e *Element) CopyFrom(eIn *Element) {
	if !e.IsConsistent(eIn) {
		panic("inconsistent inputs")
	}

	isBackwards := false
	for i := 1; i < len(eIn.Coeffs); i++ {
		if &e.Coeffs[0] == &eIn.Coeffs[i] {
			isBackwards = true
			break
		}
	}

	if isBackwards {
		for i := len(e.Coeffs) - 1; i >= 0; i-- {
			copy(e.Coeffs[i], eIn.Coeffs[i])
		}
	} else {
		for i := 0; i < len(e.Coeffs); i++ {
			copy(e.Coeffs[i], eIn.Coeffs[i])
		}
	}
	e.Form = eIn.Form
	e.Type = eIn.Type
}

// IsEqual checks if e is equal to eIn.
func (e *Element) IsEqual(eIn *Element) bool {
	if !e.IsConsistent(eIn) {
		return false
	}

	if e.Type != eIn.Type {
		return false
	}

	if e.Type == TypePoly && e.Form != eIn.Form {
		return false
	}

	for i := range e.Coeffs {
		if !slices.Equal(e.Coeffs[i], eIn.Coeffs[i]) {
			return false
		}
	}

	return true
}

// IsConsistent checks if e has the same shape as eIn.
func (e *Element) IsConsistent(eIn *Element) bool {
	if len(e.Coeffs) != len(eIn.Coeffs) {
		return false
	}

	for i := range e.Coeffs {
		if len(e.Coeffs[i]) != len(eIn.Coeffs[i]) {
			return false
		}
	}

	return true
}
