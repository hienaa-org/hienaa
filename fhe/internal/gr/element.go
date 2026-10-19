package gr

import "github.com/hienaa-org/hienaa/math/crt"

// Element represents an element in a finite field.
type Element struct {
	value *crt.Poly
}

// NewElement creates a new [Element].
func NewElement(rank int) *Element {
	return &Element{
		value: crt.NewNTTPoly(rank, 1),
	}
}

// Rank returns the rank of e.
func (e *Element) Rank() int {
	return e.value.Rank()
}

// Clear clears e.
func (e *Element) Clear() {
	e.value.Clear()
}

// Coeffs returns the underlying coefficients of e.
func (e *Element) Coeffs() []uint64 {
	return e.value.Coeffs[0]
}

// Copy returns a copy of e.
func (e *Element) Copy() *Element {
	return &Element{
		value: e.value.Copy(),
	}
}

// CopyFrom copies the coefficients from e0 to e.
//
// Panics when e and eIn are not consistent.
func (e *Element) CopyFrom(eIn *Element) {
	e.value.CopyFrom(eIn.value)
}

// IsEqual checks if e is equal to e0.
func (e *Element) IsEqual(e0 *Element) bool {
	return e.value.IsEqual(e0.value)
}
