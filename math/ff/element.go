package ff

import "github.com/hienaa-org/hienaa/math/crt"

// Element represents an element in a finite field.
type Element struct {
	poly *crt.Poly
}

// NewElement creates a new [Element].
func NewElement(rank int) *Element {
	return &Element{
		poly: crt.NewNTTPoly(rank, 1),
	}
}

// Rank returns the rank of e.
func (e *Element) Rank() int {
	return e.poly.Rank()
}

// Clear clears e.
func (e *Element) Clear() {
	e.poly.Clear()
}

// Coeffs returns the underlying coefficients of e.
func (e *Element) Coeffs() []uint64 {
	return e.poly.Coeffs[0]
}

// Copy returns a copy of e.
func (e *Element) Copy() *Element {
	return &Element{
		poly: e.poly.Copy(),
	}
}

// CopyFrom copies the coefficients from e0 to e.
// Panics when e and e0 are not consistent.
func (e *Element) CopyFrom(e0 *Element) {
	e.poly.CopyFrom(e0.poly)
}

// IsEqual checks if e is equal to e0.
func (e *Element) IsEqual(e0 *Element) bool {
	return e.poly.IsEqual(e0.poly)
}
