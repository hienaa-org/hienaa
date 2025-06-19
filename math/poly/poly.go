// Package poly implements polynomial and its operations.
package poly

// Poly represents a polynomial with "Double-CRT" representation.
//
// A polynomial can have two forms:
//
//   - Standard form, where the coefficients are stored in a [][]uint64 slice in normal order.
//   - NTT form, where the NTT coefficients are stored in a [][]uint64 slice in radix-r reversed order.
//     All coefficients in NTT form are also in Montgomery form.
type Poly struct {
	// Coeffs are the coefficients of the polynomial.
	// Ordered as [ModLen][Degree].
	//
	// All subslice of Coeffs are assumed to have the same length,
	// and the length of the first subslice is considered the degree of the polynomial.
	Coeffs [][]uint64

	// IsNTT indicates whether the polynomial is in NTT form.
	IsNTT bool
}

// NewPoly creates a new [Poly] with the given degree and modLen.
func NewPoly(deg, modLen int) *Poly {
	return NewPolyCustom(deg, modLen, false)
}

// NewNTTPoly creates a new [Poly] in NTT form with the given degree and modLen.
func NewNTTPoly(deg, modLen int) *Poly {
	return NewPolyCustom(deg, modLen, true)
}

// NewPolyCustom creates a new [Poly] with the given degree, modLen and NTT flag.
func NewPolyCustom(deg, modLen int, isNTT bool) *Poly {
	coeffs := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		coeffs[i] = make([]uint64, deg)
	}
	return &Poly{
		Coeffs: coeffs,
		IsNTT:  isNTT,
	}
}

// Degree returns the degree of p.
func (p *Poly) Degree() int {
	if len(p.Coeffs) == 0 {
		return 0
	}
	return len(p.Coeffs[0])
}

// ModLen returns the number of RNS moduli of p.
func (p *Poly) ModLen() int {
	return len(p.Coeffs)
}

// Clear clears p.
func (p *Poly) Clear() {
	for i := range p.Coeffs {
		clear(p.Coeffs[i])
	}
}

// Copy returns a copy of p.
func (p *Poly) Copy() *Poly {
	pOut := NewPolyCustom(p.Degree(), p.ModLen(), p.IsNTT)
	for i := range p.Coeffs {
		copy(pOut.Coeffs[i], p.Coeffs[i])
	}
	return pOut
}

// CopyFrom copies the coefficients from p0 to p.
// Panics when p and p0 are not consistent.
func (p *Poly) CopyFrom(p0 *Poly) {
	if !p.IsConsistent(p0) {
		panic("CopyFrom: inconsistent polynomials")
	}

	for i := range p.Coeffs {
		copy(p.Coeffs[i], p0.Coeffs[i])
	}
}

// IsEqual checks if p is equal to p0.
func (p *Poly) IsEqual(p0 *Poly) bool {
	if !p.IsConsistent(p0) {
		return false
	}

	for i := range p.Coeffs {
		for j := range p.Coeffs[i] {
			if p.Coeffs[i][j] != p0.Coeffs[i][j] {
				return false
			}
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
