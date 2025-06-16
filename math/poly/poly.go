// Package poly implements polynomial and its operations.
package poly

// Poly represents a polynomial with "Double-CRT" representation.
//
// A polynomial can have two forms:
//
//   - Standard form, where the coefficients are stored in a [][]uint64 slice.
//   - NTT form, where the NTT coefficients are stored in a [][]uint64 slice.
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

// NewPoly creates a new polynomial with the given degree and rnsLen.
func NewPoly(N, rnsLen int) *Poly {
	return NewPolyCustom(N, rnsLen, false)
}

// NewNTTPoly creates a new polynomial in NTT form with the given degree and rnsLen.
func NewNTTPoly(N, rnsLen int) *Poly {
	return NewPolyCustom(N, rnsLen, true)
}

// NewPolyCustom creates a new polynomial with the given degree, rnsLen and NTT flag.
func NewPolyCustom(N, rnsLen int, isNTT bool) *Poly {
	coeffs := make([][]uint64, rnsLen)
	for i := 0; i < rnsLen; i++ {
		coeffs[i] = make([]uint64, N)
	}
	return &Poly{
		Coeffs: coeffs,
		IsNTT:  isNTT,
	}
}

// N returns the ring degree of p.
func (p *Poly) N() int {
	if len(p.Coeffs) == 0 {
		return 0
	}
	return len(p.Coeffs[0])
}

// RNSLen returns the number of RNS moduli of p.
func (p *Poly) RNSLen() int {
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
	pOut := NewPolyCustom(p.N(), p.RNSLen(), p.IsNTT)
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
