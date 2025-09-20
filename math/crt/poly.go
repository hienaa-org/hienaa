// Package poly implements polynomial and its operations.
package crt

// Poly represents a polynomial with CRT representation.
//
// A polynomial can have Standard or NTT form.
// All coefficients in NTT form are also in Montgomery form.
type Poly struct {
	// Coeffs are the coefficients of the polynomial.
	// Ordered as [ModLen][Rank].
	//
	// All subslice of Coeffs are assumed to have the same length,
	// and the length of the first subslice is considered the rank of the polynomial.
	Coeffs [][]uint64

	// isNTT indicates whether the polynomial is in NTT form.
	isNTT bool
}

// NewPoly creates a new [Poly] in Standard form.
func NewPoly(rank, modLen int) *Poly {
	return NewPolyCustom(rank, modLen, false)
}

// NewNTTPoly creates a new [Poly] in NTT form.
func NewNTTPoly(rank, modLen int) *Poly {
	return NewPolyCustom(rank, modLen, true)
}

// NewPolyCustom creates a new [Poly].
func NewPolyCustom(rank, modLen int, isNTT bool) *Poly {
	coeffs := make([][]uint64, modLen)
	for i := 0; i < modLen; i++ {
		coeffs[i] = make([]uint64, rank)
	}

	return &Poly{
		Coeffs: coeffs,
		isNTT:  isNTT,
	}
}

// Rank returns the length of the coefficients of p.
func (p *Poly) Rank() int {
	if len(p.Coeffs) == 0 {
		return 0
	}
	return len(p.Coeffs[0])
}

// ModLen returns the number of RNS moduli of p.
func (p *Poly) ModLen() int {
	return len(p.Coeffs)
}

// IsNTT returns if p is in NTT form.
func (p *Poly) IsNTT() bool {
	return p.isNTT
}

// Clear clears p.
func (p *Poly) Clear() {
	for i := range p.Coeffs {
		clear(p.Coeffs[i])
	}
}

// Copy returns a copy of p.
func (p *Poly) Copy() *Poly {
	pOut := NewPolyCustom(p.Rank(), p.ModLen(), p.isNTT)
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
	p.isNTT = p0.isNTT
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
