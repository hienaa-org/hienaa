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

	// IsNTT indicates whether the polynomial is in NTT form.
	IsNTT bool
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
		IsNTT:  isNTT,
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

// Clear clears p.
func (p *Poly) Clear() {
	for i := range p.Coeffs {
		clear(p.Coeffs[i])
	}
}

// SetCoeff sets the i-th coefficient to c.
func (p *Poly) SetCoeff(i int, c Scalar) {
	if len(c) != p.ModLen() {
		panic("SetCoeff: inconsistent modulus length")
	}

	for j := range p.Coeffs {
		p.Coeffs[j][i] = c[j]
	}
}

// WithModIdx returns a copy of p with the given modulus indices.
//
// Panics when idx is out of range.
func (p *Poly) WithModIdx(idx ...int) *Poly {
	for _, idxi := range idx {
		if idxi < 0 || idxi >= p.ModLen() {
			panic("WithModIdx: index out of range")
		}
	}

	coeffs := make([][]uint64, len(idx))
	for i, idxi := range idx {
		coeffs[i] = p.Coeffs[idxi]
	}

	return &Poly{
		Coeffs: coeffs,
		IsNTT:  p.IsNTT,
	}
}

// Copy returns a copy of p.
func (p *Poly) Copy() *Poly {
	pOut := NewPolyCustom(p.Rank(), p.ModLen(), p.IsNTT)
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
		panic("CopyFrom: inconsistent polynomials")
	}

	for i := range p.Coeffs {
		copy(p.Coeffs[i], pIn.Coeffs[i])
	}
	p.IsNTT = pIn.IsNTT
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
