package polyutils

import "github.com/hienaa-org/hienaa/math/num"

// PolynomialType is the type of polynomial.
type PolynomialType int

const (
	Monomial PolynomialType = iota
	Chebyshev
)

// Polynomial is a generic interface for polynomials.
type Polynomial[T num.Number] interface {
	// Degree returns the degree of the polynomial.
	Degree() int
	// Type returns the type of the polynomial.
	Type() PolynomialType
	// Coeffs returns the coefficients of the polynomial.
	Coeffs() map[int]T
}
