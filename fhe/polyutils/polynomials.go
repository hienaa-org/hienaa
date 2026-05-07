package polyutils

import (
	"math"

	"github.com/hienaa-org/hienaa/math/num"
)

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

// MaxLevel returns the maximum level of the polynomial for the Paterson-Stockmeyer algorithm.
func MaxLevel(deg int) int {
	return int(math.Ceil(num.Log2(deg + 1)))
}

// MaxDeg returns the maximum degree of the polynomial for the Paterson-Stockmeyer algorithm.
func MaxDeg(deg int) int {
	return (1 << MaxLevel(deg)) - 1
}

// BabyLevel returns the baby level of the polynomial for the Paterson-Stockmeyer algorithm.
func BabyLevel(deg int) int {
	return (MaxLevel(deg) + 1) / 2
}

// BabyDeg returns the baby degree of the polynomial for the Paterson-Stockmeyer algorithm.
func BabyDeg(deg int) int {
	return 1 << BabyLevel(deg)
}
