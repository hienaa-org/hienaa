package poly

import "github.com/hienaa-org/hienaa/math/num"

// RingType is a type of the polynomial ring.
type RingType uint64

const (
	// Cyclotomic is a cyclotomic ring ZZ[X]/Phi_M(X).
	Cyclotomic RingType = iota
	// Cyclic is a cyclic ring ZZ[X]/(X^N - 1).
	Cyclic
	// Decomposition is a decomposition ring of a cyclotomic ring.
	// In other words, it is a subring of a cyclotomic ring invariant under some automorphism.
	Decomposition
)

// RingParameters contains the parameters for the ring.
type RingParameters struct {
	// CycloDegree is the degree of the underlying cyclotomic polynomial.
	// 0 if the RingType is not [Cyclotomic] or [Decomposition].
	cycloDegree int

	// Degree is the number of coefficients of the polynomial in the ring.
	degree int

	// foldFactor is the folding factor of the decomposition ring.
	// 0 if the RingType is not [Decomposition].
	foldFactor int

	// ringType is the type of the ring.
	ringType RingType
}

// NewCyclotomicParameters creates a new [RingParameters] for a cyclotomic ring.
func NewCyclotomicParameters(cycloDeg int) RingParameters {
	if cycloDeg <= 0 {
		panic("NewCyclotomicParameters: cycloDegree must be positive")
	}

	degree := int(num.EulerPhi(uint64(cycloDeg)))
	return RingParameters{
		cycloDegree: cycloDeg,
		degree:      degree,
	}
}

// CycloDegree is the degree of the underlying cyclotomic polynomial.
// 0 if the RingType is not [Cyclotomic] or [Decomposition].
func (p RingParameters) CycloDegree() int {
	return p.cycloDegree
}

// Degree is the number of coefficients of the polynomial in the ring.
func (p RingParameters) Degree() int {
	return p.degree
}

// FoldFactor is the folding factor of the decomposition ring.
// 0 if the RingType is not [Decomposition].
func (p RingParameters) FoldFactor() int {
	return p.foldFactor
}

// Type is the type of the ring.
func (p RingParameters) Type() RingType {
	return p.ringType
}
