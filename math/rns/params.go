package rns

import (
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

// RingType is a type of the polynomial ring.
type RingType uint64

const (
	// Cyclotomic is a cyclotomic ring ZZ[X]/Phi_M(X).
	Cyclotomic RingType = iota
	// Cyclic is a cyclic ring ZZ[X]/(X^N - 1).
	Cyclic
	// AutFixed is a decomposition ring of a cyclotomic ring.
	// In other words, it is a subring of a cyclotomic ring invariant under some automorphism.
	AutFixed
)

// RingParameters contains the parameters for the ring.
type RingParameters struct {
	// CycloDegree is the degree of the underlying cyclotomic polynomial.
	// 0 if the RingType is not [Cyclotomic] or [AutFixed].
	cycloDegree int

	// Degree is the number of coefficients of the polynomial in the ring.
	degree int

	// ringType is the type of the ring.
	ringType RingType
}

// NewCyclotomicParameters creates a new [RingParameters] for a cyclotomic ring.
func NewCyclotomicParameters(cycloDeg int) RingParameters {
	if cycloDeg <= 0 {
		panic("NewCyclotomicParameters: cycloDegree must be positive")
	}

	return RingParameters{
		cycloDegree: cycloDeg,
		degree:      int(num.Totient(uint64(cycloDeg))),
		ringType:    Cyclotomic,
	}
}

// NewCyclicParameters creates a new [RingParameters] for a cyclic ring.
func NewCyclicParameters(deg int) RingParameters {
	if deg <= 0 {
		panic("NewCyclicParameters: degree must be positive")
	}

	return RingParameters{
		cycloDegree: 0,
		degree:      deg,
		ringType:    Cyclic,
	}
}

// CycloDegree is the degree of the underlying cyclotomic polynomial.
// 0 if the RingType is [Cyclic].
func (p RingParameters) CycloDegree() int {
	return p.cycloDegree
}

// Degree is the number of coefficients of the polynomial in the ring.
func (p RingParameters) Degree() int {
	return p.degree
}

// Type is the type of the ring.
func (p RingParameters) Type() RingType {
	return p.ringType
}

// isNTTFriendly checks if the given modulus is NTT-friendly with respect to the ring parameters.
func isNTTFriendly(ringParams RingParameters, modulus *mod.Modulus) bool {
	switch ringParams.ringType {
	case Cyclotomic:
		return isCyclotomicNTTFriendly(uint64(ringParams.cycloDegree), uint64(ringParams.degree), modulus)
	case Cyclic:
		return isCyclicNTTFriendly(uint64(ringParams.degree), modulus)
	case AutFixed:
		return isAutFixedNTTFriendly(uint64(ringParams.cycloDegree), uint64(ringParams.degree), modulus)
	}
	return false
}

// isCyclotomicNTTFriendly checks if the given modulus is NTT-friendly for Cyclotomic rings.
func isCyclotomicNTTFriendly(cycloDeg, deg uint64, modulus *mod.Modulus) bool {
	var step uint64
	if num.IsPowerOfTwo(cycloDeg) {
		step = cycloDeg
	} else {
		bluesteinDeg := num.NextProdPower(2*cycloDeg-1, cyclicNTTFactors)

		var p uint64
		for f := range num.Factor(cycloDeg) {
			if f%2 == 1 {
				p = f
				break
			}
		}

		redDeg := cycloDeg - cycloDeg/p
		if redDeg == deg {
			step = num.LCM(cycloDeg, bluesteinDeg)
		} else {
			degNext := num.NextProdPower(deg, ambientDegreeFactors)
			diffDegNext := num.NextProdPower(2*(redDeg-deg)+1, ambientDegreeFactors)
			step = num.LCM(num.LCM(degNext, diffDegNext), num.LCM(cycloDeg, bluesteinDeg))
		}
	}

	for f := range num.Factor(modulus.Value()) {
		if f%step != 1 {
			return false
		}
	}
	return true
}

// isCyclicNTTFriendly checks if the given modulus is NTT-friendly for Cyclic rings.
func isCyclicNTTFriendly(deg uint64, modulus *mod.Modulus) bool {
	var step uint64
	if num.IsProdPowerOf(deg, cyclicNTTFactors) {
		step = deg
	} else {
		step = num.LCM(num.NextProdPower(2*deg-1, cyclicNTTFactors), deg)
	}

	for f := range num.Factor(modulus.Value()) {
		if f%step != 1 {
			return false
		}
	}
	return true
}

// isAutFixedNTTFriendly checks if the given modulus is NTT-friendly for AutFixed rings.
func isAutFixedNTTFriendly(cycloDeg, deg uint64, modulus *mod.Modulus) bool {
	var step uint64
	if num.IsProdPowerOf(deg, cyclicNTTFactors) {
		step = num.LCM(cycloDeg, deg)
	} else {
		step = num.LCM(num.NextProdPower(2*deg-1, cyclicNTTFactors), cycloDeg)
	}

	for f := range num.Factor(modulus.Value()) {
		if f%step != 1 {
			return false
		}
	}
	return true
}
