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

// cyclotomicGap finds the "gap" of the NTT-friendly modulus for cyclotomic rings.
func cyclotomicGap(cycloDeg, deg uint64) uint64 {
	var gap uint64

	if num.IsPowerOfTwo(cycloDeg) {
		gap = cycloDeg
	} else {
		bluesteinDeg := num.NextProdPower(2*cycloDeg-1, []uint64{2})

		var p uint64
		for f := range num.Factor(cycloDeg) {
			if f%2 == 1 {
				p = f
				break
			}
		}

		redDeg := cycloDeg - cycloDeg/p
		if redDeg == deg {
			gap = num.LCM(cycloDeg, bluesteinDeg)
		} else {
			degNext := num.NextProdPower(deg, []uint64{2})
			diffDegNext := num.NextProdPower(2*(redDeg-deg)+1, []uint64{2})
			gap = num.LCM(num.LCM(degNext, diffDegNext), num.LCM(cycloDeg, bluesteinDeg))
		}
	}

	return gap
}

// cyclicGap finds the "gap" of the NTT-friendly modulus for cyclic rings.
func cyclicGap(deg uint64) uint64 {
	var gap uint64

	if num.IsProdPowerOf(deg, cyclicNTTFactors) {
		gap = deg
	} else {
		gap = num.LCM(num.NextProdPower(2*deg-1, cyclicNTTFactors), deg)
	}

	return gap
}

// autFixedGap finds the "gap" of the NTT-friendly modulus for AutFixed rings.
func autFixedGap(cycloDeg, deg uint64) uint64 {
	var gap uint64

	if num.IsProdPowerOf(deg, cyclicNTTFactors) {
		gap = num.LCM(cycloDeg, deg)
	} else {
		gap = num.LCM(num.NextProdPower(2*deg-1, cyclicNTTFactors), cycloDeg)
	}

	return gap
}

// isNTTFriendly checks if the given modulus is NTT-friendly with respect to the ring parameters.
func isNTTFriendly(ringParams RingParameters, modulus *mod.Modulus) bool {
	var gap uint64

	switch ringParams.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(uint64(ringParams.cycloDegree), uint64(ringParams.degree))
	case Cyclic:
		gap = cyclicGap(uint64(ringParams.degree))
	case AutFixed:
		gap = autFixedGap(uint64(ringParams.cycloDegree), uint64(ringParams.degree))
	}

	for f := range num.Factor(modulus.Value()) {
		if f%gap != 1 {
			return false
		}
	}
	return true
}

// FindNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
func FindNTTPrimes(ringParams RingParameters, start uint64, cnt int) []*mod.Modulus {
	var gap uint64

	switch ringParams.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(uint64(ringParams.cycloDegree), uint64(ringParams.degree))
	case Cyclic:
		gap = cyclicGap(uint64(ringParams.degree))
	case AutFixed:
		gap = autFixedGap(uint64(ringParams.cycloDegree), uint64(ringParams.degree))
	}

	start = (start/gap)*gap + gap + 1
	primes := make([]*mod.Modulus, cnt)
	primes[0] = mod.NewModulus(num.NextPrime(start, gap))
	for i := 1; i < cnt; i++ {
		primes[i] = mod.NewModulus(num.NextPrime(primes[i-1].Value(), gap))
	}

	return primes
}
