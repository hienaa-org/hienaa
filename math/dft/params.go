package dft

import (
	"math"

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
	// CycloOrder is the order of the underlying cyclotomic polynomial.
	// 0 if the RingType is not [Cyclotomic] or [AutFixed].
	cycloOrd int

	// Rank is the number of coefficients of the polynomial in the ring.
	rank int

	// ringType is the type of the ring.
	ringType RingType
}

// NewCyclotomicParameters creates a new [RingParameters] for a cyclotomic ring.
func NewCyclotomicParameters(cycloOrd int) RingParameters {
	if cycloOrd < 1 {
		panic("NewCyclotomicParameters: cycloOrder must be larger or equal than 1")
	}

	return RingParameters{
		cycloOrd: cycloOrd,
		rank:     int(num.Totient(uint64(cycloOrd))),
		ringType: Cyclotomic,
	}
}

// NewCyclicParameters creates a new [RingParameters] for a cyclic ring.
func NewCyclicParameters(rank int) RingParameters {
	if rank < 1 {
		panic("NewCyclicParameters: rank must be larger or equal than 1")
	}

	return RingParameters{
		cycloOrd: 0,
		rank:     rank,
		ringType: Cyclic,
	}
}

// NewAutFixedParameters creates a new [RingParameters] for an AutFixed ring.
func NewAutFixedParameters(cycloOrd, rank int) RingParameters {
	if rank < 1 {
		panic("NewAutFixedParameters: rank must be larger or equal than 1")
	}

	if num.IsPrime(uint64(cycloOrd)) {
		if (cycloOrd-1)%rank != 0 {
			panic("NewAutFixedParameters: rank should divide cycloOrd-1 for prime cycloOrd")
		}
	} else if num.IsPowerOfTwo(uint64(cycloOrd)) {
		if cycloOrd != rank<<2 {
			panic("NewAutFixedParameters: cycloOrd must be four times the rank for power-of-two cycloOrd")
		}
	} else {
		panic("NewAutFixedParameters: cycloOrd must be a prime or a power of two")
	}

	return RingParameters{
		cycloOrd: cycloOrd,
		rank:     rank,
		ringType: AutFixed,
	}
}

// CycloOrder is the order of the underlying cyclotomic polynomial.
// 0 if the RingType is [Cyclic].
func (p RingParameters) CycloOrder() int {
	return p.cycloOrd
}

// Rank is the number of coefficients of the polynomial in the ring.
func (p RingParameters) Rank() int {
	return p.rank
}

// RingType is the type of the ring.
func (p RingParameters) RingType() RingType {
	return p.ringType
}

// cyclotomicGap finds the "gap" of the NTT-friendly modulus for cyclotomic rings.
func cyclotomicGap(cycloOrd, rank uint64) uint64 {
	var gap uint64

	if num.IsPowerOfTwo(cycloOrd) {
		gap = cycloOrd
	} else {
		bluesteinRank := num.NextProdPower(2*cycloOrd-1, []uint64{2})

		var p uint64
		primes, _ := num.Factor(cycloOrd)
		for _, f := range primes {
			if f%2 == 1 {
				p = f
				break
			}
		}

		redDeg := cycloOrd - cycloOrd/p
		if redDeg == rank {
			gap = num.LCM(cycloOrd, bluesteinRank)
		} else {
			degNext := num.NextProdPower(rank, []uint64{2})
			diffDegNext := num.NextProdPower(2*(redDeg-rank)+1, []uint64{2})
			gap = num.LCM(num.LCM(degNext, diffDegNext), num.LCM(cycloOrd, bluesteinRank))
		}
	}

	return gap
}

// cyclicGap finds the "gap" of the NTT-friendly modulus for cyclic rings.
func cyclicGap(rank uint64) uint64 {
	var gap uint64

	if num.IsProdPowerOf(rank, cyclicNTTFactors) {
		gap = rank
	} else {
		gap = num.LCM(num.NextProdPower(2*rank-1, []uint64{2}), rank)
	}

	return gap
}

// autFixedGap finds the "gap" of the NTT-friendly modulus for AutFixed rings.
func autFixedGap(cycloOrd, rank uint64) uint64 {
	var gap uint64

	if num.IsPowerOfTwo(cycloOrd) {
		gap = cycloOrd
	} else if num.IsProdPowerOf(rank, []uint64{2}) {
		gap = cycloOrd * rank
	} else {
		gap = cycloOrd * num.NextProdPower(2*rank-1, []uint64{2})
	}

	return gap
}

// IsNTTFriendly checks if the given modulus is NTT-friendly with respect to the ring parameters.
func IsNTTFriendly(ringParams RingParameters, mod *num.Modulus) bool {
	var gap uint64

	switch ringParams.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	case Cyclic:
		gap = cyclicGap(uint64(ringParams.rank))
	case AutFixed:
		gap = autFixedGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	}

	primes, _ := num.Factor(mod.Value())
	for _, f := range primes {
		if f%gap != 1 {
			return false
		}
	}
	return true
}

// FindNextNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes greater than or equal to 2^bits.
func FindNextNTTPrimes(ringParams RingParameters, bits float64, cnt int) []*num.Modulus {
	var gap uint64

	switch ringParams.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	case Cyclic:
		gap = cyclicGap(uint64(ringParams.rank))
	case AutFixed:
		gap = autFixedGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	}

	start := (uint64(math.Round(math.Exp2(bits)))/gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	prime := num.NextPrime(start, gap)
	for i := 0; i < cnt; i++ {
		primes[i] = num.NewModulus(prime)
		prime = num.NextPrime(prime, gap)
	}
	return primes
}

// FindPrevNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes less than or equal to 2^bits.
func FindPrevNTTPrimes(ringParams RingParameters, bits float64, cnt int) []*num.Modulus {
	var gap uint64

	switch ringParams.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	case Cyclic:
		gap = cyclicGap(uint64(ringParams.rank))
	case AutFixed:
		gap = autFixedGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	}

	start := (uint64(math.Floor(math.Exp2(bits)))/gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	prime := num.PrevPrime(start, gap)
	for i := 0; i < cnt; i++ {
		primes[i] = num.NewModulus(prime)
		prime = num.PrevPrime(prime, gap)
	}
	return primes
}

// FindNearestNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes nearest to 2^bits.
// Output modulus are alternating in size.
func FindNearestNTTPrimes(ringParams RingParameters, bits float64, cnt int) []*num.Modulus {
	var gap uint64

	switch ringParams.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	case Cyclic:
		gap = cyclicGap(uint64(ringParams.rank))
	case AutFixed:
		gap = autFixedGap(uint64(ringParams.cycloOrd), uint64(ringParams.rank))
	}

	start := (uint64(math.Round(math.Exp2(bits)))/gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	prime := num.NextPrime(start, gap)
	for i := 0; i < cnt; i += 2 {
		primes[i] = num.NewModulus(prime)
		prime = num.NextPrime(prime, gap)
	}

	prime = num.PrevPrime(start, gap)
	for i := 1; i < cnt; i += 2 {
		primes[i] = num.NewModulus(prime)
		prime = num.PrevPrime(prime, gap)
	}
	return primes
}
