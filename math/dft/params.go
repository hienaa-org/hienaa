package dft

import (
	"errors"
	"math"
	"slices"

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
	// Other covers arbitrary quotient rings.
	Other
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
	if cycloOrd <= 0 {
		panic("cycloOrd must be positive")
	}

	return RingParameters{
		cycloOrd: cycloOrd,
		rank:     int(num.Totient(uint64(cycloOrd))),
		ringType: Cyclotomic,
	}
}

// NewCyclicParameters creates a new [RingParameters] for a cyclic ring.
func NewCyclicParameters(rank int) RingParameters {
	if rank <= 0 {
		panic("rank must be positive")
	}

	return RingParameters{
		cycloOrd: 0,
		rank:     rank,
		ringType: Cyclic,
	}
}

// NewAutFixedParameters creates a new [RingParameters] for an autfixed ring.
func NewAutFixedParameters(cycloOrd, rank int) RingParameters {
	if rank <= 0 {
		panic("rank must be positive")
	}

	switch {
	case num.IsPowerOfTwo(cycloOrd):
		if cycloOrd != 4*rank {
			panic("cycloOrd must be four times the rank for power-of-two cycloOrd")
		}
	case num.IsPrime(cycloOrd):
		if (cycloOrd-1)%rank != 0 {
			panic("rank should divide cycloOrd-1 for prime cycloOrd")
		}
	default:
		panic("cycloOrd must be a prime or a power of two")
	}

	return RingParameters{
		cycloOrd: cycloOrd,
		rank:     rank,
		ringType: AutFixed,
	}
}

// NewOtherParameters creates a new [RingParameters] for arbitrary quotient ring.
func NewOtherParameters(modPoly []int64) RingParameters {
	if len(modPoly) == 0 {
		panic("modPoly must be non-empty")
	}

	return RingParameters{
		cycloOrd: 0,
		rank:     len(modPoly) - 1,
		ringType: Other,
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
func cyclotomicGap(cycloOrd, rank int) uint64 {
	var gap int

	if num.IsPowerOfTwo(cycloOrd) {
		gap = cycloOrd
	} else {
		bluesteinRank := num.NextProdPower(2*cycloOrd-1, []int{2})

		primes, _ := num.Factor(cycloOrd)
		var redDeg int
		if cycloOrd%2 == 1 {
			redDeg = cycloOrd - cycloOrd/primes[0]
		} else {
			redDeg = cycloOrd/2 - (cycloOrd/2)/primes[1]
		}

		if redDeg == rank {
			gap = num.LCM(cycloOrd, bluesteinRank)
		} else {
			degNext := num.NextProdPower(rank, []int{2})
			diffDegNext := num.NextProdPower(2*(redDeg-rank)+1, []int{2})
			gap = num.LCM(num.LCM(degNext, diffDegNext), num.LCM(cycloOrd, bluesteinRank))
		}
	}

	return uint64(gap)
}

// cyclicGap finds the "gap" of the NTT-friendly modulus for cyclic rings.
func cyclicGap(rank int) uint64 {
	var gap int

	if num.IsProdPowerOf(rank, cyclicNTTFactors) {
		gap = rank
	} else {
		gap = num.LCM(num.NextProdPower(2*rank-1, []int{2}), rank)
	}

	return uint64(gap)
}

// autFixedGap finds the "gap" of the NTT-friendly modulus for autfixed rings.
func autFixedGap(cycloOrd, rank int) uint64 {
	var gap int

	if num.IsPowerOfTwo(cycloOrd) {
		gap = cycloOrd
	} else if num.IsProdPowerOf(rank, []int{2}) {
		gap = cycloOrd * rank
	} else {
		gap = cycloOrd * num.NextProdPower(2*rank-1, []int{2})
	}

	return uint64(gap)
}

func RingGap(params RingParameters) (uint64, error) {
	var gap uint64

	switch params.ringType {
	case Cyclotomic:
		gap = cyclotomicGap(params.cycloOrd, params.rank)
	case Cyclic:
		gap = cyclicGap(params.rank)
	case AutFixed:
		gap = autFixedGap(params.cycloOrd, params.rank)
	default:
		return 0, errors.New("RingGap: invalid ring type")
	}

	return gap, nil
}

// IsNTTFriendly checks if the given modulus is NTT-friendly with respect to the ring parameters.
func IsNTTFriendly(params RingParameters, mod *num.Modulus) bool {
	gap, err := RingGap(params)
	if err != nil {
		return false
	}

	primes, _ := num.Factor(mod.Value())
	for _, f := range primes {
		if f%gap != 1 {
			return false
		}
	}
	return true
}

func MustFindNextNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	primes, err := FindNextNTTPrimes(params, bits, cnt)
	if err != nil {
		panic(err)
	}
	return primes
}

// FindNextNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes greater than or equal to 2^bits.
func FindNextNTTPrimes(params RingParameters, bits float64, cnt int) ([]*num.Modulus, error) {
	gap, err := RingGap(params)
	if err != nil {
		return nil, err
	}

	start := (uint64(math.Round(math.Exp2(bits)))/gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	prime, err := num.NextPrime(start, gap)
	if err != nil {
		return nil, err
	}

	for i := 0; i < cnt; i++ {
		primes[i] = num.NewModulus(prime)
		prime, err = num.NextPrime(prime, gap)
		if err != nil {
			return nil, err
		}
	}
	return primes, nil
}

func MustFindPrevNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	primes, err := FindPrevNTTPrimes(params, bits, cnt)
	if err != nil {
		panic(err)
	}
	return primes
}

// FindPrevNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes less than or equal to 2^bits.
func FindPrevNTTPrimes(params RingParameters, bits float64, cnt int) ([]*num.Modulus, error) {
	gap, err := RingGap(params)
	if err != nil {
		return nil, err
	}

	start := (uint64(math.Floor(math.Exp2(bits)))/gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	prime, err := num.PrevPrime(start, gap)
	if err != nil {
		return nil, err
	}
	for i := 0; i < cnt; i++ {
		primes[i] = num.NewModulus(prime)
		prime, err = num.PrevPrime(prime, gap)
		if err != nil {
			return nil, err
		}
	}
	return primes, nil
}

func MustFindNearestNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	primes, err := FindNearestNTTPrimes(params, bits, cnt)
	if err != nil {
		panic(err)
	}
	return primes
}

// FindNearestNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes nearest to 2^bits.
// Output moduli are alternating in size.
func FindNearestNTTPrimes(params RingParameters, bits float64, cnt int) ([]*num.Modulus, error) {
	gap, err := RingGap(params)
	if err != nil {
		return nil, err
	}

	// Sample half of the primes from the larger side.
	start := (uint64(math.Round(math.Exp2(bits)))/gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	halfcnt := cnt / 2
	prime, err := num.NextPrime(start, gap)
	for i := 0; i < halfcnt; i++ {
		if err != nil {
			halfcnt = i
		} else {
			primes[i] = num.NewModulus(prime)
			prime, err = num.NextPrime(prime, gap)
		}
	}

	// Sample the other half of the primes from the smaller side.
	start = (uint64(math.Floor(math.Exp2(bits)))/gap)*gap + 1
	prime, err = num.PrevPrime(start, gap)
	currcnt := cnt
	for i := halfcnt; i < cnt; i++ {
		if err != nil {
			currcnt = i
			break
		} else {
			primes[i] = num.NewModulus(prime)
			prime, err = num.PrevPrime(prime, gap)
		}
	}

	// If the number of sampled primes is less than cnt, sample the remaining primes from the larger side.
	if currcnt < cnt {
		if currcnt != 0 {
			start = primes[currcnt-1].Value()
		}
		prime, err = num.NextPrime(start, gap)
		for i := currcnt; i < cnt; i++ {
			if err != nil {
				return nil, err
			}
			primes[i] = num.NewModulus(prime)
			prime, err = num.NextPrime(prime, gap)
		}
	}

	// Sort the primes by value.
	slices.SortFunc(primes, func(a, b *num.Modulus) int {
		return int(a.Value() - b.Value())
	})

	return primes, nil
}
