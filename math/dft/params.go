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
	// TypeCyclotomic is a cyclotomic ring ZZ[X]/Phi_M(X).
	TypeCyclotomic RingType = iota
	// TypeCyclic is a cyclic ring ZZ[X]/(X^N - 1).
	TypeCyclic
	// TypeAutFixed is a decomposition ring of a cyclotomic ring.
	// In other words, it is a subring of a cyclotomic ring invariant under some automorphism.
	TypeAutFixed
	// TypeOther covers arbitrary quotient rings.
	TypeOther
)

// RingParameters contains the parameters for the ring.
type RingParameters struct {
	// cycloIdx is the order of the underlying cyclotomic polynomial.
	// 0 if the RingType is not [Cyclotomic] or [AutFixed].
	cycloIdx int
	// rank is the number of coefficients of the polynomial in the ring.
	rank int
	// expFac is the expansion factor of the ring.
	expFac int
	// ringType is the type of the ring.
	ringType RingType
	// modPoly is the modulus polynomial of the ring.
	modPoly []int64
}

// NewCyclotomicParameters creates a new [RingParameters] for a cyclotomic ring.
func NewCyclotomicParameters(cycloIdx int) RingParameters {
	if cycloIdx <= 0 {
		panic("cycloIdx must be positive")
	}

	cycloPoly := CyclotomicPolynomial(cycloIdx)

	return RingParameters{
		cycloIdx: cycloIdx,
		rank:     len(cycloPoly) - 1,
		expFac:   cyclotomicExpFac(cycloIdx, cycloPoly),
		ringType: TypeCyclotomic,
		modPoly:  cycloPoly,
	}
}

// NewCyclicParameters creates a new [RingParameters] for a cyclic ring.
func NewCyclicParameters(rank int) RingParameters {
	if rank <= 0 {
		panic("rank must be positive")
	}

	modPoly := make([]int64, rank+1)
	modPoly[0] = -1
	modPoly[rank] = 1

	return RingParameters{
		cycloIdx: 0,
		rank:     rank,
		expFac:   cyclicExpFac(rank),
		ringType: TypeCyclic,
		modPoly:  modPoly,
	}
}

// NewAutFixedParameters creates a new [RingParameters] for an autfixed ring.
func NewAutFixedParameters(cycloIdx, rank int) RingParameters {
	if rank <= 0 {
		panic("rank must be positive")
	}

	switch {
	case num.IsPowerOfTwo(cycloIdx):
		if cycloIdx != 4*rank {
			panic("cycloIdx must be four times the rank for power-of-two cycloIdx")
		}
	case num.IsPrime(cycloIdx):
		if (cycloIdx-1)%rank != 0 {
			panic("rank should divide cycloIdx-1 for prime cycloIdx")
		}
	default:
		panic("cycloIdx must be a prime or a power of two")
	}

	return RingParameters{
		cycloIdx: cycloIdx,
		rank:     rank,
		expFac:   autFixedExpFac(cycloIdx),
		ringType: TypeAutFixed,
		modPoly:  CyclotomicPolynomial(cycloIdx),
	}
}

// NewOtherParameters creates a new [RingParameters] for arbitrary quotient ring.
func NewOtherParameters(modPoly []int64) RingParameters {
	if len(modPoly) <= 1 {
		panic("rank must be larger than 1")
	}

	return RingParameters{
		cycloIdx: 0,
		rank:     len(modPoly) - 1,
		expFac:   otherExpFac(modPoly),
		ringType: TypeOther,
		modPoly:  modPoly,
	}
}

// CycloIndex is the order of the underlying cyclotomic polynomial.
// 0 if the RingType is [TypeCyclic].
func (p RingParameters) CycloIndex() int {
	return p.cycloIdx
}

// Rank is the number of coefficients of the polynomial in the ring.
func (p RingParameters) Rank() int {
	return p.rank
}

// ExpandFactor is the expansion factor of the ring.
func (p RingParameters) ExpandFactor() int {
	return p.expFac
}

// ModulusPoly returns the modulus polynomial of the ring.
func (p RingParameters) ModulusPoly() []int64 {
	return p.modPoly
}

// RingType is the type of the ring.
func (p RingParameters) RingType() RingType {
	return p.ringType
}

// Equal checks if two parameters are equal.
func (p RingParameters) Equal(p0 RingParameters) bool {
	eq := p.cycloIdx == p0.cycloIdx && p.rank == p0.rank && p.ringType == p0.ringType

	if p.ringType == TypeOther {
		return eq && slices.Equal(p.modPoly, p0.modPoly)
	}
	return eq
}

// cyclotomicExpFac computes the expansion factor for cyclotomic ring.
func cyclotomicExpFac(cycloIdx int, cycloPoly []int64) int {
	if num.IsPowerOfTwo(cycloIdx) {
		if cycloIdx == 1 {
			return 1
		}
		return cycloIdx >> 1
	}

	rank := len(cycloPoly) - 1
	primes, _ := num.Factor(cycloIdx)
	if len(primes) == 1 {
		return 2 * rank
	}

	var max int
	for i := 0; i < 2*rank-cycloIdx; i++ {
		var sum int
		for j := 0; j < cycloIdx-rank; j++ {
			sum += num.Abs(int(cycloPoly[j+i]))
		}
		if sum > max {
			max = sum
		}
	}

	return (max + 1) * rank
}

// cyclicExpFac computes the expansion factor of the cyclic ring.
func cyclicExpFac(rank int) int {
	return rank
}

// autFixedExpFac computes the expansion factor of the autfixed ring.
func autFixedExpFac(cycloIdx int) int {
	// TODO: Can we reduce the factor with respect to the rank?
	if num.IsPowerOfTwo(cycloIdx) {
		return cycloIdx >> 1
	}

	return 2*cycloIdx - 2
}

// otherExpFac is the expansion factor of the arbitrary quotient ring.
func otherExpFac(modPoly []int64) int {
	deg := len(modPoly)
	cBound := make([]int, 2*deg-1)
	for i := 0; i < deg; i++ {
		cBound[i] = i + 1
		cBound[2*deg-2-i] = i + 1
	}

	for i := deg - 1; i >= 0; i-- {
		for j := 0; j < deg; j++ {
			cBound[i+j] += int(num.Abs(modPoly[j])) * cBound[i+deg-1]
		}
		cBound[i+deg-1] = 0
	}

	max := 0
	for i := 0; i < deg; i++ {
		if cBound[i] > max {
			max = cBound[i]
		}
	}

	return max
}

// cyclotomicGap finds the "gap" of the NTT-friendly modulus for cyclotomic rings.
func cyclotomicGap(cycloIdx, rank int) uint64 {
	var gap int

	if num.IsPowerOfTwo(cycloIdx) {
		gap = cycloIdx
	} else {
		bluesteinRank := num.NextProdPower(2*cycloIdx-1, []int{2})

		primes, _ := num.Factor(cycloIdx)
		var redDeg int
		if cycloIdx%2 == 1 {
			redDeg = cycloIdx - cycloIdx/primes[0]
		} else {
			redDeg = cycloIdx/2 - (cycloIdx/2)/primes[1]
		}

		if redDeg == rank {
			gap = num.LCM(cycloIdx, bluesteinRank)
		} else {
			degNext := num.NextProdPower(rank, []int{2})
			diffDegNext := num.NextProdPower(2*(redDeg-rank)+1, []int{2})
			gap = num.LCM(num.LCM(degNext, diffDegNext), num.LCM(cycloIdx, bluesteinRank))
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
func autFixedGap(cycloIdx, rank int) uint64 {
	var gap int

	if num.IsPowerOfTwo(cycloIdx) {
		gap = cycloIdx
	} else if num.IsProdPowerOf(rank, []int{2}) {
		gap = cycloIdx * rank
	} else {
		gap = cycloIdx * num.NextProdPower(2*rank-1, []int{2})
	}

	return uint64(gap)
}

// NTTPrimeGap returns the "gap" of the NTT-friendly modulus for the given ring parameters.
func NTTPrimeGap(params RingParameters) (uint64, error) {
	switch params.ringType {
	case TypeCyclic:
		return cyclicGap(params.rank), nil
	case TypeCyclotomic:
		return cyclotomicGap(params.cycloIdx, params.rank), nil
	case TypeAutFixed:
		return autFixedGap(params.cycloIdx, params.rank), nil
	default:
		return 0, errors.New("unsupported parameters")
	}
}

// IsNTTFriendly checks if the given modulus is NTT-friendly with respect to the ring parameters.
func IsNTTFriendly(params RingParameters, mod *num.Modulus) bool {
	gap, err := NTTPrimeGap(params)
	if err != nil {
		return false
	}

	primes, _ := num.Factor(mod.Value())
	for _, f := range primes {
		if (f-1)%gap != 0 {
			return false
		}
	}
	return true
}

// MustFindNextNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes greater than to 2^bits.
// It panics if an error occurs.
func MustFindNextNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	primes, err := FindNextNTTPrimes(params, bits, cnt)
	if err != nil {
		panic(err)
	}
	return primes
}

// FindNextNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes greater than to 2^bits.
func FindNextNTTPrimes(params RingParameters, bits float64, cnt int) ([]*num.Modulus, error) {
	gap, err := NTTPrimeGap(params)
	if err != nil {
		return nil, err
	}

	start := (uint64(math.Floor(math.Exp2(bits))/float64(gap)))*gap + 1

	primes := make([]*num.Modulus, cnt)
	for i := 0; i < cnt; i++ {
		prime, err := num.NextPrime(start, gap)
		if err != nil {
			return nil, err
		}
		primes[i] = num.NewModulus(prime)
		start = prime
	}

	return primes, nil
}

// MustFindPrevNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes less than to 2^bits.
// It panics if an error occurs.
func MustFindPrevNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	primes, err := FindPrevNTTPrimes(params, bits, cnt)
	if err != nil {
		panic(err)
	}
	return primes
}

// FindPrevNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes less than 2^bits.
func FindPrevNTTPrimes(params RingParameters, bits float64, cnt int) ([]*num.Modulus, error) {
	gap, err := NTTPrimeGap(params)
	if err != nil {
		return nil, err
	}

	start := (uint64(math.Floor(math.Exp2(bits))/float64(gap)))*gap + 1

	primes := make([]*num.Modulus, cnt)
	for i := 0; i < cnt; i++ {
		prime, err := num.PrevPrime(start, gap)
		if err != nil {
			return nil, err
		}
		primes[i] = num.NewModulus(prime)
		start = prime
	}

	return primes, nil
}

// FindNearestNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes nearest to 2^bits.
// Output moduli are alternating in size.
// It panics if an error occurs.
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
	gap, err := NTTPrimeGap(params)
	if err != nil {
		return nil, err
	}

	start := (uint64(math.Floor(math.Exp2(bits))/float64(gap)))*gap + 1
	nextCnt := cnt >> 1
	prevCnt := cnt - nextCnt

	nextStart := start
	nextPrimes := make([]*num.Modulus, 0, nextCnt)
	for i := 0; i < nextCnt; i++ {
		prime, err := num.NextPrime(nextStart, gap)
		if err != nil {
			break
		}
		nextPrimes = append(nextPrimes, num.NewModulus(prime))
		nextStart = prime
	}

	prevStart := start
	prevPrimes := make([]*num.Modulus, 0, prevCnt)
	for i := 0; i < prevCnt; i++ {
		prime, err := num.PrevPrime(prevStart, gap)
		if err != nil {
			break
		}
		prevPrimes = append(prevPrimes, num.NewModulus(prime))
		prevStart = prime
	}

	if len(nextPrimes) < nextCnt && len(prevPrimes) < prevCnt {
		return nil, errors.New("not enough primes found")
	} else if len(nextPrimes) < nextCnt {
		prevStart := prevPrimes[prevCnt-1].Value()
		for i := 0; i < cnt-(len(nextPrimes)+len(prevPrimes)); i++ {
			prime, err := num.PrevPrime(prevStart, gap)
			if err != nil {
				return nil, err
			}
			prevPrimes = append(prevPrimes, num.NewModulus(prime))
			prevStart = prime
		}
	} else if len(prevPrimes) < prevCnt {
		nextStart := nextPrimes[nextCnt-1].Value()
		for i := 0; i < cnt-(len(nextPrimes)+len(prevPrimes)); i++ {
			prime, err := num.NextPrime(nextStart, gap)
			if err != nil {
				return nil, err
			}
			nextPrimes = append(nextPrimes, num.NewModulus(prime))
			nextStart = prime
		}
	}

	primes := append(append(make([]*num.Modulus, 0, cnt), prevPrimes...), nextPrimes...)
	slices.SortFunc(primes, num.CmpModulus)

	return primes, nil
}

// MustFindAmbientPrimes finds the list of NTT-friendly primes such that their product is at least 2^minBits.
// It panics if an error occurs.
func MustFindAmbientPrimes(params RingParameters, minBits float64) []*num.Modulus {
	primes, err := FindAmbientPrimes(params, minBits)
	if err != nil {
		panic(err)
	}
	return primes
}

// FindAmbientPrimes finds the list of NTT-friendly primes such that their product is at least 2^minBits.
func FindAmbientPrimes(params RingParameters, minBits float64) ([]*num.Modulus, error) {
	gap, err := NTTPrimeGap(params)
	if err != nil {
		return nil, err
	}

	start := (uint64(math.Floor(num.MaxModulus/float64(gap))))*gap + 1
	if start > num.MaxModulus {
		for start > num.MaxModulus {
			start -= gap
		}
	}

	bits := 0.0
	primes := make([]*num.Modulus, 0)
	for {
		prime, err := num.NextPrime(start, gap)
		if err != nil {
			break
		}
		primes = append(primes, num.NewModulus(prime))
		start = prime

		bits += num.Log2(prime)
		if bits > minBits {
			for _, p := range primes {
				if p.Value() < 1<<num.MinAmbModulusBits {
					return nil, errors.New("FindAmbientPrimes: no suitable ambient modulus exists")
				}
			}

			return primes, nil
		}
	}

	for {
		prime, err := num.PrevPrime(start, gap)
		if err != nil {
			return nil, err
		}
		primes = append(primes, num.NewModulus(prime))
		start = prime

		bits += num.Log2(prime)
		if bits > minBits {
			for _, p := range primes {
				if p.Value() < 1<<num.MinAmbModulusBits {
					return nil, errors.New("FindAmbientPrimes: no suitable ambient modulus exists")
				}
			}

			return primes, nil
		}
	}
}
