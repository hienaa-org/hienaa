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
	// cycloOrder is the order of the underlying cyclotomic polynomial.
	// 0 if the RingType is not [Cyclotomic] or [AutFixed].
	cycloOrd int

	// rank is the number of coefficients of the polynomial in the ring.
	rank int

	// expFactor is the expansion factor of the ring.
	expFactor int

	// ringType is the type of the ring.
	ringType RingType
}

// NewCyclotomicParameters creates a new [RingParameters] for a cyclotomic ring.
func NewCyclotomicParameters(cycloOrd int) RingParameters {
	if cycloOrd <= 0 {
		panic("cycloOrd must be positive")
	}

	return RingParameters{
		cycloOrd:  cycloOrd,
		rank:      int(num.Totient(uint64(cycloOrd))),
		expFactor: cyclotomicExpFactor(cycloOrd),
		ringType:  TypeCyclotomic,
	}
}

// NewCyclicParameters creates a new [RingParameters] for a cyclic ring.
func NewCyclicParameters(rank int) RingParameters {
	if rank <= 0 {
		panic("rank must be positive")
	}

	return RingParameters{
		cycloOrd:  0,
		rank:      rank,
		expFactor: cyclicExpFactor(rank),
		ringType:  TypeCyclic,
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
		cycloOrd:  cycloOrd,
		rank:      rank,
		expFactor: autFixedExpFactor(cycloOrd, rank),
		ringType:  TypeAutFixed,
	}
}

// NewOtherParameters creates a new [RingParameters] for arbitrary quotient ring.
func NewOtherParameters(modPoly []int64) RingParameters {
	if len(modPoly) == 0 {
		panic("modPoly must be non-empty")
	}

	return RingParameters{
		cycloOrd:  0,
		rank:      len(modPoly) - 1,
		expFactor: otherExpFactor(modPoly),
		ringType:  TypeOther,
	}
}

// CycloOrder is the order of the underlying cyclotomic polynomial.
// 0 if the RingType is [TypeCyclic].
func (p RingParameters) CycloOrder() int {
	return p.cycloOrd
}

// Rank is the number of coefficients of the polynomial in the ring.
func (p RingParameters) Rank() int {
	return p.rank
}

// ExpFactor is the expansion factor of the ring.
func (p RingParameters) ExpFactor() int {
	return p.expFactor
}

// RingType is the type of the ring.
func (p RingParameters) RingType() RingType {
	return p.ringType
}

// ExpFactor is the expansion factor of the ring.
func cyclotomicExpFactor(cycloOrd int) int {
	if num.IsPowerOfTwo(cycloOrd) {
		return cycloOrd >> 1
	} else {
		primes, exps := num.Factor(cycloOrd)
		if len(primes) == 1 {
			return 2 * int(num.Exp(uint64(primes[0]), uint64(exps[0]-1), nil)) * int(primes[0]-1)
		} else {
			cycloPoly := CyclotomicPolynomial(cycloOrd)
			tot := num.Totient(cycloOrd)
			max := 0
			for i := 0; i < 2*tot-cycloOrd; i++ {
				sum := 0
				for j := 0; j < cycloOrd-tot; j++ {
					sum += int(math.Abs(float64(cycloPoly[j+i])))
				}
				if sum > max {
					max = sum
				}
			}

			return (max + 1) * tot
		}
	}
}

// cyclicExpFactor is the expansion factor of the cyclic ring.
func cyclicExpFactor(rank int) int {
	return rank
}

// autFixedExpFactor is the expansion factor of the autfixed ring.
func autFixedExpFactor(cycloOrd, rank int) int {
	// TODO: Can we reduce the expansion factor with respect to the rank?
	if num.IsPowerOfTwo(cycloOrd) {
		return cycloOrd >> 1
	} else {
		return 2*cycloOrd - 2
	}
}

// otherExpFactor is the expansion factor of the arbitrary quotient ring.
func otherExpFactor(modPoly []int64) int {
	len := len(modPoly)
	cBound := make([]int, 2*len-1)
	for i := 0; i < len; i++ {
		cBound[i] = i + 1
		cBound[2*len-2-i] = i + 1
	}

	for i := len - 1; i >= 0; i-- {
		for j := 0; j < len; j++ {
			cBound[i+j] += int(math.Abs(float64(modPoly[j]))) * cBound[i+len-1]
		}
		cBound[i+len-1] = 0
	}

	max := 0
	for i := 0; i < len; i++ {
		if cBound[i] > max {
			max = cBound[i]
		}
	}

	return max
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

// NTTPrimeGap returns the "gap" of the NTT-friendly modulus for the given ring parameters.
func NTTPrimeGap(params RingParameters) (uint64, error) {
	switch params.ringType {
	case TypeCyclic:
		return cyclicGap(params.rank), nil
	case TypeCyclotomic:
		return cyclotomicGap(params.cycloOrd, params.rank), nil
	case TypeAutFixed:
		return autFixedGap(params.cycloOrd, params.rank), nil
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
		if f%gap != 1 {
			return false
		}
	}
	return true
}

// MustFindNextNTTPrimes finds a list of prime moduli that are NTT-friendly with respect to the given ring parameters.
// Specifically, it outputs the first cnt NTT-friendly primes greater than or equal to 2^bits.
// It panics if an error occurs.
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
// Specifically, it outputs the first cnt NTT-friendly primes less than or equal to 2^bits.
// It panics if an error occurs.
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
			return primes, nil
		}
	}
}
