package dft

import (
	"errors"
	"math"
	"slices"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
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
	// cycloIdx is the index of the underlying cyclotomic polynomial.
	// 0 if the RingType is not [TypeCyclotomic] or [TypeAutFixed].
	cycloIdx int
	// rank is the number of coefficients of the polynomial in the ring.
	rank int
	// expFac is the expansion factor of the ring.
	// If the expansion factor is too large, it is set to 0.
	expFac uint64
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

	if len(cycloPoly)-1 <= 1 {
		panic("rank must be larger than 1")
	}

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
	if rank <= 1 {
		panic("rank must be larger than 1")
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
	if rank <= 1 {
		panic("rank must be larger than 1")
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
	} else if modPoly[len(modPoly)-1] != 1 {
		panic("modPoly must be monic")
	}

	return RingParameters{
		cycloIdx: 0,
		rank:     len(modPoly) - 1,
		expFac:   otherExpFac(modPoly),
		ringType: TypeOther,
		modPoly:  modPoly,
	}
}

// CycloIndex is the index of the underlying cyclotomic polynomial.
// 0 if the RingType is not [TypeCyclotomic] or [TypeAutFixed].
func (p RingParameters) CycloIndex() int {
	return p.cycloIdx
}

// Rank is the number of coefficients of the polynomial in the ring.
func (p RingParameters) Rank() int {
	return p.rank
}

// ExpandFactor is the expansion factor of the ring.
// If the expansion factor is too large, it is set to 0.
func (p RingParameters) ExpandFactor() uint64 {
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
func cyclotomicExpFac(cycloIdx int, cycloPoly []int64) uint64 {
	if num.IsPowerOfTwo(cycloIdx) {
		return uint64(cycloIdx >> 1)
	}

	rank := len(cycloPoly) - 1
	primes, _ := num.Factor(cycloIdx)
	if len(primes) == 1 {
		return uint64(2*rank - cycloIdx/primes[0])
	}

	return otherExpFac(cycloPoly)
}

// cyclicExpFac computes the expansion factor of the cyclic ring.
func cyclicExpFac(rank int) uint64 {
	return uint64(rank)
}

// autFixedExpFac computes the expansion factor of the autfixed ring.
func autFixedExpFac(cycloIdx int) uint64 {
	if num.IsPowerOfTwo(cycloIdx) {
		return uint64(cycloIdx >> 1)
	}

	return uint64(2*cycloIdx - 2)
}

// otherExpFac computes the expansion factor of arbitrary quotient ring.
func otherExpFac(modPoly []int64) uint64 {
	deg := len(modPoly) - 1
	modPolyFloat := vec.Cast[float64](modPoly)

	bounds := make([]float64, deg)
	r := make([]float64, deg)
	r[0] = 1

	for m := 0; m < 2*deg-1; m++ {
		w := float64(min(m+1, 2*deg-1-m))
		for k := 0; k < deg; k++ {
			if r[k] == 0 {
				continue
			}
			bounds[k] += w * math.Abs(r[k])
		}

		if m == 2*deg-2 {
			break
		}

		lc := r[len(r)-1]
		copy(r[1:], r[:len(r)-1])
		r[0] = 0
		if lc == 0 {
			continue
		}
		for i := 0; i < deg; i++ {
			if modPolyFloat[i] == 0 {
				continue
			}
			r[i] -= lc * modPolyFloat[i]
		}
	}

	for i := range bounds {
		if math.IsNaN(bounds[i]) || bounds[i] >= math.Ldexp(1, 52) {
			return 0
		}
	}

	return uint64(math.Ceil(vec.Max(bounds)))
}

// cyclotomicGap finds the "gap" of the NTT-friendly modulus for cyclotomic rings.
func cyclotomicGap(cycloIdx, rank int) uint64 {
	if num.IsPowerOfTwo(cycloIdx) {
		return uint64(cycloIdx)
	}

	gap := cyclicGap(cycloIdx)
	primes, _ := num.Factor(cycloIdx)
	var redDeg int
	if cycloIdx%2 == 1 {
		redDeg = cycloIdx - cycloIdx/primes[0]
	} else {
		redDeg = cycloIdx/2 - (cycloIdx/2)/primes[1]
	}

	if redDeg != rank {
		degNext := num.NextProdPower((rank)+1, []int{2})
		diffDegNext := num.NextProdPower((2*(redDeg-rank)+1)+1, []int{2})
		gap = num.LCM(gap, num.LCM(uint64(degNext), uint64(diffDegNext)))
	}

	return gap
}

// cyclicGap finds the "gap" of the NTT-friendly modulus for cyclic rings.
func cyclicGap(rank int) uint64 {
	var gap uint64

	if num.IsProdPowerOf(rank, cyclicNTTFactors) {
		gap = uint64(rank)
	} else {
		gap = num.LCM(uint64(num.NextProdPower(2*rank-1, []int{2})), uint64(rank))
	}

	return gap
}

// autFixedGap finds the "gap" of the NTT-friendly modulus for autfixed rings.
func autFixedGap(cycloIdx, rank int) uint64 {
	var gap uint64

	if num.IsPowerOfTwo(cycloIdx) {
		gap = uint64(cycloIdx)
	} else if num.IsProdPowerOf(rank, []int{2}) {
		gap = uint64(cycloIdx) * uint64(rank)
	} else {
		gap = uint64(cycloIdx) * uint64(num.NextProdPower(2*rank-1, []int{2}))
	}

	return gap
}

// nttPrimeGap returns the "gap" of the NTT-friendly modulus for the given ring parameters.
func nttPrimeGap(params RingParameters) (uint64, error) {
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

// NTTPrimeGap returns the "gap" of the NTT-friendly modulus for the given ring parameters.
//
// Panics when params is [TypeOther].
func NTTPrimeGap(params RingParameters) uint64 {
	gap, err := nttPrimeGap(params)
	if err != nil {
		panic(err)
	}
	return gap
}

// IsNTTFriendly checks if the given modulus is NTT-friendly with respect to the ring parameters.
func IsNTTFriendly(params RingParameters, mod *num.Modulus) bool {
	gap, err := nttPrimeGap(params)
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

// NextNTTPrimes finds a list of prime moduli that are NTT-friendly and greater than or equal to 2^bits.
//
// Panics if no such primes exist.
func NextNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	gap := NTTPrimeGap(params)

	r := max(uint64(math.Ceil(math.Exp2(bits))), 2)
	start := num.DivCeil(r-1, gap)*gap + 1
	primes := make([]*num.Modulus, cnt)
	for i := range primes {
		primes[i] = num.NewModulus(num.NextPrime(start, gap))
		start = primes[i].Value() + gap
		if start < primes[i].Value() {
			panic("no next prime exists")
		}
	}
	return primes
}

// PrevNTTPrimes finds a list of prime moduli that are NTT-friendly and greater than or equal to 2^bits.
//
// Panics if no such primes exist.
func PrevNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	gap := NTTPrimeGap(params)

	r := max(uint64(math.Floor(math.Exp2(bits))), 2)
	start := ((r-1)/gap)*gap + 1

	primes := make([]*num.Modulus, cnt)
	for i := range primes {
		primes[i] = num.NewModulus(num.PrevPrime(start, gap))
		start = primes[i].Value() - gap
		if start > primes[i].Value() {
			panic("no previous prime exists")
		}
	}
	return primes
}

// NearestNTTPrimes finds a list of prime moduli that are NTT-friendly and alternates from 2^bits.
//
// Panics if no such primes exist.
func NearestNTTPrimes(params RingParameters, bits float64, cnt int) []*num.Modulus {
	gap := NTTPrimeGap(params)

	r := uint64(math.Floor(math.Exp2(bits)))
	nextStart := num.DivCeil(r, gap)*gap + 1
	prevStart := nextStart - gap
	if cnt > 1 && prevStart > nextStart {
		panic("no previous prime exists")
	}

	primes := make([]*num.Modulus, cnt)
	for i := 0; i < cnt; i += 2 {
		primes[i] = num.NewModulus(num.NextPrime(nextStart, gap))
		nextStart = primes[i].Value() + gap
		if nextStart < primes[i].Value() {
			panic("no next prime exists")
		}
	}
	for i := 1; i < cnt; i += 2 {
		primes[i] = num.NewModulus(num.PrevPrime(prevStart, gap))
		prevStart = primes[i].Value() - gap
		if prevStart > primes[i].Value() {
			panic("no previous prime exists")
		}
	}
	return primes
}
