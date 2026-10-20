package rlwe

import (
	"math"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// FindNTTPrimes returns distinct NTT-friendly base and auxiliary prime moduli for the given parameters.
// The logarithms of the products of baseMod and auxMod are smaller but as close as possible
// to baseModBits and auxModBits respectively.
func FindNTTPrimes(params dft.RingParameters, baseModBits, auxModBits float64) (baseMod []*num.Modulus, auxMod []*num.Modulus) {
	found := make(map[uint64]struct{})
	return prevPrimes(params, baseModBits, found), prevPrimes(params, auxModBits, found)
}

func prevPrimes(params dft.RingParameters, totalBits float64, found map[uint64]struct{}) []*num.Modulus {
	if totalBits == 0 {
		return nil
	}

	gap := dft.NTTPrimeGap(params)

	cnt := max(1, int(math.Ceil(totalBits/num.MaxModulusBits)))
	r := max(2, uint64(math.Floor(math.Exp2(totalBits/float64(cnt)))))
	start := ((r-1)/gap)*gap + 1

	primes := make([]*num.Modulus, 0)
	for len(primes) < cnt {
		p := num.PrevPrime(start, gap)
		start = p - gap
		if start > p {
			panic("no previous prime exists")
		}
		if _, ok := found[p]; ok {
			continue
		}

		primes = append(primes, num.NewModulus(p))
		found[p] = struct{}{}
	}

	return primes
}
