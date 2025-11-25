package rlwe

import (
	"math"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

func FindNTTPrimesFromBits(params dft.RingParameters, modulusBits, auxModulusBits float64) ([]*num.Modulus, []*num.Modulus) {
	modLen := int(math.Ceil(modulusBits / num.MaxModulusBits))
	auxModLen := int(math.Ceil(auxModulusBits / num.MaxModulusBits))

	var modulus, auxModulus []*num.Modulus
	var bitlen float64

	// Sample the modulus.
	for {
		modulus = make([]*num.Modulus, modLen)

		gap, err := dft.RingGap(params)
		if err != nil {
			panic(err)
		}

		start := (uint64(math.Round(math.Exp2(modulusBits/float64(modLen))))/gap)*gap + 1
		prime, err := num.PrevPrime(start, gap)
		if err != nil {
			panic(err)
		}
		bitlen = modulusBits
		for i := 0; i < modLen-1; i++ {
			modulus[i] = num.NewModulus(prime)
			bitlen -= num.Log2(modulus[i].Value())
			prime, err = num.PrevPrime(prime, gap)
			if err != nil {
				panic(err)
			}
		}

		if bitlen >= 62 {
			modLen++
		} else {
			modulus[modLen-1] = num.NewModulus(num.MustPrevPrime((uint64(math.Round(math.Exp2(bitlen)))/gap)*gap+1, gap))

			if bitlen-1 > num.Log2(modulus[modLen-1].Value()) {
				panic("FindNTTPrimesFromBits: failed to sample the modulus")
			}

			break
		}
	}

	slices.SortFunc(modulus, func(a, b *num.Modulus) int {
		return int(a.Value() - b.Value())
	})

	// Sample the auxiliary modulus.
	for {
		auxModulus = make([]*num.Modulus, auxModLen)

		gap, err := dft.RingGap(params)
		if err != nil {
			panic(err)
		}

		bitlen = auxModulusBits
		start := (uint64(math.Round(math.Exp2(auxModulusBits/float64(auxModLen))))/gap)*gap + 1
		prime := num.MustPrevPrime(start, gap)

		cnt := 0
		for cnt < auxModLen-1 {
			primemod := num.NewModulus(prime)
			for {
				_, found := slices.BinarySearchFunc(modulus, primemod, func(a, b *num.Modulus) int {
					return int(a.Value() - b.Value())
				})

				if !found {
					break
				} else {
					prime = num.MustPrevPrime(prime, gap)
					primemod = num.NewModulus(prime)
				}
			}
			auxModulus[cnt] = primemod
			bitlen -= num.Log2(primemod.Value())
			prime = num.MustPrevPrime(prime, gap)
			cnt++
		}

		if bitlen >= 62 {
			auxModLen++
		} else {
			auxModulus[auxModLen-1] = num.NewModulus(num.MustPrevPrime((uint64(math.Round(math.Exp2(bitlen)))/gap)*gap+1, gap))

			if bitlen-1 > num.Log2(auxModulus[auxModLen-1].Value()) {
				panic("FindNTTPrimesFromBits: failed to sample the auxiliary modulus")
			}

			break
		}
	}

	slices.SortFunc(auxModulus, func(a, b *num.Modulus) int {
		return int(a.Value() - b.Value())
	})

	return modulus, auxModulus
}
