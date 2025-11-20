package bfv

import (
	"math"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

func FindNTTPrimesFromBits(params dft.RingParameters, modulusBits, auxModulusBits float64) ([]*num.Modulus, []*num.Modulus) {
	modLen := int(math.Ceil(modulusBits / num.MaxModulusBits))
	auxModLen := int(math.Ceil(auxModulusBits / num.MaxModulusBits))

	modulus := make([]*num.Modulus, modLen)
	auxModulus := make([]*num.Modulus, auxModLen)

	// Sample the modulus.
	bitlen := modulusBits
	gap, err := dft.RingGap(params)
	if err != nil {
		panic(err)
	}
	start := (uint64(math.Round(math.Exp2(modulusBits/float64(modLen))))/gap)*gap + 1
	prime := num.MustPrevPrime(start, gap)
	for i := 0; i < modLen-1; i++ {
		modulus[i] = num.NewModulus(num.MustPrevPrime(prime, gap))
		bitlen -= num.Log2(modulus[i].Value())
		prime = num.MustPrevPrime(prime, gap)
	}
	// TODO: Can bitlen be larger than 62?
	modulus[modLen-1] = num.NewModulus(num.MustPrevPrime((uint64(math.Round(math.Exp2(bitlen)))/gap)*gap+1, gap))

	if bitlen-1 > num.Log2(modulus[modLen-1].Value()) {
		panic("FindNTTPrimesFromBits: failed to sample the modulus")
	}

	slices.SortFunc(modulus, func(a, b *num.Modulus) int {
		return int(a.Value() - b.Value())
	})

	// Sample the auxiliary modulus.
	bitlen = auxModulusBits
	start = (uint64(math.Round(math.Exp2(auxModulusBits/float64(auxModLen))))/gap)*gap + 1
	cnt := 0
	prime = num.MustPrevPrime(start, gap)
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

	auxModulus[cnt] = num.NewModulus(num.MustPrevPrime((uint64(math.Round(math.Exp2(bitlen)))/gap)*gap+1, gap))

	if bitlen-1 > num.Log2(auxModulus[cnt].Value()) {
		panic("FindNTTPrimesFromBits: failed to sample the auxiliary modulus")
	}

	slices.SortFunc(auxModulus, func(a, b *num.Modulus) int {
		return int(a.Value() - b.Value())
	})

	return modulus, auxModulus
}
