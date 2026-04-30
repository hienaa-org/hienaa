package rlwe

import (
	"cmp"
	"math"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// FindNTTPrimes returns NTT primes for the given parameters,
// where the product of modulus and auxModulus is approximately 2^modulusBits and 2^auxModulusBits respectively.
func FindNTTPrimes(params dft.RingParameters, baseModBits, auxModBits float64) ([]*num.Modulus, []*num.Modulus) {
	modLen := int(math.Ceil(baseModBits / num.MaxModulusBits))
	auxModLen := int(math.Ceil(auxModBits / num.MaxModulusBits))

	var modulus, auxModulus []*num.Modulus
	var bitlen float64

	for {
		modulus = make([]*num.Modulus, modLen)

		gap, err := dft.NTTPrimeGap(params)
		if err != nil {
			panic(err)
		}

		limbBit := math.Ceil(baseModBits / float64(modLen))
		if limbBit >= 62 {
			limbBit = 61
		}

		start := (uint64(math.Round(math.Exp2(limbBit)))/gap)*gap + 1
		prime, err := num.PrevPrime(start, gap)
		if err != nil {
			panic(err)
		}
		bitlen = baseModBits
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
			slices.SortFunc(modulus[:modLen-1], func(a, b *num.Modulus) int {
				return cmp.Compare(a.Value(), b.Value())
			})

			start := (uint64(math.Round(math.Exp2(bitlen)))/gap)*gap + 1
			prime := num.MustPrevPrime(start, gap)
			primemod := num.NewModulus(prime)
			for {
				_, ok := slices.BinarySearchFunc(modulus[:modLen-1], primemod, func(a, b *num.Modulus) int {
					return cmp.Compare(a.Value(), b.Value())
				})

				if !ok {
					break
				} else {
					prime = num.MustPrevPrime(prime, gap)
					primemod = num.NewModulus(prime)
				}
			}
			modulus[modLen-1] = primemod

			if bitlen-1 > num.Log2(modulus[modLen-1].Value()) {
				panic("FindNTTPrimesFromBits: failed to sample the modulus")
			}

			break
		}
	}

	slices.SortFunc(modulus, func(a, b *num.Modulus) int {
		return cmp.Compare(a.Value(), b.Value())
	})

	// Sample the auxiliary modulus.
	if auxModLen != 0 {
		for {
			auxModulus = make([]*num.Modulus, auxModLen)

			gap, err := dft.NTTPrimeGap(params)
			if err != nil {
				panic(err)
			}

			bitlen = auxModBits
			start := (uint64(math.Round(math.Exp2(auxModBits/float64(auxModLen))))/gap)*gap + 1
			prime := num.MustPrevPrime(start, gap)

			cnt := 0
			for cnt < auxModLen-1 {
				primemod := num.NewModulus(prime)
				for {
					_, ok := slices.BinarySearchFunc(modulus, primemod, func(a, b *num.Modulus) int {
						return cmp.Compare(a.Value(), b.Value())
					})

					if !ok {
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
				slices.SortFunc(auxModulus[:auxModLen-1], func(a, b *num.Modulus) int {
					return cmp.Compare(a.Value(), b.Value())
				})

				start := (uint64(math.Round(math.Exp2(bitlen)))/gap)*gap + 1
				prime := num.MustPrevPrime(start, gap)
				primemod := num.NewModulus(prime)
				for {
					_, okMod := slices.BinarySearchFunc(modulus, primemod, func(a, b *num.Modulus) int {
						return cmp.Compare(a.Value(), b.Value())
					})

					_, okAux := slices.BinarySearchFunc(auxModulus[:auxModLen-1], primemod, func(a, b *num.Modulus) int {
						return cmp.Compare(a.Value(), b.Value())
					})

					if !(okMod || okAux) {
						break
					} else {
						prime = num.MustPrevPrime(prime, gap)
						primemod = num.NewModulus(prime)
					}
				}
				auxModulus[auxModLen-1] = primemod

				if bitlen-1 > num.Log2(auxModulus[auxModLen-1].Value()) {
					panic("FindNTTPrimesFromBits: failed to sample the auxiliary modulus")
				}

				break
			}
		}

		slices.SortFunc(auxModulus, func(a, b *num.Modulus) int {
			return int(a.Value() - b.Value())
		})
	}

	return modulus, auxModulus
}
