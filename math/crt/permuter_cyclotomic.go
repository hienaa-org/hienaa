package crt

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

type cyclotomicPow2Permuter struct {
	params dft.RingParameters
	mod    []*num.Modulus

	// buf is a buffer for the permuter.
	buf []uint64
}

func newCyclotomicPow2Permuter(params dft.RingParameters, mod []*num.Modulus) *cyclotomicPow2Permuter {
	return &cyclotomicPow2Permuter{
		params: params,
		mod:    mod,

		buf: make([]uint64, params.Rank()),
	}
}

func (p *cyclotomicPow2Permuter) modLen() int {
	return len(p.mod)
}

func (p *cyclotomicPow2Permuter) permute(idx uint64, poly Poly) {
	idx = idx % uint64(p.params.CycloOrder())

	if idx&1 == 0 {
		panic("Permute: idx must be odd")
	} else if poly.ModLen() != p.modLen() {
		panic("Permute: poly.ModLen() != p.modLen()")
	} else if idx != 0 {
		for i := 0; i < p.modLen(); i++ {
			if poly.IsNTT() && dft.IsNTTFriendly(p.params, p.mod[i]) {
				copy(p.buf, poly.Coeffs[i])

				mask := p.params.CycloOrder() - 1
				revShiftBits := 64 - num.Log2(uint64(p.params.Rank()))

				for j := 0; j < p.params.Rank(); j++ {
					oldIdx := int(bits.Reverse64(uint64(((j<<1+1)*int(idx)&mask-1)>>1)) >> revShiftBits)
					newIdx := int(bits.Reverse64(uint64(j)) >> revShiftBits)

					poly.Coeffs[i][newIdx] = p.buf[oldIdx]
				}
			} else {
				clear(p.buf)
				mask := p.params.CycloOrder() - 1
				rank := p.params.Rank()

				for j := 0; j < rank; j++ {
					largeIdx := (j * int(idx)) & mask
					if largeIdx >= rank {
						p.buf[largeIdx-rank] = num.Neg(poly.Coeffs[i][j], p.mod[i])
					} else {
						p.buf[largeIdx] = poly.Coeffs[i][j]
					}
				}

				copy(poly.Coeffs[i], p.buf)
			}
		}
	}
}

type cyclotomicAnyPermuter struct {
	params dft.RingParameters
	mod    []*num.Modulus

	// reducer is the reducer for the cyclotomic polynomial.
	reducer []reducer

	// primeExpMods is the prime power factors of the cyclotomic order.
	primeExpMods []*num.Modulus
	// rootExps is the generators modulo prime power factors of the cyclotomic order.
	rootExps [][]uint64
	// dims is the dimension of the hypercube structure.
	dims []int

	// buf is a buffer for the permuter.
	buf []uint64
}

func newCyclotomicAnyPermuter(params dft.RingParameters, mod []*num.Modulus) *cyclotomicAnyPermuter {
	reducer := make([]reducer, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(params, mod[i]) {
			reducer[i] = newCyclotomicReducerNTTModulus(params, mod[i])
		} else {
			reducer[i] = newCyclotomicReducerAnyModulus(params, mod[i])
		}
	}

	primes, exps := num.Factor(uint64(params.CycloOrder()))
	pExpMods := make([]*num.Modulus, len(primes))
	rootExps := make([][]uint64, len(primes))
	dims := make([]int, len(primes))
	for i := range rootExps {
		pExp := uint64(1)
		for j := 0; j < int(exps[i]); j++ {
			pExp *= primes[i]
		}
		pExpMods[i] = num.NewModulus(pExp)
		dims[i] = int(pExp - pExp/primes[i])
		rootExps[i] = make([]uint64, dims[i])
		root := uint64(1)
		if pExp != 2 {
			root = num.GeneratorsWithFactors(pExpMods[i], []uint64{primes[i]}, []uint64{exps[i]})[0]
		}

		if primes[i] == 2 && exps[i] > 2 {
			root = 5
		}

		rootExps[i][0] = 1
		for j := 1; j < len(rootExps[i]); j++ {
			rootExps[i][j] = num.Mul(root, rootExps[i][j-1], pExpMods[i])
		}
	}

	if primes[0] == 2 {
		if exps[0] == 1 {
			rootExps = rootExps[1:]
			pExpMods = pExpMods[1:]
			dims = dims[1:]
		} else if exps[0] > 2 {
			rootExps = append([][]uint64{rootExps[0][:len(rootExps[0])/2], {1, pExpMods[0].Value() - 1}}, rootExps[1:]...)
			pExpMods = append([]*num.Modulus{pExpMods[0], pExpMods[0]}, pExpMods[1:]...)
			dims = append([]int{dims[0] / 2, 2}, dims[1:]...)
		}
	}

	return &cyclotomicAnyPermuter{
		params: params,
		mod:    mod,

		reducer: reducer,

		primeExpMods: pExpMods,
		rootExps:     rootExps,
		dims:         dims,

		buf: make([]uint64, params.CycloOrder()),
	}
}

func (p *cyclotomicAnyPermuter) modLen() int {
	return len(p.mod)
}

func (p *cyclotomicAnyPermuter) permute(idx uint64, poly Poly) {
	idx = idx % uint64(p.params.CycloOrder())

	if num.GCD(idx, uint64(p.params.CycloOrder())) != 1 {
		panic("Permute: idx must be coprime with the cyclotomic order")
	} else if poly.ModLen() != p.modLen() {
		panic("Permute: poly.ModLen() != p.ModLen()")
	} else if idx != 0 {
		cycloOrdMod := num.NewModulus(uint64(p.params.CycloOrder()))
		for i := 0; i < p.modLen(); i++ {
			if poly.IsNTT() && dft.IsNTTFriendly(p.params, p.mod[i]) {
				idxDigits := make([]uint64, len(p.primeExpMods))
				cnt := 0

				for cnt < len(idxDigits) {
					if p.primeExpMods[cnt].Value()&7 == 0 {
						idxRed := num.Reduce(idx, p.primeExpMods[cnt])
						if idxRed&3 == 1 {
							idxDigits[cnt+1] = 0
						} else {
							idxDigits[cnt+1] = 1
							idxRed = p.primeExpMods[cnt].Value() - idxRed
						}
						idxj := slices.Index(p.rootExps[cnt], idxRed)
						idxDigits[cnt] = uint64(idxj)

						cnt += 2
					} else {
						idxRed := num.Reduce(idx, p.primeExpMods[cnt])
						idxj := slices.Index(p.rootExps[cnt], idxRed)
						idxDigits[cnt] = uint64(idxj)

						cnt++
					}
				}

				copy(p.buf, poly.Coeffs[i])
				newIdxDigits := make([]uint64, len(idxDigits))
				for j := 0; j < p.params.Rank(); j++ {
					oldIdx := j
					for k := 0; k < len(newIdxDigits); k++ {
						newIdxDigits[k] = uint64((oldIdx - int(idxDigits[k]) + p.dims[k]) % p.dims[k])
						oldIdx /= int(p.dims[k])
					}

					newIdx := int(newIdxDigits[len(newIdxDigits)-1])
					for k := len(newIdxDigits) - 2; k >= 0; k-- {
						newIdx *= p.dims[k]
						newIdx += int(newIdxDigits[k])
					}

					poly.Coeffs[i][newIdx] = p.buf[j]
				}
			} else {
				clear(p.buf)

				for j := 0; j < p.params.Rank(); j++ {
					newIdx := int(num.Mul(uint64(j), idx, cycloOrdMod))
					p.buf[newIdx] = poly.Coeffs[i][j]
				}

				p.reducer[i].reduceTo(poly.Coeffs[i], p.buf)
			}
		}
	}
}
