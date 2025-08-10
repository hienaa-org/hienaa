package perm

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

type autFixedPow2Permuter struct {
	params dft.RingParameters
	mod    []*num.Modulus

	buf []uint64
}

func newAutFixedPow2Permuter(params dft.RingParameters, mod []*num.Modulus) *autFixedPow2Permuter {
	return &autFixedPow2Permuter{
		params: params,
		mod:    mod,

		buf: make([]uint64, params.Rank()),
	}
}

func (p *autFixedPow2Permuter) ModLen() int {
	return len(p.mod)
}

func (p *autFixedPow2Permuter) Permute(idx uint64, poly crt.Poly) {
	idx = idx % uint64(p.params.CycloOrder())

	if idx%4 != 1 {
		panic("Permute: idx must be 1 mod 4")
	} else if poly.ModLen() != p.ModLen() {
		panic("Permute: poly.ModLen() != p.ModLen()")
	} else if idx != 0 {
		for i := 0; i < p.ModLen(); i++ {
			if poly.IsNTT() && dft.IsNTTFriendly(p.params, p.mod[i]) {
				copy(p.buf, poly.Coeffs[i])

				mask := p.params.CycloOrder() - 1
				rank := p.params.Rank()
				revShiftBits := 64 - (num.Log2(uint64(p.params.Rank())) + 1)

				for j := 0; j < rank; j++ {
					oldIdx := int(bits.Reverse64(uint64(((j<<1+1)*int(idx)&mask-1)>>1)) >> revShiftBits)
					if oldIdx >= rank {
						oldIdx = rank<<1 - oldIdx - 1
					}
					newIdx := int(bits.Reverse64(uint64(j)) >> revShiftBits)
					if newIdx >= rank {
						newIdx = rank<<1 - newIdx - 1
					}

					poly.Coeffs[i][newIdx] = p.buf[oldIdx]
				}
			} else {
				clear(p.buf)
				mask := p.params.CycloOrder() - 1
				rank := p.params.Rank()

				for j := 0; j < rank; j++ {
					largeIdx := (j * int(idx)) & mask
					if largeIdx >= 3*rank {
						p.buf[rank<<2-largeIdx] = poly.Coeffs[i][j]
					} else if largeIdx >= rank<<1 {
						p.buf[largeIdx-rank<<1] = num.Neg(poly.Coeffs[i][j], p.mod[i])
					} else if largeIdx >= p.params.Rank() {
						p.buf[rank<<1-largeIdx] = num.Neg(poly.Coeffs[i][j], p.mod[i])
					} else {
						p.buf[largeIdx] = poly.Coeffs[i][j]
					}
				}

				copy(poly.Coeffs[i], p.buf)
			}
		}
	}
}

type autFixedPrimePermuter struct {
	params dft.RingParameters
	mod    []*num.Modulus

	rootPow    []uint64
	rootPowInv []uint64

	buf []uint64
}

func newAutFixedPrimePermuter(params dft.RingParameters, mod []*num.Modulus) *autFixedPrimePermuter {
	cycloOrdMod := num.NewModulus(uint64(params.CycloOrder()))
	root := num.Generators(cycloOrdMod)[0]
	rootInv := num.Inv(root, cycloOrdMod)

	rootPow := make([]uint64, params.Rank())
	rootPowInv := make([]uint64, params.Rank())

	rootPow[0] = 1
	rootPowInv[0] = 1

	for i := 1; i < params.Rank(); i++ {
		rootPow[i] = num.Mul(rootPow[i-1], root, cycloOrdMod)
		rootPowInv[i] = num.Mul(rootPowInv[i-1], rootInv, cycloOrdMod)
	}

	return &autFixedPrimePermuter{
		params: params,
		mod:    mod,

		rootPow:    rootPow,
		rootPowInv: rootPowInv,

		buf: make([]uint64, params.Rank()),
	}
}

func (p *autFixedPrimePermuter) ModLen() int {
	return len(p.mod)
}

func (p *autFixedPrimePermuter) Permute(idx uint64, poly crt.Poly) {
	if poly.ModLen() != p.ModLen() {
		panic("Permute: poly.ModLen() != p.ModLen()")
	}

	idx = idx % uint64(p.params.CycloOrder())
	if idx != 0 {
		rotIdx, isFound := slices.BinarySearch(p.rootPow, idx)

		if !isFound {
			rotIdx, isFound = slices.BinarySearch(p.rootPowInv, idx)

			if !isFound {
				panic("Permute: idx is not a valid index")
			}

			rotIdx = (p.params.Rank() - rotIdx) % p.params.Rank()
		}

		for i := 0; i < p.ModLen(); i++ {
			if poly.IsNTT() && dft.IsNTTFriendly(p.params, p.mod[i]) {
				copy(p.buf[:p.params.Rank()-rotIdx], poly.Coeffs[i][rotIdx:])
				copy(p.buf[p.params.Rank()-rotIdx:], poly.Coeffs[i][:rotIdx])
				copy(poly.Coeffs[i], p.buf)
			} else {
				copy(p.buf[:rotIdx], poly.Coeffs[i][p.params.Rank()-rotIdx:])
				copy(p.buf[rotIdx:], poly.Coeffs[i][:p.params.Rank()-rotIdx])
				copy(poly.Coeffs[i], p.buf)
			}
		}
	}
}
