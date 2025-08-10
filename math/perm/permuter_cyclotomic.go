package perm

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

type cyclotomicPow2Permuter struct {
	params dft.RingParameters
	mod    []*num.Modulus

	buf []uint64
}

func newCyclotomicPow2Permuter(params dft.RingParameters, mod []*num.Modulus) *cyclotomicPow2Permuter {
	return &cyclotomicPow2Permuter{
		params: params,
		mod:    mod,

		buf: make([]uint64, params.Rank()),
	}
}

func (p *cyclotomicPow2Permuter) ModLen() int {
	return len(p.mod)
}

func (p *cyclotomicPow2Permuter) Permute(idx uint64, poly crt.Poly) {
	idx = idx % uint64(p.params.CycloOrder())

	if idx&1 == 0 {
		panic("Permute: idx must be odd")
	} else if poly.ModLen() != p.ModLen() {
		panic("Permute: poly.ModLen() != p.ModLen()")
	} else if idx != 0 {
		for i := 0; i < p.ModLen(); i++ {
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
