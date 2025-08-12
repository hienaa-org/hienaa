package crt

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

type Permuter struct {
	permuter permuter
}

func NewPermuter(ringParams dft.RingParameters, mod []*num.Modulus) *Permuter {
	switch ringParams.RingType() {
	case dft.Cyclotomic:
		switch {
		case num.IsPowerOfTwo(uint64(ringParams.CycloOrder())):
			return &Permuter{
				permuter: newCyclotomicPow2Permuter(ringParams, mod),
			}
		default:
			return &Permuter{
				permuter: newCyclotomicAnyPermuter(ringParams, mod),
			}
		}
	case dft.AutFixed:
		switch {
		case num.IsPowerOfTwo(uint64(ringParams.CycloOrder())):
			return &Permuter{
				permuter: newAutFixedPow2Permuter(ringParams, mod),
			}
		case num.IsPrime(uint64(ringParams.CycloOrder())):
			return &Permuter{
				permuter: newAutFixedPrimePermuter(ringParams, mod),
			}
		}
	}

	panic("NewPermuter: unsupported ring type or parameters")
}

func (p *Permuter) Permute(idx uint64, poly Poly) {
	p.permuter.permute(idx, poly)
}

type permuter interface {
	permute(idx uint64, poly Poly)
}
