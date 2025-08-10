package perm

import (
	"github.com/hienaa-org/hienaa/math/crt"
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

type permuter interface {
	ModLen() int
	Permute(idx uint64, poly crt.Poly)
}
