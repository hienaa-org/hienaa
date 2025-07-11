package rns

import (
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclicFactorDegree factors the degree of a native cyclic NTT.
func cyclicFactorDegree(deg int, factors []uint64) []int {
	degs := make([]int, len(factors))
	for i, f := range factors {
		degs[i] = 1
		for deg%int(f) == 0 {
			deg /= int(f)
			degs[i] *= int(f)
		}
	}
	return degs
}

// cyclicTwiddleFactor computes the (inverse) twiddle factor for Cyclic (Inv)NTT.
func cyclicTwiddleFactor(deg, radix int, root uint64, modulus *mod.Modulus) (tw, twInv []uint64) {
	if deg == 1 {
		return []uint64{1}, []uint64{1}
	}

	g := num.NthRoot(deg, root, modulus)
	gInv := mod.Inv(g, modulus)

	tw = make([]uint64, deg)
	twInv = make([]uint64, deg)
	tw[0] = 1
	twInv[0] = 1
	for i := 1; i < deg/radix; i++ {
		tw[i] = mod.Mul(tw[i-1], g, modulus)
		twInv[i] = mod.Mul(twInv[i-1], gInv, modulus)
	}

	vec.RadixReverseInPlace(tw[:deg/radix], radix)
	vec.RadixReverseInPlace(twInv[:deg/radix], radix)

	t := deg / radix
	for i := 0; i < t; i++ {
		tw[i+t] = tw[i]
		twInv[i+t] = twInv[i]
		for j := 2 * t; j < radix*t; j += t {
			tw[i+j] = mod.Mul(tw[i+j-t], tw[i], modulus)
			twInv[i+j] = mod.Mul(twInv[i+j-t], twInv[i], modulus)
		}
	}

	for t > 1 {
		t /= radix
		for i := 0; i < t; i++ {
			tw[t+i] = tw[radix*t+i]
			twInv[t+i] = twInv[radix*t+i]
		}
		for j := 2; j < radix; j++ {
			for i := 0; i < t; i++ {
				tw[j*t+i] = mod.Mul(tw[t*(j-1)+i], tw[t+i], modulus)
				twInv[j*t+i] = mod.Mul(twInv[t*(j-1)+i], twInv[t+i], modulus)
			}
		}
	}

	return tw, twInv
}
