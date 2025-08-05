package dft

import (
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclicTwiddleFactor computes the (inverse) twiddle factor for Cyclic (Inv)NTT.
func cyclicTwiddleFactor(rank, radix int, root []uint64, mod *num.Modulus) (tw, twInv []uint64) {
	if rank == 1 {
		return []uint64{1}, []uint64{1}
	}

	z := num.NthRoot(rank, root, mod)
	zInv := num.Inv(z, mod)

	tw = make([]uint64, rank)
	twInv = make([]uint64, rank)
	tw[0] = 1
	twInv[0] = 1
	for i := 1; i < rank/radix; i++ {
		tw[i] = num.Mul(tw[i-1], z, mod)
		twInv[i] = num.Mul(twInv[i-1], zInv, mod)
	}

	vec.RadixReverseInPlace(tw[:rank/radix], radix)
	vec.RadixReverseInPlace(twInv[:rank/radix], radix)

	t := rank / radix
	for i := 0; i < t; i++ {
		tw[i+t] = tw[i]
		twInv[i+t] = twInv[i]
		for j := 2 * t; j < radix*t; j += t {
			tw[i+j] = num.Mul(tw[i+j-t], tw[i], mod)
			twInv[i+j] = num.Mul(twInv[i+j-t], twInv[i], mod)
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
				tw[j*t+i] = num.Mul(tw[t*(j-1)+i], tw[t+i], mod)
				twInv[j*t+i] = num.Mul(twInv[t*(j-1)+i], twInv[t+i], mod)
			}
		}
	}

	return tw, twInv
}
