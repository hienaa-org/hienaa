package dftops

import (
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// CyclicPow2Transformer is a transformer for power-of-two ranks.
// Internally used as ambient space.
type CyclicPow2Transformer struct {
	rank int
	mod  *num.Modulus

	// tw is the twiddle factor for NTT.
	tw []uint64
	// twS is the Shoup form of tw.
	twS []uint64
	// twInv is the twiddle factor for InvNTT.
	twInv []uint64
	// twInvS is the Shoup form of twInv.
	twInvS []uint64

	// rankInv is the modular inverse of the rank.
	rankInv uint64
}

// NewCyclicPow2Transformer creates a new [CyclicPow2Transformer].
func NewCyclicPow2Transformer(rank int, mod *num.Modulus) *CyclicPow2Transformer {
	root := num.Generators(mod)

	tw := make([]uint64, rank/2)
	twInv := make([]uint64, rank/2)
	tw[0], tw[1] = 1, num.NthRoot(rank, root, mod)
	twInv[0], twInv[1] = 1, num.Inv(tw[1], mod)
	for i := 2; i < rank/2; i++ {
		tw[i] = num.Mul(tw[i-1], tw[1], mod)
		twInv[i] = num.Mul(twInv[i-1], twInv[1], mod)
	}
	vec.RadixReverseInPlace(tw, 2)
	vec.RadixReverseInPlace(twInv, 2)

	twS := make([]uint64, rank/2)
	twInvS := make([]uint64, rank/2)
	for i := 0; i < rank/2; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	rankInv := num.InvMForm(num.Inv(uint64(rank), mod), mod)

	return &CyclicPow2Transformer{
		rank: rank,
		mod:  mod,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: rankInv,
	}
}

func (ntt *CyclicPow2Transformer) ForwardInPlace(coeffs []uint64) {
	NTTInPlacePow2(coeffs, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(coeffs, coeffs, ntt.mod)
}

func (ntt *CyclicPow2Transformer) InverseInPlace(coeffs []uint64) {
	INTTInPlacePow2(coeffs, ntt.twInv, ntt.twInvS, ntt.mod.Value())
	vec.ScalarMulTo(coeffs, coeffs, ntt.rankInv, ntt.mod)
}
