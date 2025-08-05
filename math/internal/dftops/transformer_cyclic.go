package dftops

import (
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// CyclicPow2Transformer is a transformer for power-of-two ranks.
// Internally used as ambient space.
type CyclicPow2Transformer struct {
	Rank int
	Mod  *num.Modulus

	// Tw is the twiddle factor for NTT.
	Tw []uint64
	// TwS is the Shoup form of tw.
	TwS []uint64
	// TwInv is the twiddle factor for InvNTT.
	TwInv []uint64
	// TwInvS is the Shoup form of twInv.
	TwInvS []uint64

	// RankInv is the modular inverse of the rank.
	RankInv uint64
}

// NewCyclicPow2Transformer creates a new [CyclicPow2Transformer].
func NewCyclicPow2Transformer(rank int, mod *num.Modulus) *CyclicPow2Transformer {
	root := num.Generators(mod)

	tw := make([]uint64, rank)
	twInv := make([]uint64, rank)
	tw[0], tw[1] = 1, num.NthRoot(rank, root, mod)
	twInv[0], twInv[1] = 1, num.Inv(tw[1], mod)
	for i := 2; i < rank/2; i++ {
		tw[i] = num.Mul(tw[i-1], tw[1], mod)
		twInv[i] = num.Mul(twInv[i-1], twInv[1], mod)
	}
	vec.RadixReverseInPlace(tw[:rank/2], 2)
	vec.RadixReverseInPlace(twInv[:rank/2], 2)

	t := rank / 2
	for i := 0; i < t; i++ {
		tw[i+t], twInv[i+t] = tw[i], twInv[i]
	}
	for t > 1 {
		t /= 2
		for i := 0; i < t; i++ {
			tw[t+i], twInv[t+i] = tw[2*t+i], twInv[2*t+i]
		}
	}

	twS := make([]uint64, rank)
	twInvS := make([]uint64, rank)
	for i := 0; i < rank; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	rankInv := num.InvMForm(num.Inv(uint64(rank), mod), mod)

	return &CyclicPow2Transformer{
		Rank: rank,
		Mod:  mod,

		Tw:     tw,
		TwS:    twS,
		TwInv:  twInv,
		TwInvS: twInvS,

		RankInv: rankInv,
	}
}

func (ntt *CyclicPow2Transformer) ForwardInPlace(coeffs []uint64) {
	NTTInPlacePow2(coeffs, ntt.Tw, ntt.TwS, ntt.Mod.Value())
	vec.MFormTo(coeffs, coeffs, ntt.Mod)
}

func (ntt *CyclicPow2Transformer) InverseInPlace(coeffs []uint64) {
	INTTInPlacePow2(coeffs, ntt.TwInv, ntt.TwInvS, ntt.Mod.Value())
	vec.ScalarMulTo(coeffs, coeffs, ntt.RankInv, ntt.Mod)
}
