package dft

import (
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclotomicPow2Transformer is a transformer for power-of-two cyclotomic NTT.
type cyclotomicPow2Transformer struct {
	params RingParameters
	mod    *num.Modulus

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

// newCyclotomicPow2Transformer creates a new [pow2CyclotomicTransformer].
func newCyclotomicPow2Transformer(ringParams RingParameters, mod *num.Modulus) *cyclotomicPow2Transformer {
	root := num.Generators(mod)

	tw := make([]uint64, ringParams.rank)
	twInv := make([]uint64, ringParams.rank)
	tw[0], tw[1] = 1, num.NthRoot(ringParams.cycloOrd, root, mod)
	twInv[0], twInv[1] = 1, num.Inv(tw[1], mod)
	for i := 2; i < ringParams.rank; i++ {
		tw[i] = num.Mul(tw[i-1], tw[1], mod)
		twInv[i] = num.Mul(twInv[i-1], twInv[1], mod)
	}
	vec.RadixReverseInPlace(tw, 2)
	vec.RadixReverseInPlace(twInv, 2)

	twS := make([]uint64, ringParams.rank)
	twInvS := make([]uint64, ringParams.rank)
	for i := 0; i < ringParams.rank; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	rankInv := num.InvMForm(num.Inv(uint64(ringParams.rank), mod), mod)

	return &cyclotomicPow2Transformer{
		params: ringParams,
		mod:    mod,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: rankInv,
	}
}

func (ntt *cyclotomicPow2Transformer) ForwardInPlace(coeffs []uint64) {
	nttInPlacePow2(coeffs, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(coeffs, coeffs, ntt.mod)
}

func (ntt *cyclotomicPow2Transformer) InverseInPlace(coeffs []uint64) {
	inttInPlacePow2(coeffs, ntt.twInv, ntt.twInvS, ntt.mod.Value())
	vec.ScalarMulTo(coeffs, coeffs, ntt.rankInv, ntt.mod)
}

func (ntt *cyclotomicPow2Transformer) SafeCopy() Transformer {
	return ntt
}
