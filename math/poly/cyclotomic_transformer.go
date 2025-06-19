package poly

import (
	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// cyclotomicPow2Transformer is a transformer for power-of-two cyclotomic NTT.
type cyclotomicPow2Transformer struct {
	ringParams RingParameters
	mod        *mod.Modulus

	// tw is the twiddle factor for NTT.
	tw []uint64
	// twS is the Shoup form of tw.
	twS []uint64
	// twInv is the twiddle factor for InvNTT.
	twInv []uint64
	// twInvS is the Shoup form of twInv.
	twInvS []uint64

	// degInv is the modular inverse of the degree.
	degInv uint64
}

// newCyclotomicPow2Transformer creates a new [pow2CyclotomicTransformer] for the given ringParams and modulus.
func newCyclotomicPow2Transformer(ringParams RingParameters, modulus *mod.Modulus) *cyclotomicPow2Transformer {
	if (modulus.Value()-1)%uint64(ringParams.cycloDegree) != 0 {
		panic("newCyclotomicPow2Transformer: modulus not NTT-friendly")
	}

	g := num.PrimitiveRoot(modulus)

	tw := make([]uint64, ringParams.degree)
	twInv := make([]uint64, ringParams.degree)
	tw[1] = num.NthRoot(ringParams.cycloDegree, g, modulus)
	twInv[1] = mod.Inv(tw[1], modulus)
	for i := 2; i < ringParams.degree; i++ {
		tw[i] = mod.Mul(tw[i-1], tw[1], modulus)
		twInv[i] = mod.Mul(twInv[i-1], twInv[1], modulus)
	}
	vec.RadixReverseInPlace(tw, 2)
	vec.RadixReverseInPlace(twInv, 2)

	twS := make([]uint64, ringParams.degree)
	twInvS := make([]uint64, ringParams.degree)
	for i := 0; i < ringParams.degree; i++ {
		twS[i] = mod.SForm(tw[i], modulus)
		twInvS[i] = mod.SForm(twInv[i], modulus)
	}

	degInv := mod.InvMForm(mod.Inv(uint64(ringParams.degree), modulus), modulus)

	return &cyclotomicPow2Transformer{
		ringParams: ringParams,
		mod:        modulus,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		degInv: degInv,
	}
}

func (ntt *cyclotomicPow2Transformer) RingParameters() RingParameters {
	return ntt.ringParams
}

func (ntt *cyclotomicPow2Transformer) NTTInPlace(coeffs []uint64) {
	nttInPlacePow2(coeffs, ntt.tw, ntt.twS, ntt.mod.Value())
	mod.MFormVecTo(coeffs, ntt.mod, coeffs)
}

func (ntt *cyclotomicPow2Transformer) InvNTTInPlace(coeffs []uint64) {
	inttInPlacePow2(coeffs, ntt.twInv, ntt.twInvS, ntt.mod.Value())
	mod.ScalarMulVecTo(coeffs, ntt.degInv, ntt.mod, coeffs)
}

func (ntt *cyclotomicPow2Transformer) Type() TransformType {
	return PowerOfTwo
}
