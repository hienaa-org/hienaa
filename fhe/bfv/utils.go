package bfv

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// computeScalingFactor computes the scaling factors for the given base and message moduli.
func computeScalingFactor(baseMod []*num.Modulus, msgMod *num.Modulus) []*rlwe.Element {
	scFacs := make([]*rlwe.Element, len(baseMod))

	msgModBig := big.NewInt(int64(msgMod.Value()))
	baseModBig := big.NewInt(1)
	scFacBig := big.NewInt(0)
	for i := range baseMod {
		baseModBig.Mul(baseModBig, big.NewInt(int64(baseMod[i].Value())))
		scFacBig.Quo(baseModBig, msgModBig)
		scFacs[i] = rlwe.NewElementFrom(crt.NewScalarFrom(scFacBig, baseMod), 0)
	}

	return scFacs
}
