package crt

import (
	"math/big"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
)

// BigConverter converts scalar and polynomial [Element] values to *[big.Int] values.
type BigConverter struct {
	Op *Operator

	bigMod     *big.Int
	bigModHalf *big.Int
	gadget     []*big.Int

	pool *pool.Pool[*big.Int]
}

// NewBigConverter creates a new [BigConverter].
func NewBigConverter(op *Operator) *BigConverter {
	mod := make([]*big.Int, len(op.mod))
	bigMod := big.NewInt(1)
	for i := range op.mod {
		mod[i] = new(big.Int).SetUint64(op.mod[i].Value())
		bigMod.Mul(bigMod, mod[i])
	}
	bigModHalf := new(big.Int).Rsh(bigMod, 1)

	gadget := make([]*big.Int, len(op.mod))
	for i := range op.mod {
		qStar := new(big.Int).Div(bigMod, mod[i])
		qStarInv := new(big.Int).ModInverse(qStar, mod[i])
		gadget[i] = new(big.Int).Mul(qStar, qStarInv)
		gadget[i].Mod(gadget[i], bigMod)
	}

	return &BigConverter{
		Op: op,

		bigMod:     bigMod,
		bigModHalf: bigModHalf,
		gadget:     gadget,

		pool: pool.NewPool(func() *big.Int {
			return new(big.Int)
		}),
	}
}

// Convert converts e to a vector of *[big.Int].
func (cvt *BigConverter) Convert(e *Element) []*big.Int {
	var vOut []*big.Int
	if e.Type == TypeScalar {
		vOut = make([]*big.Int, 1)
	} else {
		vOut = make([]*big.Int, cvt.Op.params.Rank())
	}

	for i := range vOut {
		vOut[i] = new(big.Int)
	}
	cvt.ConvertTo(vOut, e)
	return vOut
}

// ConvertTo converts e to vOut.
func (cvt *BigConverter) ConvertTo(vOut []*big.Int, e *Element) {
	checkUnaryOperable(cvt.Op.params.Rank(), len(cvt.Op.mod), e, e)

	mul := cvt.pool.Get()
	defer cvt.pool.Put(mul)

	if e.Type == TypeScalar {
		if len(vOut) != 1 {
			panic("inconsistent inputs")
		}

		vOut[0].SetUint64(0)
		for l := range cvt.Op.mod {
			mul.SetUint64(e.Coeffs[l][0])
			mul.Mul(mul, cvt.gadget[l])
			vOut[0].Add(vOut[0], mul)
		}
		vOut[0].Mod(vOut[0], cvt.bigMod)
		if vOut[0].Cmp(cvt.bigModHalf) > 0 {
			vOut[0].Sub(vOut[0], cvt.bigMod)
		}
	} else {
		if e.Form == dft.FormNTT {
			panic("inputs must be in Coeff form")
		} else if len(vOut) != cvt.Op.params.Rank() {
			panic("inconsistent inputs")
		}

		for i := range cvt.Op.params.Rank() {
			vOut[i].SetUint64(0)
			for l := range cvt.Op.mod {
				mul.SetUint64(e.Coeffs[l][i])
				mul.Mul(mul, cvt.gadget[l])
				vOut[i].Add(vOut[i], mul)
			}
			vOut[i].Mod(vOut[i], cvt.bigMod)
			if vOut[i].Cmp(cvt.bigModHalf) > 0 {
				vOut[i].Sub(vOut[i], cvt.bigMod)
			}
		}
	}
}
