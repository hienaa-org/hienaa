package crt

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
)

// halfProductMixedRadix returns the mixed-radix digits of floor(prod(mod) / 2).
func halfProductMixedRadix(mod []*num.Modulus) []uint64 {
	half := make([]uint64, len(mod))
	isOdd := true
	for i := len(mod) - 1; i >= 0; i-- {
		qv := mod[i].Value()
		if isOdd {
			half[i] = qv >> 1
		}
		if qv&1 == 0 {
			isOdd = false
		}
	}
	return half
}

// embedToModOut returns balanced x mod qOut for x in [0, qIn).
func embedToModOut(x uint64, qOut, qOutDivHi, qIn, halfQIn uint64) uint64 {
	if x <= halfQIn {
		return modops.BMod64(x, qOut, qOutDivHi)
	}
	return modops.Neg(modops.BMod64(qIn-x, qOut, qOutDivHi), qOut)
}

// isMixedRadixNegative checks if i-th index of v is negative in signed mixed radix representation.
func isMixedRadixNegative[T *[embedBatch]uint64 | []uint64](v []T, i int, half []uint64) uint64 {
	for j := len(v) - 1; j >= 0; j-- {
		x := v[j][i]
		qHalf := half[j]

		if x > qHalf {
			return 1
		}
		if x < qHalf {
			return 0
		}
	}
	return 0
}

// getPolyFromPool fetches an *[Poly] with given rank, modLen and [dft.Form] from pool.
func getPolyFromPool(pool *pool.Pool[*[]uint64], rank, modLen int, outForm dft.Form) (eOut *Poly, put func()) {
	eOut = &Poly{
		Coeffs: make([][]uint64, modLen),
		Form:   outForm,
	}

	coeffsPtr := make([]*[]uint64, modLen)
	for i := 0; i < modLen; i++ {
		vPtr := pool.Get()
		coeffsPtr[i] = vPtr
		eOut.Coeffs[i] = (*vPtr)[:rank]
	}

	put = func() {
		for i := range coeffsPtr {
			pool.Put(coeffsPtr[i])
		}
	}

	return eOut, put
}
