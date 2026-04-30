package crt

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/num"
)

// embedToModOut returns sign(x) mod qOut for x in [0, qIn).
func embedToModOut(x uint64, qOut *num.Modulus, qIn, halfQIn uint64) uint64 {
	if x <= halfQIn {
		return num.Reduce(x, qOut)
	}
	return num.Neg(num.Reduce(qIn-x, qOut), qOut)
}

// isMixedRadixNegative checks if i-th index of v is negative in signed mixed radix representation.
func isMixedRadixNegative[T *[embedBatch]uint64 | []uint64](v []T, i int, modInHalf []uint64) uint64 {
	for j := len(v) - 1; j >= 0; j-- {
		x := v[j][i]
		qHalf := modInHalf[j]

		if x > qHalf {
			return 1
		}
		if x < qHalf {
			return 0
		}
	}
	return 0
}

// getElementFromPool fetches an [Element] with given modLen and isNTT flag from pool.
func getElementFromPool(pool *pool.Pool[*[]uint64], modLen int, isNTT bool) (eOut *Element, put func()) {
	eOut = &Element{
		Coeffs: make([][]uint64, modLen),
		IsNTT:  isNTT,
	}

	coeffsPtr := make([]*[]uint64, modLen)
	for i := 0; i < modLen; i++ {
		vPtr := pool.Get()
		coeffsPtr[i] = vPtr
		eOut.Coeffs[i] = *vPtr
	}

	put = func() {
		for i := range coeffsPtr {
			pool.Put(coeffsPtr[i])
		}
	}

	return eOut, put
}
