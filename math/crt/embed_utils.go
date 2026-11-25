package crt

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/internal/modops"
)

// embedToModOut returns balanced x mod qOut for x in [0, qIn).
func embedToModOut(x uint64, qOut, qOutDivHi, qIn, halfQIn uint64) uint64 {
	if x <= halfQIn {
		return modops.BMod64(x, qOut, qOutDivHi)
	}
	return modops.Neg(modops.BMod64(qIn-x, qOut, qOutDivHi), qOut)
}

// embedToModOutNoReduce returns balanced x mod qOut for x in [0, qIn), when qOut > qIn.
func embedToModOutNoReduce(x uint64, qOut, qIn, halfQIn uint64) uint64 {
	if x > halfQIn {
		x += qOut - qIn
	}
	return x
}

// getPolyFromPool fetches a polynomial [Element] with given rank, modLen and [dft.Form] from pool.
func getPolyFromPool(pool *pool.Pool[*[]uint64], rank, modLen int, outForm dft.Form) (eOut *Element, put func()) {
	eOut = &Element{
		Coeffs: make([][]uint64, modLen),
		Form:   outForm,
		Type:   TypePoly,
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
