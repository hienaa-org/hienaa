// Package vec implements vector operations acting on slices.
//
// Operations usually take two forms: for example,
//   - Add(v0, v1) adds v0, v1, allocates a new vector to store the result and returns it.
//   - AddTo(vOut, v0, v1) adds v0, v1 and writes the result to pre-allocated vOut without returning.
//
// Note that in most cases, v0, v1, and vOut can overlap.
// However, for operations that cannot, InPlace methods are implemented separately.
//
// For performance reasons, most functions in this package don't implement bound checks.
// If length mismatch happens, it may panic or produce wrong results.
package vec

import (
	"math"

	"github.com/hienaa-org/hienaa/math/num"
)

// Cast casts vector v of type []T to []TOut.
func Cast[TOut, T num.Real](v []T) []TOut {
	vOut := make([]TOut, len(v))
	CastTo(vOut, v)
	return vOut
}

// CastTo casts v of type []T to vOut of type []TOut.
func CastTo[TIn, TOut num.Real](vOut []TOut, vIn []TIn) {
	for i := range vIn {
		vOut[i] = TOut(vIn[i])
	}
}

// RadixReverseInPlace computes the radix-r reverse of v in-place.
// Assumes len(v) is a power of r.
func RadixReverseInPlace(v []uint64, r int) {
	if r == 2 {
		var bit, j int
		for i := 1; i < len(v); i++ {
			bit = len(v) >> 1
			for j >= bit {
				j -= bit
				bit >>= 1
			}
			j += bit
			if i < j {
				v[i], v[j] = v[j], v[i]
			}
		}
		return
	}

	logN := int(math.Round(math.Log(float64(len(v))) / math.Log(float64(r))))
	for i := 0; i < len(v); i++ {
		idx, j := i, 0
		for t := 0; t < logN; t++ {
			j = j*r + (idx % r)
			idx /= r
		}
		if i < j {
			v[i], v[j] = v[j], v[i]
		}
	}
}
