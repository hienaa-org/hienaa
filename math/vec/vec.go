// Package vec implements vector operations acting on slices.
//
// Operations usually take two forms: for example,
//   - Add(v0, v1) adds v0, v1, allocates a new vector to store the result and returns it.
//   - AddTo(v0, v1, vOut) adds v0, v1 and writes the result to pre-allocated vOut without returning.
//
// Note that in most cases, v0, v1, and vOut can overlap.
// However, for operations that cannot, InPlace methods are implemented separately.
//
// For performance reasons, most functions in this package don't implement bound checks.
// If length mismatch happens, it may panic or produce wrong results.
package vec

import "github.com/hienaa-org/hienaa/math/num"

// Cast casts vector v of type []T1 to []T2.
func Cast[T1, T2 num.Real](v []T1) []T2 {
	vOut := make([]T2, len(v))
	CastTo(v, vOut)
	return vOut
}

// CastTo casts v of type []T1 to vOut of type []T2.
func CastTo[T1, T2 num.Real](v []T1, vOut []T2) {
	for i := range vOut {
		vOut[i] = T2(v[i])
	}
}

// BitReverseInPlace reorders v into bit-reversal order in-place.
func BitReverseInPlace[T any](v []T) {
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
}
