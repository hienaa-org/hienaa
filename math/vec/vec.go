// Package vec implements vector operations acting on slices.
//
// Operations usually take two forms: for example,
//   - Add(v0, v1) adds v0, v1, allocates a new vector to store the result and returns it.
//   - AddTo(vOut, v0, v1) adds v0, v1 and writes the result to pre-allocated vOut without returning.
//
// Note that in most cases, v0, v1, and vOut can be identical,
// but partially overlapping slices may produce wrong results.
package vec

import (
	"math"

	"github.com/hienaa-org/hienaa/math/num"
)

// checkLength checks if all vectors have the same length,
// and panics if not.
func checkLength(xs ...int) {
	if len(xs) == 0 {
		return
	}

	for i := 1; i < len(xs); i++ {
		if xs[i] != xs[0] {
			panic("inconsistent input(s)")
		}
	}
}

// Gather returns a new vector with elements of given indices.
func Gather[T any](v []T, idx ...int) []T {
	vOut := make([]T, len(idx))
	for i := range idx {
		vOut[i] = v[idx[i]]
	}
	return vOut
}

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

// Range returns a vector containing [start, end).
func Range[T num.Integer](start, end T) []T {
	v := make([]T, uint64(end)-uint64(start))
	for i := range v {
		v[i] = start + T(i)
	}
	return v
}

// RadixReverseInPlace computes the radix-r reverse of v in-place.
// Assumes len(v) is a power of r.
func RadixReverseInPlace[T any](v []T, r int) {
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

	logN := int(math.Round(num.Log2(len(v)) / num.Log2(r)))
	for i := range v {
		idx, j := i, 0
		for k := 0; k < logN; k++ {
			j = j*r + (idx % r)
			idx /= r
		}
		if i < j {
			v[i], v[j] = v[j], v[i]
		}
	}
}

// Max returns max(v).
// If len(v) == 0, it returns 0.
func Max[T num.Real](v []T) T {
	switch len(v) {
	case 0:
		return 0
	case 1:
		return v[0]
	}

	r := v[0]
	for i := 1; i < len(v); i++ {
		if v[i] > r {
			r = v[i]
		}
	}
	return r
}

// Min returns min(v).
// If len(v) == 0, it returns 0.
func Min[T num.Real](v []T) T {
	switch len(v) {
	case 0:
		return 0
	case 1:
		return v[0]
	}

	r := v[0]
	for i := 1; i < len(v); i++ {
		if v[i] < r {
			r = v[i]
		}
	}
	return r
}
