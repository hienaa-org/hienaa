// Package modops implements modular arithmetic operations for internal usage.
package modops

import (
	"math/bits"
)

// Unsigned represents the unsigned Integer type.
type Unsigned interface {
	~uint | ~uint8 | ~uint16 | ~uint32 | ~uint64
}

// Integer represents the Integer type.
type Integer interface {
	Unsigned | ~int | ~int8 | ~int16 | ~int32 | ~int64
}

// Add returns x0 + x1 mod q.
func Add(x0, x1, q uint64) uint64 {
	xOut := x0 + x1
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// Sub returns x0 - x1 mod q.
func Sub(x0, x1, q uint64) uint64 {
	xOut := x0 - x1
	if xOut >= q {
		xOut += q
	}
	return xOut
}

// Neg returns -x mod q.
func Neg(x, q uint64) uint64 {
	if x == 0 {
		return 0
	}
	return q - x
}

// BMod64 returns x mod q using Barrett reduction.
func BMod64(x, q, div, log uint64) uint64 {
	quo, _ := bits.Mul64(x>>log, div)
	xOut := x - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// BMod returns x mod q using Barrett reduction.
func BMod[T Integer](x T, q, div, log uint64) uint64 {
	if x < 0 {
		return Neg(BMod64(uint64(-int64(x)), q, div, log), q)
	}
	return BMod64(uint64(x), q, div, log)
}

// SForm transforms x into Shoup form.
func SForm(x, q uint64) uint64 {
	xS, _ := bits.Div64(x, 0, q)
	return xS
}

// SMul returns x0 * x1 mod q using Shoup multiplication.
func SMul(x0, x1, x1S, q uint64) uint64 {
	quo, _ := bits.Mul64(x0, x1S)

	xOut := x0*x1 - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// Reduce2Q reduces x assuming it is in [0, 2q).
func Reduce2Q(x, q uint64) uint64 {
	if x >= q {
		x -= q
	}
	return x
}

// Reduce4Q reduces x assuming it is in [0, 4q).
func Reduce4Q(x, q, twoQ uint64) uint64 {
	if x >= twoQ {
		x -= twoQ
	}
	if x >= q {
		x -= q
	}
	return x
}
