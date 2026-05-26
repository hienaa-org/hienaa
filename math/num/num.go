// Package num implements various utility functions for modular and number arithmetic.
package num

import (
	"math"
	"math/bits"
)

// absUint64 returns the abs(x) in uint64.
func absUint64[T Integer](x T) uint64 {
	if x < 0 {
		return uint64(^x) + 1
	}
	return uint64(x)
}

// Abs returns the absolute value of x.
func Abs[T Integer](x T) T {
	return T(absUint64(x))
}

// IsPowerOfTwo returns whether x is a power of two.
func IsPowerOfTwo[T Integer](x T) bool {
	return (x > 0) && (x&(x-1)) == 0
}

// Log2 returns Log2(x). Panics if x <= 0.
func Log2[T Real](x T) float64 {
	if x <= 0 {
		panic("x must be positive")
	}

	return math.Log2(float64(x))
}

// GCD returns the greatest common divisor of x0 and x1.
func GCD[T Integer](x0, x1 T) T {
	return T(gcdUint64(absUint64(x0), absUint64(x1)))
}

// gcdUint64 returns the greatest common divisor of x0 and x1.
func gcdUint64(x0, x1 uint64) uint64 {
	switch {
	case x0 == 0:
		return x1
	case x1 == 0:
		return x0
	}

	i0 := bits.TrailingZeros64(x0)
	i1 := bits.TrailingZeros64(x1)
	k := min(i0, i1)

	x0 >>= i0
	x1 >>= i1

	for {
		if x0 > x1 {
			x0, x1 = x1, x0
		}

		x1 -= x0
		if x1 == 0 {
			return x0 << k
		}
		x1 >>= bits.TrailingZeros64(x1)
	}
}

// LCM returns the least common multiple of x0 and x1.
func LCM[T Integer](x0, x1 T) T {
	if x0 == 0 || x1 == 0 {
		return 0
	}

	return T(absUint64(x0/GCD(x0, x1)) * absUint64(x1))
}

// DivCeil returns ceil(x/y).
func DivCeil[T Integer](x, y T) T {
	q := x / y
	r := x % y

	if r != 0 && ((x < 0) == (y < 0)) {
		q++
	}
	return q
}

// DivRound returns round(x/y).
func DivRound[T Integer](x, y T) T {
	q := x / y
	r := x % y

	rr := absUint64(r)
	yy := absUint64(y)

	if rr >= (yy>>1)+(yy&1) {
		if (x < 0) == (y < 0) {
			q++
		} else {
			q--
		}
	}
	return q
}
