// Package num implements various utility functions for arithmetic.
package num

import (
	"math/bits"
)

// IsPowerOfTwo returns whether x is a power of two.
func IsPowerOfTwo(x uint64) bool {
	return (x > 0) && (x&(x-1)) == 0
}

// Log2 returns floor(log2(x)). Panics if x <= 0.
func Log2(x uint64) int {
	if x <= 0 {
		panic("log2: x must be positive")
	}

	return int(bits.Len64(x)) - 1
}

// GCD returns the greatest common divisor of x0 and x1.
func GCD(x0, x1 uint64) uint64 {
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
func LCM(x0, x1 uint64) uint64 {
	if x0 == 0 || x1 == 0 {
		return 0
	}

	return (x0 / GCD(x0, x1)) * x1
}

// Exp returns x**e.
func Exp(x, e uint64) uint64 {
	switch x {
	case 0:
		return 0
	case 1:
		return 1
	}

	switch e {
	case 0:
		return 1
	case 1:
		return x
	}

	r := uint64(1)
	for e > 0 {
		if e&1 == 1 {
			r *= x
		}
		x *= x
		e >>= 1
	}
	return r
}
