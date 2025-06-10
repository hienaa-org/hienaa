// Package num implements various utility functions for arithmetic.
package num

import "math/bits"

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
