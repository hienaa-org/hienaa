//go:build !(arm64 || ppc64 || ppc64le || s390x || riscv64 || loong64)

package modops

import "math/bits"

// Mul returns x0 * x1 mod q using Barrett reduction.
func Mul(x0, x1, q, div, log uint64, qf, qfInv float64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

	quo, _ := bits.Mul64(xOutHi<<(64-log)+xOutLo>>log, div)

	xOut := xOutLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}
