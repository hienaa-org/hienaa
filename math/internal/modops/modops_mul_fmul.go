//go:build arm64 || ppc64 || ppc64le || s390x || riscv64 || loong64

package modops

import (
	"math"
)

// In these CPUs, math.FMA and math.Floor is natively supported
// without runtime CPU flag checks.
// See: https://github.com/golang/go/issues/36351

// Mul returns x0 * x1 mod q using float64 reduction.
func Mul(x0, x1, q, divHi, divLo uint64, qf, qfInv float64) uint64 {
	x0f := float64(x0)
	x1f := float64(x1)

	hi := x0f * x1f
	lo := math.FMA(x0f, x1f, -hi)

	quo := math.Floor(hi * qfInv)
	xOutf := math.FMA(-quo, qf, hi) + lo + qf

	xOut := uint64(xOutf)
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// MulLazy returns x0 * x1 mod q using float64 reduction,
// but the result is in [0, 2q).
func MulLazy(x0, x1, q, divHi, divLo uint64, qf, qfInv float64) uint64 {
	x0f := float64(x0)
	x1f := float64(x1)

	hi := x0f * x1f
	lo := math.FMA(x0f, x1f, -hi)

	quo := math.Floor(hi * qfInv)
	xOutf := math.FMA(-quo, qf, hi) + lo + qf

	return uint64(xOutf)
}
