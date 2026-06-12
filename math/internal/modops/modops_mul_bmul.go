//go:build !(arm64 || ppc64 || ppc64le || s390x || riscv64 || loong64)

package modops

import "math/bits"

// Mul returns x0 * x1 mod q using Barrett reduction.
func Mul(x0, x1, q, divHi, divLo uint64, qf, qfInv float64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

	quo := xOutHi * divHi

	quoLo, _ := bits.Mul64(xOutLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xOutLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xOutHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	xOut := xOutLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// MulLazy returns x0 * x1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func MulLazy(x0, x1, q, divHi, divLo uint64, qf, qfInv float64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

	quo := xOutHi * divHi

	quoLo, _ := bits.Mul64(xOutLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xOutLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xOutHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	return xOutLo - quo*q
}
