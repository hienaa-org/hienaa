package mod

import "math/bits"

// add returns x0 + x1 mod q.
func add(x0, x1, q uint64) uint64 {
	xOut := x0 + x1
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// sub returns x0 - x1 mod q.
func sub(x0, x1, q uint64) uint64 {
	xOut := x0 - x1
	if xOut >= q {
		xOut += q
	}
	return xOut
}

// bMul returns x0 * x1 mod q using Barrett reduction.
func bMul(x0, x1, q, divHi, divLo uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

	quo := xOutHi * divHi

	quoLo, _ := bits.Mul64(xOutLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xOutLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xOutHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo += quoMidCarry

	xOut := xOutLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// bMulLazy returns x0 * x1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func bMulLazy(x0, x1, q, divHi, divLo uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

	quo := xOutHi * divHi

	quoLo, _ := bits.Mul64(xOutLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xOutLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xOutHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo += quoMidCarry

	return xOutLo - quo*q
}

// bMod128 returns x mod q using Barrett reduction.
func bMod128(xHi, xLo, q, divHi, divLo uint64) uint64 {
	quo := xHi * divHi

	quoLo, _ := bits.Mul64(xLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo += quoMidCarry

	xOut := xLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// bMod64 returns x mod q using Barrett reduction.
func bMod64(x, q, divHi uint64) uint64 {
	quo, _ := bits.Mul64(x, divHi)
	xOut := x - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// bMod128Lazy returns x mod q using Barrett reduction,
// but the result is in [0, 2q).
func bMod128Lazy(xHi, xLo, q, divHi, divLo uint64) uint64 {
	quo := xHi * divHi

	quoLo, _ := bits.Mul64(xLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo += quoMidCarry

	return xLo - quo*q
}

// mForm transforms x into Montgomery form.
func mForm(x, q, divHi, divLo uint64) uint64 {
	xM, _ := bits.Mul64(x, divLo)
	xM += x * divHi

	xMOut := -xM * q
	if xMOut >= q {
		xMOut -= q
	}
	return xMOut
}

// invMForm transforms xM to Normal form.
func invMForm(xM, q, inv uint64) uint64 {
	x, _ := bits.Mul64(xM*inv, q)

	xOut := q - x
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// mMul returns x0 * x1 mod q in Montgomery form.
func mMul(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)
	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	xOutM := xOutMHi - wHi + q
	if xOutM >= q {
		xOutM -= q
	}
	return xOutM
}

// mMulLazy returns x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
func mMulLazy(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)
	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	return xOutMHi - wHi + q
}

// sForm transforms x into Shoup form.
func sForm(x, q uint64) uint64 {
	xS, _ := bits.Div64(x, 0, q)
	return xS
}

// sMul returns x0 * x1 mod q using Shoup multiplication.
func sMul(x0, x1, x1S, q uint64) uint64 {
	quo, _ := bits.Mul64(x0, x1S)

	xOut := x0*x1 - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// sMulLazy returns x0 * x1 mod q using Shoup multiplication,
// but the result is in [0, 2q).
func sMulLazy(x0, x1, x1S, q uint64) uint64 {
	quo, _ := bits.Mul64(x0, x1S)

	return x0*x1 - quo*q
}
