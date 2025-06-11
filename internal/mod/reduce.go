package mod

import "math/bits"

// Add computes x0 + x1 mod q.
func Add(x0, x1, q uint64) uint64 {
	xOut := x0 + x1
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// Sub computes x0 - x1 mod q.
func Sub(x0, x1, q uint64) uint64 {
	xOut := x0 - x1
	if xOut >= q {
		xOut += q
	}
	return xOut
}

// BMul computes x * y mod q using Barrett reduction.
func BMul(x0, y0, q, divHi, divLo uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, y0)

	quo := xOutHi * divHi

	quoLo, _ := bits.Mul64(xOutLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xOutLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xOutHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	xOut := xOutLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// BMulLazy computes x * y mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulLazy(x0, y0, q, divHi, divLo uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, y0)

	quo := xOutHi * divHi

	quoLo, _ := bits.Mul64(xOutLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xOutLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xOutHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	return xOutLo - quo*q
}

// BMod computes x mod q using Barrett reduction.
func BMod(xHi, xLo, q, divHi, divLo uint64) uint64 {
	quo := xHi * divHi

	quoLo, _ := bits.Mul64(xLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	xOut := xLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// BModLazy computes x mod q using Barrett reduction,
// but the result is in [0, 2q).
func BModLazy(xHi, xLo, q, divHi, divLo uint64) uint64 {
	quo := xHi * divHi

	quoLo, _ := bits.Mul64(xLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	return xLo - quo*q
}

// MForm transforms x into Montgomery form.
func MForm(x, q, divHi, divLo uint64) uint64 {
	xM, _ := bits.Mul64(x, divLo)
	xM += x * divHi

	xMOut := -xM * q
	if xMOut >= q {
		xMOut -= q
	}
	return xMOut
}

// InvMForm transforms xM to Normal form.
func InvMForm(xM, q, inv uint64) uint64 {
	x, _ := bits.Mul64(xM*inv, q)
	xOut := q - x
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// MMul computes x0 * x1 mod q in Montgomery form.
func MMul(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)

	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	xOutM := xOutMHi - wHi + q
	if xOutM >= q {
		xOutM -= q
	}
	return xOutM
}

// MMulLazy computes x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)

	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	return xOutMHi - wHi + q
}
