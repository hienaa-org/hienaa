// Package modops implements modular arithmetic operations for internal usage.
package modops

import (
	"math/bits"
)

// Unsigned represents the unsigned Integer type.
type Unsigned interface {
	uint | uint8 | uint16 | uint32 | uint64
}

// Signed represents the signed Integer type.
type Signed interface {
	int | int8 | int16 | int32 | int64
}

// Integer represents the Integer type.
type Integer interface {
	Unsigned | Signed
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

// BMul returns x0 * x1 mod q using Barrett reduction.
func BMul(x0, x1, q, divHi, divLo uint64) uint64 {
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

// BMulLazy returns x0 * x1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulLazy(x0, x1, q, divHi, divLo uint64) uint64 {
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

// BMod128 returns x mod q using Barrett reduction.
func BMod128(xHi, xLo, q, divHi, divLo uint64) uint64 {
	quo := xHi * divHi

	quoLo, _ := bits.Mul64(xLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	xOut := xLo - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// BMod128Lazy returns x mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMod128Lazy(xHi, xLo, q, divHi, divLo uint64) uint64 {
	quo := xHi * divHi

	quoLo, _ := bits.Mul64(xLo, divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoMid1Lo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	_, quoMidCarry = bits.Add64(quoMidSum, quoLo, 0)
	quo, _ = bits.Add64(quo, 0, quoMidCarry)

	return xLo - quo*q
}

// BMod64 returns x mod q using Barrett reduction.
func BMod64(x, q, divHi uint64) uint64 {
	quo, _ := bits.Mul64(x, divHi)
	xOut := x - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// BMod returns x mod q using Barrett reduction.
func BMod[T Integer](x T, q, divHi uint64) uint64 {
	if x < 0 {
		return Neg(BMod64(uint64(-int64(x)), q, divHi), q)
	}
	return BMod64(uint64(x), q, divHi)
}

// MForm transforms x into Montgomery form.
func MForm(x, q, divHi, divLo uint64) uint64 {
	xM, _ := bits.Mul64(x, divLo)
	xM += x * divHi

	xOutM := -xM * q
	if xOutM >= q {
		xOutM -= q
	}
	return xOutM
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

// MMul returns x0 * x1 mod q in Montgomery form.
// x0M and x1M must be a valid Montgomery form.
func MMul(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)
	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	xOutM := xOutMHi - wHi + q
	if xOutM >= q {
		xOutM -= q
	}
	return xOutM
}

// MMulLazy returns x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
// x0M and x1M must be a valid Montgomery form.
func MMulLazy(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)
	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	return xOutMHi - wHi + q
}

// SForm transforms x into Shoup form.
func SForm(x, q uint64) uint64 {
	xS, _ := bits.Div64(x, 0, q)
	return xS
}

// SMul returns x0 * x1 mod q using Shoup multiplication.
// x1S must be a valid Shoup form.
func SMul(x0, x1, x1S, q uint64) uint64 {
	quo, _ := bits.Mul64(x0, x1S)

	xOut := x0*x1 - quo*q
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// SMulLazy returns x0 * x1 mod q using Shoup multiplication,
// but the result is in [0, 2q).
// x1S must be a valid Shoup form.
func SMulLazy(x0, x1, x1S, q uint64) uint64 {
	quo, _ := bits.Mul64(x0, x1S)

	return x0*x1 - quo*q
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

// AnyAdd returns x0 + x1 mod q for arbitrary uint64 inputs.
func AnyAdd(x0, x1, q uint64) uint64 {
	xOut, carry := bits.Add64(x0, x1, 0)
	return bits.Rem64(carry, xOut, q)
}

// AnySub returns x0 - x1 mod q for arbitrary uint64 inputs.
func AnySub(x0, x1, q uint64) uint64 {
	xOut, borrow := bits.Sub64(x0, x1, 0)
	return bits.Rem64(borrow*(q-1), xOut, q)
}

// AnyNeg returns -x mod q for arbitrary uint64 inputs.
func AnyNeg(x, q uint64) uint64 {
	return Neg(x%q, q)
}

// AnyMul returns x0 * x1 mod q for arbitrary uint64 inputs.
func AnyMul(x0, x1, q uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)
	return bits.Rem64(xOutHi, xOutLo, q)
}
