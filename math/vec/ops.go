package vec

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/math/num"
)

// Add returns v0 + v1 mod q.
func Add(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	AddTo(v0, v1, q, vOut)
	return vOut
}

// AddLazy returns v0 + v1.
func AddLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	AddLazyTo(v0, v1, vOut)
	return vOut
}

// addMod returns x0 + x1 mod q.
// See also [num.Add].
func addMod(x0, x1 uint64, q uint64) uint64 {
	xOut := x0 + x1
	if xOut >= q {
		xOut -= q
	}
	return xOut
}

// Sub returns v0 - v1 mod q.
func Sub(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SubTo(v0, v1, q, vOut)
	return vOut
}

// SubLazy returns v0 - v1.
func SubLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	SubLazyTo(v0, v1, vOut)
	return vOut
}

// subMod returns x0 - x1 mod q.
// See also [num.Sub].
func subMod(x0, x1 uint64, q uint64) uint64 {
	xOut := x0 - x1
	if xOut >= q {
		xOut += q
	}
	return xOut
}

// mMul returns x0 * x1 mod q in Montgomery form.
// See also [num.MMul].
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
// See also [num.MMulLazy].
func mMulLazy(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)

	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	return xOutMHi - wHi + q
}
