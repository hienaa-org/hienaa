//go:build !(amd64 && !purego)

package vec

import (
	"math/bits"
	"unsafe"

	"github.com/hienaa-org/hienaa/math/num"
)

// AddTo computes vOut = v0 + v1 mod q.
func AddTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = addMod(w0[0], w1[0], qv)
		wOut[1] = addMod(w0[1], w1[1], qv)
		wOut[2] = addMod(w0[2], w1[2], qv)
		wOut[3] = addMod(w0[3], w1[3], qv)

		wOut[4] = addMod(w0[4], w1[4], qv)
		wOut[5] = addMod(w0[5], w1[5], qv)
		wOut[6] = addMod(w0[6], w1[6], qv)
		wOut[7] = addMod(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = addMod(v0[i], v1[i], qv)
	}
}

// AddLazyTo computes vOut = v0 + v1.
func AddLazyTo(v0, v1, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = w0[0] + w1[0]
		wOut[1] = w0[1] + w1[1]
		wOut[2] = w0[2] + w1[2]
		wOut[3] = w0[3] + w1[3]

		wOut[4] = w0[4] + w1[4]
		wOut[5] = w0[5] + w1[5]
		wOut[6] = w0[6] + w1[6]
		wOut[7] = w0[7] + w1[7]
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = v0[i] + v1[i]
	}
}

// SubTo computes vOut = v0 - v1 mod q.
func SubTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = subMod(w0[0], w1[0], qv)
		wOut[1] = subMod(w0[1], w1[1], qv)
		wOut[2] = subMod(w0[2], w1[2], qv)
		wOut[3] = subMod(w0[3], w1[3], qv)

		wOut[4] = subMod(w0[4], w1[4], qv)
		wOut[5] = subMod(w0[5], w1[5], qv)
		wOut[6] = subMod(w0[6], w1[6], qv)
		wOut[7] = subMod(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = subMod(v0[i], v1[i], qv)
	}
}

// SubLazyTo computes vOut = v0 - v1.
func SubLazyTo(v0, v1, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = w0[0] - w1[0]
		wOut[1] = w0[1] - w1[1]
		wOut[2] = w0[2] - w1[2]
		wOut[3] = w0[3] - w1[3]

		wOut[4] = w0[4] - w1[4]
		wOut[5] = w0[5] - w1[5]
		wOut[6] = w0[6] - w1[6]
		wOut[7] = w0[7] - w1[7]
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = v0[i] - v1[i]
	}
}

// BMulTo computes vOut = v0 * v1 mod q using Barrett reduction.
func BMulTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = bMul(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] = bMul(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] = bMul(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] = bMul(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] = bMul(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] = bMul(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] = bMul(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] = bMul(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = bMul(v0[i], v1[i], qv, divHi, divLo)
	}
}

// BMulAddTo computes vOut += v0 * v1 mod q using Barrett reduction.
func BMulAddTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = addMod(wOut[0], bMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = addMod(wOut[1], bMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = addMod(wOut[2], bMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = addMod(wOut[3], bMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = addMod(wOut[4], bMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = addMod(wOut[5], bMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = addMod(wOut[6], bMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = addMod(wOut[7], bMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = addMod(vOut[i], bMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
}

// BMulSubTo computes vOut -= v0 * v1 mod q using Barrett reduction.
func BMulSubTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = subMod(wOut[0], bMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = subMod(wOut[1], bMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = subMod(wOut[2], bMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = subMod(wOut[3], bMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = subMod(wOut[4], bMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = subMod(wOut[5], bMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = subMod(wOut[6], bMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = subMod(wOut[7], bMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = subMod(vOut[i], bMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
}

// bMul returns x0 * x1 mod q using Barrett reduction.
// See also [num.BMul].
func bMul(x0, x1, q, divHi, divLo uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

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

// BMulLazyTo computes vOut = v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = bMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] = bMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] = bMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] = bMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] = bMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] = bMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] = bMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] = bMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = bMulLazy(v0[i], v1[i], qv, divHi, divLo)
	}
}

// BMulAddLazyTo computes vOut += v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
func BMulAddLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += bMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += bMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += bMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += bMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += bMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += bMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += bMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += bMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += bMulLazy(v0[i], v1[i], qv, divHi, divLo)
	}
}

// BMulSubLazyTo computes vOut -= v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
func BMulSubLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += qv<<1 - bMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += qv<<1 - bMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += qv<<1 - bMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += qv<<1 - bMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += qv<<1 - bMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += qv<<1 - bMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += qv<<1 - bMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += qv<<1 - bMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += qv<<1 - bMulLazy(v0[i], v1[i], qv, divHi, divLo)
	}
}

// bMulLazy returns x0 * x1 mod q using Barrett reduction,
// but the result is in [0, 2q).
// See also [num.BMulLazy].
func bMulLazy(x0, x1, q, divHi, divLo uint64) uint64 {
	xOutHi, xOutLo := bits.Mul64(x0, x1)

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

// MMulTo computes vOut = v0 * v1 mod q in Montgomery form,
func MMulTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mMul(w0[0], w1[0], qv, inv)
		wOut[1] = mMul(w0[1], w1[1], qv, inv)
		wOut[2] = mMul(w0[2], w1[2], qv, inv)
		wOut[3] = mMul(w0[3], w1[3], qv, inv)

		wOut[4] = mMul(w0[4], w1[4], qv, inv)
		wOut[5] = mMul(w0[5], w1[5], qv, inv)
		wOut[6] = mMul(w0[6], w1[6], qv, inv)
		wOut[7] = mMul(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mMul(v0[i], v1[i], qv, inv)
	}
}

// MMulLazyTo computes vOut = v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] = mMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] = mMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] = mMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] = mMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] = mMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] = mMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] = mMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mMulLazy(v0[i], v1[i], qv, inv)
	}
}
