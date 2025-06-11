package vec

import (
	"math/bits"
	"unsafe"

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

// MMul returns v0 * v1 mod q in Montgomery form.
func MMul(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	MMulTo(v0, v1, q, vOut)
	return vOut
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

// MMulAddTo computes vOut += v0 * v1 mod q in Montgomery form.
func MMulAddTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = addMod(wOut[0], mMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = addMod(wOut[1], mMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = addMod(wOut[2], mMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = addMod(wOut[3], mMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = addMod(wOut[4], mMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = addMod(wOut[5], mMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = addMod(wOut[6], mMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = addMod(wOut[7], mMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = addMod(vOut[i], mMul(v0[i], v1[i], qv, inv), qv)
	}
}

// MMulSubTo computes vOut -= v0 * v1 mod q in Montgomery form.
func MMulSubTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = subMod(wOut[0], mMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = subMod(wOut[1], mMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = subMod(wOut[2], mMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = subMod(wOut[3], mMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = subMod(wOut[4], mMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = subMod(wOut[5], mMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = subMod(wOut[6], mMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = subMod(wOut[7], mMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = subMod(vOut[i], mMul(v0[i], v1[i], qv, inv), qv)
	}
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

// MMulLazy returns v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	MMulLazyTo(v0, v1, q, vOut)
	return vOut
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

// MMulAddLazyTo computes vOut += v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulAddLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += mMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] += mMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] += mMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] += mMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] += mMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] += mMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] += mMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] += mMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += mMulLazy(v0[i], v1[i], qv, inv)
	}
}

// MMulSubLazyTo computes vOut -= v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulSubLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += qv<<1 - mMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] += qv<<1 - mMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] += qv<<1 - mMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] += qv<<1 - mMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] += qv<<1 - mMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] += qv<<1 - mMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] += qv<<1 - mMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] += qv<<1 - mMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += qv<<1 - mMulLazy(v0[i], v1[i], qv, inv)
	}
}

// mMulLazy returns x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
// See also [num.MMulLazy].
func mMulLazy(x0M, x1M, q, inv uint64) uint64 {
	xOutMHi, xOutMLo := bits.Mul64(x0M, x1M)

	wHi, _ := bits.Mul64(xOutMLo*inv, q)

	return xOutMHi - wHi + q
}
