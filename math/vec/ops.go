package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/internal/mod"
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

// MForm returns v in Montgomery form.
func MForm(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MFormTo(v, q, vOut)
	return vOut
}

// MFormTo computes vOut as v in Montgomery form.
func MFormTo(v []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mod.MForm(w0[0], qv, divHi, divLo)
		wOut[1] = mod.MForm(w0[1], qv, divHi, divLo)
		wOut[2] = mod.MForm(w0[2], qv, divHi, divLo)
		wOut[3] = mod.MForm(w0[3], qv, divHi, divLo)

		wOut[4] = mod.MForm(w0[4], qv, divHi, divLo)
		wOut[5] = mod.MForm(w0[5], qv, divHi, divLo)
		wOut[6] = mod.MForm(w0[6], qv, divHi, divLo)
		wOut[7] = mod.MForm(w0[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.MForm(v[i], qv, divHi, divLo)
	}
}

// InvMForm transforms vM to Normal form.
func InvMForm(vM []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(vM))
	InvMFormTo(vM, q, vOut)
	return vOut
}

// InvMFormTo computes vOut as vM in Normal form.
func InvMFormTo(vM []uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mod.InvMForm(w0[0], qv, inv)
		wOut[1] = mod.InvMForm(w0[1], qv, inv)
		wOut[2] = mod.InvMForm(w0[2], qv, inv)
		wOut[3] = mod.InvMForm(w0[3], qv, inv)

		wOut[4] = mod.InvMForm(w0[4], qv, inv)
		wOut[5] = mod.InvMForm(w0[5], qv, inv)
		wOut[6] = mod.InvMForm(w0[6], qv, inv)
		wOut[7] = mod.InvMForm(w0[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.InvMForm(vM[i], qv, inv)
	}
}

// ScalarBMul returns c * v mod q using Barrett multiplication.
func ScalarBMul(v0 []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	ScalarBMulTo(v0, c, q, vOut)
	return vOut
}

// ScalarBMulLazy returns c * v mod q using Barrett multiplication,
// but the result is in [0, 2q).
func ScalarBMulLazy(v0 []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	ScalarBMulLazyTo(v0, c, q, vOut)
	return vOut
}

// ScalarMMul returns c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, the output is the same form as v.
func ScalarMMul(v0 []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	ScalarMMulTo(v0, cM, q, vOut)
	return vOut
}

// ScalarMMulTo computes vOut = c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulTo(v0 []uint64, cM uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mod.MMul(w0[0], cM, qv, inv)
		wOut[1] = mod.MMul(w0[1], cM, qv, inv)
		wOut[2] = mod.MMul(w0[2], cM, qv, inv)
		wOut[3] = mod.MMul(w0[3], cM, qv, inv)

		wOut[4] = mod.MMul(w0[4], cM, qv, inv)
		wOut[5] = mod.MMul(w0[5], cM, qv, inv)
		wOut[6] = mod.MMul(w0[6], cM, qv, inv)
		wOut[7] = mod.MMul(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.MMul(v0[i], cM, qv, inv)
	}
}

// ScalarMMulAddTo computes vOut += c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulAddTo(v0 []uint64, cM uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mod.Add(wOut[0], mod.MMul(w0[0], cM, qv, inv), qv)
		wOut[1] = mod.Add(wOut[1], mod.MMul(w0[1], cM, qv, inv), qv)
		wOut[2] = mod.Add(wOut[2], mod.MMul(w0[2], cM, qv, inv), qv)
		wOut[3] = mod.Add(wOut[3], mod.MMul(w0[3], cM, qv, inv), qv)

		wOut[4] = mod.Add(wOut[4], mod.MMul(w0[4], cM, qv, inv), qv)
		wOut[5] = mod.Add(wOut[5], mod.MMul(w0[5], cM, qv, inv), qv)
		wOut[6] = mod.Add(wOut[6], mod.MMul(w0[6], cM, qv, inv), qv)
		wOut[7] = mod.Add(wOut[7], mod.MMul(w0[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.Add(vOut[i], mod.MMul(v0[i], cM, qv, inv), qv)
	}
}

// ScalarMMulSubTo computes vOut -= c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulSubTo(v0 []uint64, cM uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mod.Sub(wOut[0], mod.MMul(w0[0], cM, qv, inv), qv)
		wOut[1] = mod.Sub(wOut[1], mod.MMul(w0[1], cM, qv, inv), qv)
		wOut[2] = mod.Sub(wOut[2], mod.MMul(w0[2], cM, qv, inv), qv)
		wOut[3] = mod.Sub(wOut[3], mod.MMul(w0[3], cM, qv, inv), qv)

		wOut[4] = mod.Sub(wOut[4], mod.MMul(w0[4], cM, qv, inv), qv)
		wOut[5] = mod.Sub(wOut[5], mod.MMul(w0[5], cM, qv, inv), qv)
		wOut[6] = mod.Sub(wOut[6], mod.MMul(w0[6], cM, qv, inv), qv)
		wOut[7] = mod.Sub(wOut[7], mod.MMul(w0[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.Sub(vOut[i], mod.MMul(v0[i], cM, qv, inv), qv)
	}
}

// ScalarMMulLazy returns c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, the output is the same form as v.
func ScalarMMulLazy(v0 []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	ScalarMMulLazyTo(v0, cM, q, vOut)
	return vOut
}

// ScalarMMulLazyTo computes vOut = c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulLazyTo(v0 []uint64, cM uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = mod.MMulLazy(w0[0], cM, qv, inv)
		wOut[1] = mod.MMulLazy(w0[1], cM, qv, inv)
		wOut[2] = mod.MMulLazy(w0[2], cM, qv, inv)
		wOut[3] = mod.MMulLazy(w0[3], cM, qv, inv)

		wOut[4] = mod.MMulLazy(w0[4], cM, qv, inv)
		wOut[5] = mod.MMulLazy(w0[5], cM, qv, inv)
		wOut[6] = mod.MMulLazy(w0[6], cM, qv, inv)
		wOut[7] = mod.MMulLazy(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.MMulLazy(v0[i], cM, qv, inv)
	}
}

// ScalarMMulAddLazyTo computes vOut += c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulAddLazyTo(v0 []uint64, cM uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += mod.MMulLazy(w0[0], cM, qv, inv)
		wOut[1] += mod.MMulLazy(w0[1], cM, qv, inv)
		wOut[2] += mod.MMulLazy(w0[2], cM, qv, inv)
		wOut[3] += mod.MMulLazy(w0[3], cM, qv, inv)

		wOut[4] += mod.MMulLazy(w0[4], cM, qv, inv)
		wOut[5] += mod.MMulLazy(w0[5], cM, qv, inv)
		wOut[6] += mod.MMulLazy(w0[6], cM, qv, inv)
		wOut[7] += mod.MMulLazy(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += mod.MMulLazy(v0[i], cM, qv, inv)
	}
}

// ScalarMMulSubLazyTo computes vOut -= c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulSubLazyTo(v0 []uint64, cM uint64, q *num.Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	cMNeg := qv - cM

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += mod.MMulLazy(w0[0], cMNeg, qv, inv)
		wOut[1] += mod.MMulLazy(w0[1], cMNeg, qv, inv)
		wOut[2] += mod.MMulLazy(w0[2], cMNeg, qv, inv)
		wOut[3] += mod.MMulLazy(w0[3], cMNeg, qv, inv)

		wOut[4] += mod.MMulLazy(w0[4], cMNeg, qv, inv)
		wOut[5] += mod.MMulLazy(w0[5], cMNeg, qv, inv)
		wOut[6] += mod.MMulLazy(w0[6], cMNeg, qv, inv)
		wOut[7] += mod.MMulLazy(w0[7], cMNeg, qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += mod.MMulLazy(v0[i], cMNeg, qv, inv)
	}
}

// MMul returns v0 * v1 mod q in Montgomery form.
func MMul(v0M, v1M []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0M))
	MMulTo(v0M, v1M, q, vOut)
	return vOut
}

// MMulTo computes vOut = v0 * v1 mod q in Montgomery form,
func MMulTo(v0M, v1M []uint64, q *num.Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mod.MMul(w0[0], w1[0], qv, inv)
		wOut[1] = mod.MMul(w0[1], w1[1], qv, inv)
		wOut[2] = mod.MMul(w0[2], w1[2], qv, inv)
		wOut[3] = mod.MMul(w0[3], w1[3], qv, inv)

		wOut[4] = mod.MMul(w0[4], w1[4], qv, inv)
		wOut[5] = mod.MMul(w0[5], w1[5], qv, inv)
		wOut[6] = mod.MMul(w0[6], w1[6], qv, inv)
		wOut[7] = mod.MMul(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mod.MMul(v0M[i], v1M[i], qv, inv)
	}
}

// MMulAddTo computes vOut += v0 * v1 mod q in Montgomery form.
func MMulAddTo(v0M, v1M []uint64, q *num.Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mod.Add(wOut[0], mod.MMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = mod.Add(wOut[1], mod.MMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = mod.Add(wOut[2], mod.MMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = mod.Add(wOut[3], mod.MMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = mod.Add(wOut[4], mod.MMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = mod.Add(wOut[5], mod.MMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = mod.Add(wOut[6], mod.MMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = mod.Add(wOut[7], mod.MMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mod.Add(vOutM[i], mod.MMul(v0M[i], v1M[i], qv, inv), qv)
	}
}

// MMulSubTo computes vOut -= v0 * v1 mod q in Montgomery form.
func MMulSubTo(v0M, v1M []uint64, q *num.Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mod.Sub(wOut[0], mod.MMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = mod.Sub(wOut[1], mod.MMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = mod.Sub(wOut[2], mod.MMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = mod.Sub(wOut[3], mod.MMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = mod.Sub(wOut[4], mod.MMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = mod.Sub(wOut[5], mod.MMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = mod.Sub(wOut[6], mod.MMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = mod.Sub(wOut[7], mod.MMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mod.Sub(vOutM[i], mod.MMul(v0M[i], v1M[i], qv, inv), qv)
	}
}

// MMulLazy returns v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(v0M, v1M []uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	MMulLazyTo(v0M, v1M, q, vOutM)
	return vOutM
}

// MMulLazyTo computes vOut = v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazyTo(v0M, v1M []uint64, q *num.Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mod.MMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] = mod.MMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] = mod.MMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] = mod.MMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] = mod.MMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] = mod.MMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] = mod.MMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] = mod.MMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mod.MMulLazy(v0M[i], v1M[i], qv, inv)
	}
}

// MMulAddLazyTo computes vOut += v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulAddLazyTo(v0M, v1M []uint64, q *num.Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] += mod.MMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] += mod.MMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] += mod.MMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] += mod.MMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] += mod.MMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] += mod.MMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] += mod.MMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] += mod.MMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += mod.MMulLazy(v0M[i], v1M[i], qv, inv)
	}
}

// MMulSubLazyTo computes vOut -= v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulSubLazyTo(v0M, v1M []uint64, q *num.Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] += mod.MMulLazy(qv-w0[0], w1[0], qv, inv)
		wOut[1] += mod.MMulLazy(qv-w0[1], w1[1], qv, inv)
		wOut[2] += mod.MMulLazy(qv-w0[2], w1[2], qv, inv)
		wOut[3] += mod.MMulLazy(qv-w0[3], w1[3], qv, inv)

		wOut[4] += mod.MMulLazy(qv-w0[4], w1[4], qv, inv)
		wOut[5] += mod.MMulLazy(qv-w0[5], w1[5], qv, inv)
		wOut[6] += mod.MMulLazy(qv-w0[6], w1[6], qv, inv)
		wOut[7] += mod.MMulLazy(qv-w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += mod.MMulLazy(qv-v0M[i], v1M[i], qv, inv)
	}
}
