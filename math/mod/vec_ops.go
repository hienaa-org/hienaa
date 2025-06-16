package mod

import (
	"unsafe"
)

// AddVec returns v0 + v1 mod q.
func AddVec(v0, v1 []uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	AddVecTo(v0, v1, q, vOut)
	return vOut
}

// AddVecLazy returns v0 + v1.
func AddVecLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	AddLazyVecTo(v0, v1, vOut)
	return vOut
}

// SubVec returns v0 - v1 mod q.
func SubVec(v0, v1 []uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SubVecTo(v0, v1, q, vOut)
	return vOut
}

// SubVecLazy returns v0 - v1.
func SubVecLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	SubLazyVecTo(v0, v1, vOut)
	return vOut
}

// MFormVec returns v in Montgomery form.
func MFormVec(v []uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MFormVecTo(v, q, vOut)
	return vOut
}

// MFormVecTo computes vOutM as v in Montgomery form.
func MFormVecTo(v []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mForm(w0[0], qv, divHi, divLo)
		wOut[1] = mForm(w0[1], qv, divHi, divLo)
		wOut[2] = mForm(w0[2], qv, divHi, divLo)
		wOut[3] = mForm(w0[3], qv, divHi, divLo)

		wOut[4] = mForm(w0[4], qv, divHi, divLo)
		wOut[5] = mForm(w0[5], qv, divHi, divLo)
		wOut[6] = mForm(w0[6], qv, divHi, divLo)
		wOut[7] = mForm(w0[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mForm(v[i], qv, divHi, divLo)
	}
}

// InvMFormVec transforms vM to Normal form.
func InvMFormVec(vM []uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(vM))
	InvMFormVecTo(vM, q, vOut)
	return vOut
}

// InvMFormVecTo computes vOut as vM in Normal form.
func InvMFormVecTo(vM []uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = invMForm(w0[0], qv, inv)
		wOut[1] = invMForm(w0[1], qv, inv)
		wOut[2] = invMForm(w0[2], qv, inv)
		wOut[3] = invMForm(w0[3], qv, inv)

		wOut[4] = invMForm(w0[4], qv, inv)
		wOut[5] = invMForm(w0[5], qv, inv)
		wOut[6] = invMForm(w0[6], qv, inv)
		wOut[7] = invMForm(w0[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = invMForm(vM[i], qv, inv)
	}
}

// ScalarMulVec returns c * v mod q using Shoup multiplication.
func ScalarMulVec(v0 []uint64, c uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	ScalarMulVecTo(v0, c, q, vOut)
	return vOut
}

// ScalarMulVecTo computes vOut = c * v mod q using Shoup multiplication.
func ScalarMulVecTo(v0 []uint64, c uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := sForm(c, qv)

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = sMul(w0[0], c, cS, qv)
		wOut[1] = sMul(w0[1], c, cS, qv)
		wOut[2] = sMul(w0[2], c, cS, qv)
		wOut[3] = sMul(w0[3], c, cS, qv)

		wOut[4] = sMul(w0[4], c, cS, qv)
		wOut[5] = sMul(w0[5], c, cS, qv)
		wOut[6] = sMul(w0[6], c, cS, qv)
		wOut[7] = sMul(w0[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = sMul(v0[i], c, cS, qv)
	}
}

// ScalarMulAddVecTo computes vOut += c * v mod q using Shoup multiplication.
func ScalarMulAddVecTo(v0 []uint64, c uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := sForm(c, qv)

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = add(wOut[0], sMul(w0[0], c, cS, qv), qv)
		wOut[1] = add(wOut[1], sMul(w0[1], c, cS, qv), qv)
		wOut[2] = add(wOut[2], sMul(w0[2], c, cS, qv), qv)
		wOut[3] = add(wOut[3], sMul(w0[3], c, cS, qv), qv)

		wOut[4] = add(wOut[4], sMul(w0[4], c, cS, qv), qv)
		wOut[5] = add(wOut[5], sMul(w0[5], c, cS, qv), qv)
		wOut[6] = add(wOut[6], sMul(w0[6], c, cS, qv), qv)
		wOut[7] = add(wOut[7], sMul(w0[7], c, cS, qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = add(vOut[i], sMul(v0[i], c, cS, qv), qv)
	}
}

// ScalarMulSubVecTo computes vOut -= c * v mod q using Shoup multiplication.
func ScalarMulSubVecTo(v0 []uint64, c uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := sForm(c, qv)

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = sub(wOut[0], sMul(w0[0], c, cS, qv), qv)
		wOut[1] = sub(wOut[1], sMul(w0[1], c, cS, qv), qv)
		wOut[2] = sub(wOut[2], sMul(w0[2], c, cS, qv), qv)
		wOut[3] = sub(wOut[3], sMul(w0[3], c, cS, qv), qv)

		wOut[4] = sub(wOut[4], sMul(w0[4], c, cS, qv), qv)
		wOut[5] = sub(wOut[5], sMul(w0[5], c, cS, qv), qv)
		wOut[6] = sub(wOut[6], sMul(w0[6], c, cS, qv), qv)
		wOut[7] = sub(wOut[7], sMul(w0[7], c, cS, qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = sub(vOut[i], sMul(v0[i], c, cS, qv), qv)
	}
}

// ScalarMulLazyVec returns c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
func ScalarMulLazyVec(v0 []uint64, c uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	ScalarMulLazyVecTo(v0, c, q, vOut)
	return vOut
}

// ScalarMulLazyVecTo computes vOut = c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
func ScalarMulLazyVecTo(v0 []uint64, c uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := sForm(c, qv)

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = sMulLazy(w0[0], c, cS, qv)
		wOut[1] = sMulLazy(w0[1], c, cS, qv)
		wOut[2] = sMulLazy(w0[2], c, cS, qv)
		wOut[3] = sMulLazy(w0[3], c, cS, qv)

		wOut[4] = sMulLazy(w0[4], c, cS, qv)
		wOut[5] = sMulLazy(w0[5], c, cS, qv)
		wOut[6] = sMulLazy(w0[6], c, cS, qv)
		wOut[7] = sMulLazy(w0[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = sMulLazy(v0[i], c, cS, qv)
	}
}

// ScalarMulAddVecTo computes vOut += c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
func ScalarMulAddLazyVecTo(v0 []uint64, c uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := sForm(c, qv)

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += sMulLazy(w0[0], c, cS, qv)
		wOut[1] += sMulLazy(w0[1], c, cS, qv)
		wOut[2] += sMulLazy(w0[2], c, cS, qv)
		wOut[3] += sMulLazy(w0[3], c, cS, qv)

		wOut[4] += sMulLazy(w0[4], c, cS, qv)
		wOut[5] += sMulLazy(w0[5], c, cS, qv)
		wOut[6] += sMulLazy(w0[6], c, cS, qv)
		wOut[7] += sMulLazy(w0[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += sMulLazy(v0[i], c, cS, qv)
	}
}

// ScalarMulSubLazyVecTo computes vOut -= c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
func ScalarMulSubLazyVecTo(v0 []uint64, c uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cNeg := qv - c
	cNegS := sForm(cNeg, qv)

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += sMulLazy(w0[0], cNeg, cNegS, qv)
		wOut[1] += sMulLazy(w0[1], cNeg, cNegS, qv)
		wOut[2] += sMulLazy(w0[2], cNeg, cNegS, qv)
		wOut[3] += sMulLazy(w0[3], cNeg, cNegS, qv)

		wOut[4] += sMulLazy(w0[4], cNeg, cNegS, qv)
		wOut[5] += sMulLazy(w0[5], cNeg, cNegS, qv)
		wOut[6] += sMulLazy(w0[6], cNeg, cNegS, qv)
		wOut[7] += sMulLazy(w0[7], cNeg, cNegS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += sMulLazy(v0[i], cNeg, cNegS, qv)
	}
}

// ScalarMMulVec returns c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, the output is the same form as v.
func ScalarMMulVec(v0M []uint64, cM uint64, q *Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	ScalarMMulVecTo(v0M, cM, q, vOutM)
	return vOutM
}

// ScalarMMulVecTo computes vOut = c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulVecTo(v0M []uint64, cM uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mMul(w0[0], cM, qv, inv)
		wOut[1] = mMul(w0[1], cM, qv, inv)
		wOut[2] = mMul(w0[2], cM, qv, inv)
		wOut[3] = mMul(w0[3], cM, qv, inv)

		wOut[4] = mMul(w0[4], cM, qv, inv)
		wOut[5] = mMul(w0[5], cM, qv, inv)
		wOut[6] = mMul(w0[6], cM, qv, inv)
		wOut[7] = mMul(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mMul(v0M[i], cM, qv, inv)
	}
}

// ScalarMMulAddVecTo computes vOut += c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulAddVecTo(v0M []uint64, cM uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = add(wOut[0], mMul(w0[0], cM, qv, inv), qv)
		wOut[1] = add(wOut[1], mMul(w0[1], cM, qv, inv), qv)
		wOut[2] = add(wOut[2], mMul(w0[2], cM, qv, inv), qv)
		wOut[3] = add(wOut[3], mMul(w0[3], cM, qv, inv), qv)

		wOut[4] = add(wOut[4], mMul(w0[4], cM, qv, inv), qv)
		wOut[5] = add(wOut[5], mMul(w0[5], cM, qv, inv), qv)
		wOut[6] = add(wOut[6], mMul(w0[6], cM, qv, inv), qv)
		wOut[7] = add(wOut[7], mMul(w0[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = add(vOutM[i], mMul(v0M[i], cM, qv, inv), qv)
	}
}

// ScalarMMulSubVecTo computes vOut -= c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulSubVecTo(v0M []uint64, cM uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = sub(wOut[0], mMul(w0[0], cM, qv, inv), qv)
		wOut[1] = sub(wOut[1], mMul(w0[1], cM, qv, inv), qv)
		wOut[2] = sub(wOut[2], mMul(w0[2], cM, qv, inv), qv)
		wOut[3] = sub(wOut[3], mMul(w0[3], cM, qv, inv), qv)

		wOut[4] = sub(wOut[4], mMul(w0[4], cM, qv, inv), qv)
		wOut[5] = sub(wOut[5], mMul(w0[5], cM, qv, inv), qv)
		wOut[6] = sub(wOut[6], mMul(w0[6], cM, qv, inv), qv)
		wOut[7] = sub(wOut[7], mMul(w0[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = sub(vOutM[i], mMul(v0M[i], cM, qv, inv), qv)
	}
}

// ScalarMMulLazyVec returns c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, the output is the same form as v.
func ScalarMMulLazyVec(v0M []uint64, cM uint64, q *Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	ScalarMMulLazyVecTo(v0M, cM, q, vOutM)
	return vOutM
}

// ScalarMMulLazyVecTo computes vOut = c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulLazyVecTo(v0M []uint64, cM uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mMulLazy(w0[0], cM, qv, inv)
		wOut[1] = mMulLazy(w0[1], cM, qv, inv)
		wOut[2] = mMulLazy(w0[2], cM, qv, inv)
		wOut[3] = mMulLazy(w0[3], cM, qv, inv)

		wOut[4] = mMulLazy(w0[4], cM, qv, inv)
		wOut[5] = mMulLazy(w0[5], cM, qv, inv)
		wOut[6] = mMulLazy(w0[6], cM, qv, inv)
		wOut[7] = mMulLazy(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mMulLazy(v0M[i], cM, qv, inv)
	}
}

// ScalarMMulAddLazyVecTo computes vOut += c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulAddLazyVecTo(v0M []uint64, cM uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] += mMulLazy(w0[0], cM, qv, inv)
		wOut[1] += mMulLazy(w0[1], cM, qv, inv)
		wOut[2] += mMulLazy(w0[2], cM, qv, inv)
		wOut[3] += mMulLazy(w0[3], cM, qv, inv)

		wOut[4] += mMulLazy(w0[4], cM, qv, inv)
		wOut[5] += mMulLazy(w0[5], cM, qv, inv)
		wOut[6] += mMulLazy(w0[6], cM, qv, inv)
		wOut[7] += mMulLazy(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += mMulLazy(v0M[i], cM, qv, inv)
	}
}

// ScalarMMulSubLazyVecTo computes vOut -= c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulSubLazyVecTo(v0M []uint64, cM uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	cMNeg := qv - cM

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] += mMulLazy(w0[0], cMNeg, qv, inv)
		wOut[1] += mMulLazy(w0[1], cMNeg, qv, inv)
		wOut[2] += mMulLazy(w0[2], cMNeg, qv, inv)
		wOut[3] += mMulLazy(w0[3], cMNeg, qv, inv)

		wOut[4] += mMulLazy(w0[4], cMNeg, qv, inv)
		wOut[5] += mMulLazy(w0[5], cMNeg, qv, inv)
		wOut[6] += mMulLazy(w0[6], cMNeg, qv, inv)
		wOut[7] += mMulLazy(w0[7], cMNeg, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += mMulLazy(v0M[i], cMNeg, qv, inv)
	}
}

// MMulVec returns v0 * v1 mod q in Montgomery form.
func MMulVec(v0M, v1M []uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v0M))
	MMulVecTo(v0M, v1M, q, vOut)
	return vOut
}

// MMulVecTo computes vOut = v0 * v1 mod q in Montgomery form,
func MMulVecTo(v0M, v1M []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mMul(w0[0], w1[0], qv, inv)
		wOut[1] = mMul(w0[1], w1[1], qv, inv)
		wOut[2] = mMul(w0[2], w1[2], qv, inv)
		wOut[3] = mMul(w0[3], w1[3], qv, inv)

		wOut[4] = mMul(w0[4], w1[4], qv, inv)
		wOut[5] = mMul(w0[5], w1[5], qv, inv)
		wOut[6] = mMul(w0[6], w1[6], qv, inv)
		wOut[7] = mMul(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mMul(v0M[i], v1M[i], qv, inv)
	}
}

// MMulAddVecTo computes vOut += v0 * v1 mod q in Montgomery form.
func MMulAddVecTo(v0M, v1M []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = add(wOut[0], mMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = add(wOut[1], mMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = add(wOut[2], mMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = add(wOut[3], mMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = add(wOut[4], mMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = add(wOut[5], mMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = add(wOut[6], mMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = add(wOut[7], mMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = add(vOutM[i], mMul(v0M[i], v1M[i], qv, inv), qv)
	}
}

// MMulSubVecTo computes vOut -= v0 * v1 mod q in Montgomery form.
func MMulSubVecTo(v0M, v1M []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = sub(wOut[0], mMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = sub(wOut[1], mMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = sub(wOut[2], mMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = sub(wOut[3], mMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = sub(wOut[4], mMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = sub(wOut[5], mMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = sub(wOut[6], mMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = sub(wOut[7], mMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = sub(vOutM[i], mMul(v0M[i], v1M[i], qv, inv), qv)
	}
}

// MMulLazyVec returns v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazyVec(v0M, v1M []uint64, q *Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	MMulLazyVecTo(v0M, v1M, q, vOutM)
	return vOutM
}

// MMulLazyVecTo computes vOut = v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazyVecTo(v0M, v1M []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] = mMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] = mMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] = mMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] = mMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] = mMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] = mMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] = mMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] = mMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = mMulLazy(v0M[i], v1M[i], qv, inv)
	}
}

// MMulAddLazyVecTo computes vOut += v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulAddLazyVecTo(v0M, v1M []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] += mMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] += mMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] += mMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] += mMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] += mMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] += mMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] += mMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] += mMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += mMulLazy(v0M[i], v1M[i], qv, inv)
	}
}

// MMulSubLazyVecTo computes vOut -= v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulSubLazyVecTo(v0M, v1M []uint64, q *Modulus, vOutM []uint64) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))

		wOut[0] += mMulLazy(qv-w0[0], w1[0], qv, inv)
		wOut[1] += mMulLazy(qv-w0[1], w1[1], qv, inv)
		wOut[2] += mMulLazy(qv-w0[2], w1[2], qv, inv)
		wOut[3] += mMulLazy(qv-w0[3], w1[3], qv, inv)

		wOut[4] += mMulLazy(qv-w0[4], w1[4], qv, inv)
		wOut[5] += mMulLazy(qv-w0[5], w1[5], qv, inv)
		wOut[6] += mMulLazy(qv-w0[6], w1[6], qv, inv)
		wOut[7] += mMulLazy(qv-w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += mMulLazy(qv-v0M[i], v1M[i], qv, inv)
	}
}

// ReduceVec returns v mod q.
func ReduceVec(v []uint64, q *Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ReduceVecTo(v, q, vOut)
	return vOut
}

// ReduceVecTo computes vOut = v0 mod q.
func ReduceVecTo(v0 []uint64, q *Modulus, vOut []uint64) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[1] = bMod64(w0[1], qv, divHi)
		wOut[0] = bMod64(w0[0], qv, divHi)
		wOut[2] = bMod64(w0[2], qv, divHi)
		wOut[3] = bMod64(w0[3], qv, divHi)

		wOut[4] = bMod64(w0[4], qv, divHi)
		wOut[5] = bMod64(w0[5], qv, divHi)
		wOut[6] = bMod64(w0[6], qv, divHi)
		wOut[7] = bMod64(w0[7], qv, divHi)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = bMod64(v0[i], qv, divHi)
	}
}
