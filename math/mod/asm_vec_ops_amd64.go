//go:build amd64 && !purego

package mod

import (
	"unsafe"

	"golang.org/x/sys/cpu"
)

// AddVecTo computes vOut = v0 + v1 mod q.
func AddVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		addVecToAVX2(v0, v1, q.Value(), vOut)
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = add(w0[0], w1[0], qv)
		wOut[1] = add(w0[1], w1[1], qv)
		wOut[2] = add(w0[2], w1[2], qv)
		wOut[3] = add(w0[3], w1[3], qv)

		wOut[4] = add(w0[4], w1[4], qv)
		wOut[5] = add(w0[5], w1[5], qv)
		wOut[6] = add(w0[6], w1[6], qv)
		wOut[7] = add(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = add(v0[i], v1[i], qv)
	}
}

// AddLazyVecTo computes vOut = v0 + v1.
func AddLazyVecTo(v0, v1, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		addLazyVecToAVX2(v0, v1, vOut)
		return
	}

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

// SubVecTo computes vOut = v0 - v1 mod q.
func SubVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		subVecToAVX2(v0, v1, q.Value(), vOut)
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = sub(w0[0], w1[0], qv)
		wOut[1] = sub(w0[1], w1[1], qv)
		wOut[2] = sub(w0[2], w1[2], qv)
		wOut[3] = sub(w0[3], w1[3], qv)

		wOut[4] = sub(w0[4], w1[4], qv)
		wOut[5] = sub(w0[5], w1[5], qv)
		wOut[6] = sub(w0[6], w1[6], qv)
		wOut[7] = sub(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = sub(v0[i], v1[i], qv)
	}
}

// SubLazyVecTo computes vOut = v0 - v1.
func SubLazyVecTo(v0, v1, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		subLazyVecToAVX2(v0, v1, vOut)
		return
	}

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

// MulVecTo computes vOut = v0 * v1 mod q using Barrett reduction.
func MulVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasBMI2 {
		divHi, divLo := q.Div()
		mulVecToX86(v0, v1, q.Value(), divHi, divLo, vOut)
		return
	}

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

// MulAddVecTo computes vOut += v0 * v1 mod q using Barrett reduction.
func MulAddVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasBMI2 {
		divHi, divLo := q.Div()
		mulAddVecToX86(v0, v1, q.Value(), divHi, divLo, vOut)
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = add(wOut[0], bMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = add(wOut[1], bMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = add(wOut[2], bMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = add(wOut[3], bMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = add(wOut[4], bMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = add(wOut[5], bMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = add(wOut[6], bMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = add(wOut[7], bMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = add(vOut[i], bMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
}

// MulSubVecTo computes vOut -= v0 * v1 mod q using Barrett reduction.
func MulSubVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasBMI2 {
		divHi, divLo := q.Div()
		mulSubVecToX86(v0, v1, q.Value(), divHi, divLo, vOut)
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] = sub(wOut[0], bMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = sub(wOut[1], bMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = sub(wOut[2], bMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = sub(wOut[3], bMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = sub(wOut[4], bMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = sub(wOut[5], bMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = sub(wOut[6], bMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = sub(wOut[7], bMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = sub(vOut[i], bMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
}

// MulLazyVecTo computes vOut = v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func MulLazyVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasBMI2 {
		divHi, divLo := q.Div()
		mulLazyVecToX86(v0, v1, q.Value(), divHi, divLo, vOut)
		return
	}

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

// MulAddLazyVecTo computes vOut += v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
func MulAddLazyVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasBMI2 {
		divHi, divLo := q.Div()
		mulAddLazyVecToX86(v0, v1, q.Value(), divHi, divLo, vOut)
		return
	}

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

// MulSubLazyVecTo computes vOut -= v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
func MulSubLazyVecTo(v0, v1 []uint64, q *Modulus, vOut []uint64) {
	if cpu.X86.HasBMI2 {
		divHi, divLo := q.Div()
		mulSubLazyVecToX86(v0, v1, q.Value(), divHi, divLo, vOut)
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[0] += bMulLazy(qv-w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += bMulLazy(qv-w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += bMulLazy(qv-w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += bMulLazy(qv-w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += bMulLazy(qv-w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += bMulLazy(qv-w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += bMulLazy(qv-w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += bMulLazy(qv-w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += bMulLazy(qv-v0[i], v1[i], qv, divHi, divLo)
	}
}
