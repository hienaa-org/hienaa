//go:build amd64 && !purego

package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/num"
	"golang.org/x/sys/cpu"
)

// AddTo computes vOut = v0 + v1 mod q.
func AddTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		addToAVX2(v0, v1, q.Value(), vOut)
		return
	}

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
	if cpu.X86.HasAVX2 {
		addLazyToAVX2(v0, v1, vOut)
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

// SubTo computes vOut = v0 - v1 mod q.
func SubTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		subToAVX2(v0, v1, q.Value(), vOut)
		return
	}

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
	if cpu.X86.HasAVX2 {
		subLazyToAVX2(v0, v1, vOut)
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

// BMulTo computes vOut = v0 * v1 mod q using Barrett reduction.
func BMulTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	divHi, divLo := q.Div()
	bMulToX86(v0, v1, q.Value(), divHi, divLo, vOut)
}

// BMulAddTo computes vOut += v0 * v1 mod q using Barrett reduction.
func BMulAddTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	divHi, divLo := q.Div()
	bMulAddToX86(v0, v1, q.Value(), divHi, divLo, vOut)
}

// BMulSubTo computes vOut -= v0 * v1 mod q using Barrett reduction.
func BMulSubTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	divHi, divLo := q.Div()
	bMulSubToX86(v0, v1, q.Value(), divHi, divLo, vOut)
}

// BMulLazyTo computes vOut = v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	divHi, divLo := q.Div()
	bMulLazyToX86(v0, v1, q.Value(), divHi, divLo, vOut)
}

// BMulAddLazyTo computes vOut += v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulAddLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	divHi, divLo := q.Div()
	bMulAddLazyToX86(v0, v1, q.Value(), divHi, divLo, vOut)
}

// BMulSubLazyTo computes vOut -= v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulSubLazyTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	divHi, divLo := q.Div()
	bMulSubLazyToX86(v0, v1, q.Value(), divHi, divLo, vOut)
}
