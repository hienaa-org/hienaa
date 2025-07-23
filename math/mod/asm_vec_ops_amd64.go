//go:build amd64 && !purego

package mod

import (
	"unsafe"

	"golang.org/x/sys/cpu"
)

// AddVecTo computes vOut = v0 + v1 mod q.
func AddVecTo(vOut, v0, v1 []uint64, q *Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		addVecToAVX512(vOut, v0, v1, q.Value())
		return
	case cpu.X86.HasAVX2:
		addVecToAVX2(vOut, v0, v1, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
func AddLazyVecTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		addLazyVecToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX2:
		addLazyVecToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
func SubVecTo(vOut, v0, v1 []uint64, q *Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		subVecToAVX512(vOut, v0, v1, q.Value())
		return
	case cpu.X86.HasAVX2:
		subVecToAVX2(vOut, v0, v1, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
func SubLazyVecTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		subLazyVecToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX2:
		subLazyVecToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
