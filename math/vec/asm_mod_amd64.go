//go:build amd64 && !purego

package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
	"golang.org/x/sys/cpu"
)

// AddTo computes vOut = v0 + v1 mod q.
func AddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		addToAVX512(vOut, v0, v1, q.Value())
		return
	case cpu.X86.HasAVX2:
		addToAVX2(vOut, v0, v1, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

		wOut[0] = modops.Add(w0[0], w1[0], qv)
		wOut[1] = modops.Add(w0[1], w1[1], qv)
		wOut[2] = modops.Add(w0[2], w1[2], qv)
		wOut[3] = modops.Add(w0[3], w1[3], qv)

		wOut[4] = modops.Add(w0[4], w1[4], qv)
		wOut[5] = modops.Add(w0[5], w1[5], qv)
		wOut[6] = modops.Add(w0[6], w1[6], qv)
		wOut[7] = modops.Add(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Add(v0[i], v1[i], qv)
	}
}

// AddLazyTo computes vOut = v0 + v1.
func AddLazyTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		addLazyToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX2:
		addLazyToAVX2(vOut, v0, v1)
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

// SubTo computes vOut = v0 - v1 mod q.
func SubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		subToAVX512(vOut, v0, v1, q.Value())
		return
	case cpu.X86.HasAVX2:
		subToAVX2(vOut, v0, v1, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

		wOut[0] = modops.Sub(w0[0], w1[0], qv)
		wOut[1] = modops.Sub(w0[1], w1[1], qv)
		wOut[2] = modops.Sub(w0[2], w1[2], qv)
		wOut[3] = modops.Sub(w0[3], w1[3], qv)

		wOut[4] = modops.Sub(w0[4], w1[4], qv)
		wOut[5] = modops.Sub(w0[5], w1[5], qv)
		wOut[6] = modops.Sub(w0[6], w1[6], qv)
		wOut[7] = modops.Sub(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Sub(v0[i], v1[i], qv)
	}
}

// SubLazyTo computes vOut = v0 - v1.
func SubLazyTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		subLazyToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX2:
		subLazyToAVX2(vOut, v0, v1)
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
