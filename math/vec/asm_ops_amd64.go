//go:build amd64 && !purego

package vec

import (
	"github.com/hienaa-org/hienaa/math/num"
	"golang.org/x/sys/cpu"
)

// AddTo computes vOut = v0 + v1 mod q.
func AddTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		addToAVX2(v0, v1, q.Value(), vOut)
		return
	}

	N := len(vOut)
	M := (N >> 3) << 3

	var w0, w1, wOut []uint64
	for i := 0; i < M; i += 8 {
		w0 = v0[i : i+8 : i+8]
		w1 = v1[i : i+8 : i+8]
		wOut = vOut[i : i+8 : i+8]

		wOut[0] = num.Add(w0[0], w1[0], q)
		wOut[1] = num.Add(w0[1], w1[1], q)
		wOut[2] = num.Add(w0[2], w1[2], q)
		wOut[3] = num.Add(w0[3], w1[3], q)

		wOut[4] = num.Add(w0[4], w1[4], q)
		wOut[5] = num.Add(w0[5], w1[5], q)
		wOut[6] = num.Add(w0[6], w1[6], q)
		wOut[7] = num.Add(w0[7], w1[7], q)
	}

	for i := M; i < N; i++ {
		vOut[i] = num.Add(v0[i], v1[i], q)
	}
}

// AddLazyTo computes vOut = v0 + v1.
func AddLazyTo(v0, v1, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		addLazyToAVX2(v0, v1, vOut)
		return
	}

	N := len(vOut)
	M := (N >> 3) << 3

	var w0, w1, wOut []uint64
	for i := 0; i < M; i += 8 {
		w0 = v0[i : i+8 : i+8]
		w1 = v1[i : i+8 : i+8]
		wOut = vOut[i : i+8 : i+8]

		wOut[0] = w0[0] + w1[0]
		wOut[1] = w0[1] + w1[1]
		wOut[2] = w0[2] + w1[2]
		wOut[3] = w0[3] + w1[3]

		wOut[4] = w0[4] + w1[4]
		wOut[5] = w0[5] + w1[5]
		wOut[6] = w0[6] + w1[6]
		wOut[7] = w0[7] + w1[7]
	}

	for i := M; i < N; i++ {
		vOut[i] = v0[i] + v1[i]
	}
}

// SubTo computes vOut = v0 - v1 mod q.
func SubTo(v0, v1 []uint64, q *num.Modulus, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		subToAVX2(v0, v1, q.Value(), vOut)
		return
	}

	N := len(vOut)
	M := (N >> 3) << 3

	var w0, w1, wOut []uint64
	for i := 0; i < M; i += 8 {
		w0 = v0[i : i+8 : i+8]
		w1 = v1[i : i+8 : i+8]
		wOut = vOut[i : i+8 : i+8]

		wOut[0] = num.Sub(w0[0], w1[0], q)
		wOut[1] = num.Sub(w0[1], w1[1], q)
		wOut[2] = num.Sub(w0[2], w1[2], q)
		wOut[3] = num.Sub(w0[3], w1[3], q)

		wOut[4] = num.Sub(w0[4], w1[4], q)
		wOut[5] = num.Sub(w0[5], w1[5], q)
		wOut[6] = num.Sub(w0[6], w1[6], q)
		wOut[7] = num.Sub(w0[7], w1[7], q)
	}

	for i := M; i < N; i++ {
		vOut[i] = num.Sub(v0[i], v1[i], q)
	}
}

// SubLazyTo computes vOut = v0 - v1.
func SubLazyTo(v0, v1, vOut []uint64) {
	if cpu.X86.HasAVX2 {
		subLazyToAVX2(v0, v1, vOut)
		return
	}

	N := len(vOut)
	M := (N >> 3) << 3

	var w0, w1, wOut []uint64
	for i := 0; i < M; i += 8 {
		w0 = v0[i : i+8 : i+8]
		w1 = v1[i : i+8 : i+8]
		wOut = vOut[i : i+8 : i+8]

		wOut[0] = w0[0] - w1[0]
		wOut[1] = w0[1] - w1[1]
		wOut[2] = w0[2] - w1[2]
		wOut[3] = w0[3] - w1[3]

		wOut[4] = w0[4] - w1[4]
		wOut[5] = w0[5] - w1[5]
		wOut[6] = w0[6] - w1[6]
		wOut[7] = w0[7] - w1[7]
	}

	for i := M; i < N; i++ {
		vOut[i] = v0[i] - v1[i]
	}
}
