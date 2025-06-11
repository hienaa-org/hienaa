//go:build !(amd64 && !purego)

package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/internal/mod"
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

		wOut[0] = mod.Add(w0[0], w1[0], qv)
		wOut[1] = mod.Add(w0[1], w1[1], qv)
		wOut[2] = mod.Add(w0[2], w1[2], qv)
		wOut[3] = mod.Add(w0[3], w1[3], qv)

		wOut[4] = mod.Add(w0[4], w1[4], qv)
		wOut[5] = mod.Add(w0[5], w1[5], qv)
		wOut[6] = mod.Add(w0[6], w1[6], qv)
		wOut[7] = mod.Add(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.Add(v0[i], v1[i], qv)
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

		wOut[0] = mod.Sub(w0[0], w1[0], qv)
		wOut[1] = mod.Sub(w0[1], w1[1], qv)
		wOut[2] = mod.Sub(w0[2], w1[2], qv)
		wOut[3] = mod.Sub(w0[3], w1[3], qv)

		wOut[4] = mod.Sub(w0[4], w1[4], qv)
		wOut[5] = mod.Sub(w0[5], w1[5], qv)
		wOut[6] = mod.Sub(w0[6], w1[6], qv)
		wOut[7] = mod.Sub(w0[7], w1[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.Sub(v0[i], v1[i], qv)
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

		wOut[0] = mod.BMul(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] = mod.BMul(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] = mod.BMul(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] = mod.BMul(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] = mod.BMul(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] = mod.BMul(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] = mod.BMul(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] = mod.BMul(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.BMul(v0[i], v1[i], qv, divHi, divLo)
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

		wOut[0] = mod.Add(wOut[0], mod.BMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = mod.Add(wOut[1], mod.BMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = mod.Add(wOut[2], mod.BMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = mod.Add(wOut[3], mod.BMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = mod.Add(wOut[4], mod.BMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = mod.Add(wOut[5], mod.BMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = mod.Add(wOut[6], mod.BMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = mod.Add(wOut[7], mod.BMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.Add(vOut[i], mod.BMul(v0[i], v1[i], qv, divHi, divLo), qv)
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

		wOut[0] = mod.Sub(wOut[0], mod.BMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = mod.Sub(wOut[1], mod.BMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = mod.Sub(wOut[2], mod.BMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = mod.Sub(wOut[3], mod.BMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = mod.Sub(wOut[4], mod.BMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = mod.Sub(wOut[5], mod.BMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = mod.Sub(wOut[6], mod.BMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = mod.Sub(wOut[7], mod.BMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.Sub(vOut[i], mod.BMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
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

		wOut[0] = mod.BMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] = mod.BMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] = mod.BMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] = mod.BMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] = mod.BMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] = mod.BMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] = mod.BMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] = mod.BMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = mod.BMulLazy(v0[i], v1[i], qv, divHi, divLo)
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

		wOut[0] += mod.BMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += mod.BMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += mod.BMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += mod.BMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += mod.BMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += mod.BMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += mod.BMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += mod.BMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += mod.BMulLazy(v0[i], v1[i], qv, divHi, divLo)
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

		wOut[0] += qv<<1 - mod.BMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += qv<<1 - mod.BMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += qv<<1 - mod.BMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += qv<<1 - mod.BMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += qv<<1 - mod.BMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += qv<<1 - mod.BMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += qv<<1 - mod.BMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += qv<<1 - mod.BMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += qv<<1 - mod.BMulLazy(v0[i], v1[i], qv, divHi, divLo)
	}
}
