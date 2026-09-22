package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
)

// vecScalarOpOut outputs vOut for vec-scalar operations.
func vecScalarOpOut[T0, T1 uint64 | []uint64](v0 T0, v1 T1) []uint64 {
	v0Vec, v0IsVec := any(v0).([]uint64)
	v1Vec, v1IsVec := any(v1).([]uint64)

	if !v0IsVec && !v1IsVec {
		panic("inconsistent input(s)")
	}

	if v0IsVec {
		return make([]uint64, len(v0Vec))
	}
	return make([]uint64, len(v1Vec))
}

// Add returns vOut = v0 + v1 mod q.
// If q is nil, then it returns vOut = v0 + v1.
//
// Panics when v0, v1 are scalars or have different lengths.
func Add[T0, T1 uint64 | []uint64](v0 T0, v1 T1, q *num.Modulus) []uint64 {
	vOut := vecScalarOpOut(v0, v1)
	AddTo(vOut, v0, v1, q)
	return vOut
}

// Sub returns vOut = v0 - v1 mod q.
// If q is nil, then it returns vOut = v0 - v1.
//
// Panics when v0, v1 are scalars or have different lengths.
func Sub[T0, T1 uint64 | []uint64](v0 T0, v1 T1, q *num.Modulus) []uint64 {
	vOut := vecScalarOpOut(v0, v1)
	SubTo(vOut, v0, v1, q)
	return vOut
}

// Neg returns -v mod q.
// If q is nil, then it returns -v.
func Neg(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	NegTo(vOut, v, q)
	return vOut
}

// NegTo computes -v mod q.
// If q is nil, then it returns -v.
func NegTo(vOut, v []uint64, q *num.Modulus) {
	if q != nil {
		negTo(vOut, v, q)
		return
	}
	negWordTo(vOut, v)
}

// Mul returns vOut = v0 * v1 mod q.
// If q is nil, then it returns vOut = v0 * v1.
//
// Panics when v0, v1 are scalars or have different lengths.
func Mul[T0, T1 uint64 | []uint64](v0 T0, v1 T1, q *num.Modulus) []uint64 {
	vOut := vecScalarOpOut(v0, v1)
	MulTo(vOut, v0, v1, q)
	return vOut
}

// MulTo computes vOut = v0 * v1 mod q.
// If q is nil, then it computes vOut = v0 * v1.
//
// Panics when v0, v1 are scalars or have different lengths.
func MulTo[T0, T1 uint64 | []uint64](vOut []uint64, v0 T0, v1 T1, q *num.Modulus) {
	switch v0 := any(v0).(type) {
	case uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			if q == nil {
				mulScalarWordTo(vOut, v1, v0)
				return
			}
			v0 = num.Reduce(v0, q)
			sMulScalarTo(vOut, v1, v0, num.SForm(v0, q), q)
		}

	case []uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			if q == nil {
				mulScalarWordTo(vOut, v0, v1)
				return
			}
			v1 = num.Reduce(v1, q)
			sMulScalarTo(vOut, v0, v1, num.SForm(v1, q), q)

		case []uint64:
			if q == nil {
				mulWordTo(vOut, v0, v1)
				return
			}
			mulTo(vOut, v0, v1, q)
		}
	}
}

// mulTo computes vOut = v0 * v1 mod q using Barrett reduction.
func mulTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	checkLength(len(vOut), len(v0), len(v1))

	qv := q.Value()
	divHi, divLo := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r0 := unsafe.Pointer(unsafe.SliceData(v0))
	r1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w0 := (*[8]uint64)(unsafe.Add(r0, uintptr(i)*L))
		w1 := (*[8]uint64)(unsafe.Add(r1, uintptr(i)*L))

		wOut[0] = modops.BMul(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] = modops.BMul(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] = modops.BMul(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] = modops.BMul(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] = modops.BMul(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] = modops.BMul(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] = modops.BMul(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] = modops.BMul(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.BMul(v0[i], v1[i], qv, divHi, divLo)
	}
}

// MulAddTo computes vOut += v0 * v1 mod q.
// If q is nil, then it computes vOut += v0 * v1.
//
// Panics when v0, v1 are scalars or have different lengths.
func MulAddTo[T0, T1 uint64 | []uint64](vOut []uint64, v0 T0, v1 T1, q *num.Modulus) {
	switch v0 := any(v0).(type) {
	case uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			if q == nil {
				mulAddScalarWordTo(vOut, v1, v0)
				return
			}
			v0 = num.Reduce(v0, q)
			sMulAddScalarTo(vOut, v1, v0, num.SForm(v0, q), q)
		}

	case []uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			if q == nil {
				mulAddScalarWordTo(vOut, v0, v1)
				return
			}
			v1 = num.Reduce(v1, q)
			sMulAddScalarTo(vOut, v0, v1, num.SForm(v1, q), q)

		case []uint64:
			if q == nil {
				mulAddWordTo(vOut, v0, v1)
				return
			}
			mulAddTo(vOut, v0, v1, q)
		}
	}
}

// mulAddTo computes vOut += v0 * v1 mod q using Barrett reduction.
func mulAddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	checkLength(len(vOut), len(v0), len(v1))

	qv := q.Value()
	divHi, divLo := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r0 := unsafe.Pointer(unsafe.SliceData(v0))
	r1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w0 := (*[8]uint64)(unsafe.Add(r0, uintptr(i)*L))
		w1 := (*[8]uint64)(unsafe.Add(r1, uintptr(i)*L))

		wOut[0] = modops.Add(wOut[0], modops.BMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = modops.Add(wOut[1], modops.BMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = modops.Add(wOut[2], modops.BMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = modops.Add(wOut[3], modops.BMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = modops.Add(wOut[4], modops.BMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = modops.Add(wOut[5], modops.BMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = modops.Add(wOut[6], modops.BMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = modops.Add(wOut[7], modops.BMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Add(vOut[i], modops.BMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
}

// MulSubTo computes vOut -= v0 * v1 mod q.
// If q is nil, then it computes vOut -= v0 * v1.
//
// Panics when v0, v1 are scalars or have different lengths.
func MulSubTo[T0, T1 uint64 | []uint64](vOut []uint64, v0 T0, v1 T1, q *num.Modulus) {
	switch v0 := any(v0).(type) {
	case uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			if q == nil {
				mulSubScalarWordTo(vOut, v1, v0)
				return
			}
			v0 = num.Reduce(v0, q)
			sMulSubScalarTo(vOut, v1, v0, num.SForm(v0, q), q)
		}

	case []uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			if q == nil {
				mulSubScalarWordTo(vOut, v0, v1)
				return
			}
			v1 = num.Reduce(v1, q)
			sMulSubScalarTo(vOut, v0, v1, num.SForm(v1, q), q)

		case []uint64:
			if q == nil {
				mulSubWordTo(vOut, v0, v1)
				return
			}
			mulSubTo(vOut, v0, v1, q)
		}
	}
}

// mulSubTo computes vOut -= v0 * v1 mod q using Barrett reduction.
func mulSubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	checkLength(len(vOut), len(v0), len(v1))

	qv := q.Value()
	divHi, divLo := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r0 := unsafe.Pointer(unsafe.SliceData(v0))
	r1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w0 := (*[8]uint64)(unsafe.Add(r0, uintptr(i)*L))
		w1 := (*[8]uint64)(unsafe.Add(r1, uintptr(i)*L))

		wOut[0] = modops.Sub(wOut[0], modops.BMul(w0[0], w1[0], qv, divHi, divLo), qv)
		wOut[1] = modops.Sub(wOut[1], modops.BMul(w0[1], w1[1], qv, divHi, divLo), qv)
		wOut[2] = modops.Sub(wOut[2], modops.BMul(w0[2], w1[2], qv, divHi, divLo), qv)
		wOut[3] = modops.Sub(wOut[3], modops.BMul(w0[3], w1[3], qv, divHi, divLo), qv)

		wOut[4] = modops.Sub(wOut[4], modops.BMul(w0[4], w1[4], qv, divHi, divLo), qv)
		wOut[5] = modops.Sub(wOut[5], modops.BMul(w0[5], w1[5], qv, divHi, divLo), qv)
		wOut[6] = modops.Sub(wOut[6], modops.BMul(w0[6], w1[6], qv, divHi, divLo), qv)
		wOut[7] = modops.Sub(wOut[7], modops.BMul(w0[7], w1[7], qv, divHi, divLo), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Sub(vOut[i], modops.BMul(v0[i], v1[i], qv, divHi, divLo), qv)
	}
}

// MulLazy returns vOut = v0 * v1 mod q,
// but the result is in [0, 2q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MulLazy[T0, T1 uint64 | []uint64](v0 T0, v1 T1, q *num.Modulus) []uint64 {
	vOut := vecScalarOpOut(v0, v1)
	MulLazyTo(vOut, v0, v1, q)
	return vOut
}

// MulLazyTo computes vOut = v0 * v1 mod q,
// but the result is in [0, 2q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MulLazyTo[T0, T1 uint64 | []uint64](vOut []uint64, v0 T0, v1 T1, q *num.Modulus) {
	switch v0 := any(v0).(type) {
	case uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			v0 = num.Reduce(v0, q)
			sMulScalarLazyTo(vOut, v1, v0, num.SForm(v0, q), q)
		}

	case []uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			v1 = num.Reduce(v1, q)
			sMulScalarLazyTo(vOut, v0, v1, num.SForm(v1, q), q)

		case []uint64:
			mulLazyTo(vOut, v0, v1, q)
		}
	}
}

// mulLazyTo computes vOut = v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
func mulLazyTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	checkLength(len(vOut), len(v0), len(v1))

	qv := q.Value()
	divHi, divLo := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r0 := unsafe.Pointer(unsafe.SliceData(v0))
	r1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w0 := (*[8]uint64)(unsafe.Add(r0, uintptr(i)*L))
		w1 := (*[8]uint64)(unsafe.Add(r1, uintptr(i)*L))

		wOut[0] = modops.BMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] = modops.BMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] = modops.BMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] = modops.BMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] = modops.BMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] = modops.BMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] = modops.BMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] = modops.BMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.BMulLazy(v0[i], v1[i], qv, divHi, divLo)
	}
}

// MulAddLazyTo computes vOut += v0 * v1 mod q,
// but the result is in [0, 3q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MulAddLazyTo[T0, T1 uint64 | []uint64](vOut []uint64, v0 T0, v1 T1, q *num.Modulus) {
	switch v0 := any(v0).(type) {
	case uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			v0 = num.Reduce(v0, q)
			sMulAddScalarLazyTo(vOut, v1, v0, num.SForm(v0, q), q)
		}

	case []uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			v1 = num.Reduce(v1, q)
			sMulAddScalarLazyTo(vOut, v0, v1, num.SForm(v1, q), q)

		case []uint64:
			mulAddLazyTo(vOut, v0, v1, q)
		}
	}
}

// mulAddLazyTo computes vOut += v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
func mulAddLazyTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	checkLength(len(vOut), len(v0), len(v1))

	qv := q.Value()
	divHi, divLo := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r0 := unsafe.Pointer(unsafe.SliceData(v0))
	r1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w0 := (*[8]uint64)(unsafe.Add(r0, uintptr(i)*L))
		w1 := (*[8]uint64)(unsafe.Add(r1, uintptr(i)*L))

		wOut[0] += modops.BMulLazy(w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += modops.BMulLazy(w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += modops.BMulLazy(w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += modops.BMulLazy(w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += modops.BMulLazy(w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += modops.BMulLazy(w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += modops.BMulLazy(w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += modops.BMulLazy(w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.BMulLazy(v0[i], v1[i], qv, divHi, divLo)
	}
}

// MulSubLazyTo computes vOut -= v0 * v1 mod q,
// but the result is in [0, 3q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MulSubLazyTo[T0, T1 uint64 | []uint64](vOut []uint64, v0 T0, v1 T1, q *num.Modulus) {
	switch v0 := any(v0).(type) {
	case uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			v0 = num.Reduce(v0, q)
			sMulSubScalarLazyTo(vOut, v1, v0, num.SForm(v0, q), q)
		}

	case []uint64:
		switch v1 := any(v1).(type) {
		case uint64:
			v1 = num.Reduce(v1, q)
			sMulSubScalarLazyTo(vOut, v0, v1, num.SForm(v1, q), q)

		case []uint64:
			mulSubLazyTo(vOut, v0, v1, q)
		}
	}
}

// mulSubLazyTo computes vOut -= v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
func mulSubLazyTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	checkLength(len(vOut), len(v0), len(v1))

	qv := q.Value()
	divHi, divLo := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r0 := unsafe.Pointer(unsafe.SliceData(v0))
	r1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w0 := (*[8]uint64)(unsafe.Add(r0, uintptr(i)*L))
		w1 := (*[8]uint64)(unsafe.Add(r1, uintptr(i)*L))

		wOut[0] += modops.BMulLazy(qv-w0[0], w1[0], qv, divHi, divLo)
		wOut[1] += modops.BMulLazy(qv-w0[1], w1[1], qv, divHi, divLo)
		wOut[2] += modops.BMulLazy(qv-w0[2], w1[2], qv, divHi, divLo)
		wOut[3] += modops.BMulLazy(qv-w0[3], w1[3], qv, divHi, divLo)

		wOut[4] += modops.BMulLazy(qv-w0[4], w1[4], qv, divHi, divLo)
		wOut[5] += modops.BMulLazy(qv-w0[5], w1[5], qv, divHi, divLo)
		wOut[6] += modops.BMulLazy(qv-w0[6], w1[6], qv, divHi, divLo)
		wOut[7] += modops.BMulLazy(qv-w0[7], w1[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.BMulLazy(qv-v0[i], v1[i], qv, divHi, divLo)
	}
}

// MForm returns v in Montgomery form.
//
// Panics if q is even or nil.
func MForm(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MFormTo(vOut, v, q)
	return vOut
}

// InvMForm transforms vM to Normal form.
//
// Panics if q is even or nil.
func InvMForm(vM []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(vM))
	InvMFormTo(vOut, vM, q)
	return vOut
}

// MMul returns vOut = v0 * v1 mod q using Montgomery multiplication.
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMul[T0, T1 uint64 | []uint64](v0M T0, v1M T1, q *num.Modulus) []uint64 {
	vOutM := vecScalarOpOut(v0M, v1M)
	MMulTo(vOutM, v0M, v1M, q)
	return vOutM
}

// MMulTo computes vOut = v0 * v1 mod q using Montgomery multiplication.
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulTo[T0, T1 uint64 | []uint64](vOutM []uint64, v0M T0, v1M T1, q *num.Modulus) {
	switch v0 := any(v0M).(type) {
	case uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			mMulScalarTo(vOutM, v1, v0, q)
		}

	case []uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			mMulScalarTo(vOutM, v0, v1, q)

		case []uint64:
			mMulTo(vOutM, v0, v1, q)
		}
	}
}

// MMulAddTo computes vOut += v0 * v1 mod q using Montgomery multiplication.
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulAddTo[T0, T1 uint64 | []uint64](vOutM []uint64, v0M T0, v1M T1, q *num.Modulus) {
	switch v0 := any(v0M).(type) {
	case uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			mMulAddScalarTo(vOutM, v1, v0, q)
		}

	case []uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			mMulAddScalarTo(vOutM, v0, v1, q)

		case []uint64:
			mMulAddTo(vOutM, v0, v1, q)
		}
	}
}

// MMulSubTo computes vOut -= v0 * v1 mod q using Montgomery multiplication.
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulSubTo[T0, T1 uint64 | []uint64](vOutM []uint64, v0M T0, v1M T1, q *num.Modulus) {
	switch v0 := any(v0M).(type) {
	case uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			mMulSubScalarTo(vOutM, v1, v0, q)
		}

	case []uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			mMulSubScalarTo(vOutM, v0, v1, q)

		case []uint64:
			mMulSubTo(vOutM, v0, v1, q)
		}
	}
}

// MMulLazy returns vOut = v0 * v1 mod q using Montgomery multiplication,
// but the result is in [0, 2q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulLazy[T0, T1 uint64 | []uint64](v0M T0, v1M T1, q *num.Modulus) []uint64 {
	vOutM := vecScalarOpOut(v0M, v1M)
	MMulLazyTo(vOutM, v0M, v1M, q)
	return vOutM
}

// MMulLazyTo computes vOut = v0 * v1 mod q using Montgomery multiplication,
// but the result is in [0, 2q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulLazyTo[T0, T1 uint64 | []uint64](vOutM []uint64, v0M T0, v1M T1, q *num.Modulus) {
	switch v0 := any(v0M).(type) {
	case uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			mMulScalarLazyTo(vOutM, v1, v0, q)
		}

	case []uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			mMulScalarLazyTo(vOutM, v0, v1, q)

		case []uint64:
			mMulLazyTo(vOutM, v0, v1, q)
		}
	}
}

// MMulAddLazyTo computes vOut += v0 * v1 mod q using Montgomery multiplication,
// but the result is in [0, 3q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulAddLazyTo[T0, T1 uint64 | []uint64](vOutM []uint64, v0M T0, v1M T1, q *num.Modulus) {
	switch v0 := any(v0M).(type) {
	case uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			mMulAddScalarLazyTo(vOutM, v1, v0, q)
		}

	case []uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			mMulAddScalarLazyTo(vOutM, v0, v1, q)

		case []uint64:
			mMulAddLazyTo(vOutM, v0, v1, q)
		}
	}
}

// MMulSubLazyTo computes vOut -= v0 * v1 mod q using Montgomery multiplication,
// but the result is in [0, 3q).
//
// Panics when v0, v1 are scalars or have different lengths, or q is nil.
func MMulSubLazyTo[T0, T1 uint64 | []uint64](vOutM []uint64, v0M T0, v1M T1, q *num.Modulus) {
	switch v0 := any(v0M).(type) {
	case uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			panic("inconsistent input(s)")

		case []uint64:
			mMulSubScalarLazyTo(vOutM, v1, v0, q)
		}

	case []uint64:
		switch v1 := any(v1M).(type) {
		case uint64:
			mMulSubScalarLazyTo(vOutM, v0, v1, q)

		case []uint64:
			mMulSubLazyTo(vOutM, v0, v1, q)
		}
	}
}

// SForm returns v in Shoup form.
//
// Panics if q is nil.
func SForm(v []uint64, q *num.Modulus) []uint64 {
	vOutS := make([]uint64, len(v))
	SFormTo(vOutS, v, q)
	return vOutS
}

// SFormTo transforms v to Shoup form to vOutS.
//
// Panics if q is nil.
func SFormTo(vOutS, v []uint64, q *num.Modulus) {
	checkLength(len(vOutS), len(v))

	qv := q.Value()

	M := (len(vOutS) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOutS))
	r := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(i)*L))

		wOut[0] = modops.SForm(w[0], qv)
		wOut[1] = modops.SForm(w[1], qv)
		wOut[2] = modops.SForm(w[2], qv)
		wOut[3] = modops.SForm(w[3], qv)

		wOut[4] = modops.SForm(w[4], qv)
		wOut[5] = modops.SForm(w[5], qv)
		wOut[6] = modops.SForm(w[6], qv)
		wOut[7] = modops.SForm(w[7], qv)
	}

	for i := M; i < len(vOutS); i++ {
		vOutS[i] = modops.SForm(v[i], qv)
	}
}

// SMul returns vOut = v0 * v1 mod q using Shoup multiplication.
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMul[T1 uint64 | []uint64](v0 []uint64, v1, v1S T1, q *num.Modulus) []uint64 {
	vOut := vecScalarOpOut(v0, v1)
	SMulTo(vOut, v0, v1, v1S, q)
	return vOut
}

// SMulTo computes vOut = v0 * v1 mod q using Shoup multiplication.
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulTo[T1 uint64 | []uint64](vOut, v0 []uint64, v1, v1S T1, q *num.Modulus) {
	switch v1 := any(v1).(type) {
	case uint64:
		switch v1S := any(v1S).(type) {
		case uint64:
			sMulScalarTo(vOut, v0, v1, v1S, q)
		}

	case []uint64:
		switch v1S := any(v1S).(type) {
		case []uint64:
			sMulTo(vOut, v0, v1, v1S, q)
		}
	}
}

// SMulAddTo computes vOut += v0 * v1 mod q using Shoup multiplication.
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulAddTo[T1 uint64 | []uint64](vOut, v0 []uint64, v1, v1S T1, q *num.Modulus) {
	switch v1 := any(v1).(type) {
	case uint64:
		switch v1S := any(v1S).(type) {
		case uint64:
			sMulAddScalarTo(vOut, v0, v1, v1S, q)
		}

	case []uint64:
		switch v1S := any(v1S).(type) {
		case []uint64:
			sMulAddTo(vOut, v0, v1, v1S, q)
		}
	}
}

// SMulSubTo computes vOut -= v0 * v1 mod q using Shoup multiplication.
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulSubTo[T1 uint64 | []uint64](vOut, v0 []uint64, v1, v1S T1, q *num.Modulus) {
	switch v1 := any(v1).(type) {
	case uint64:
		switch v1S := any(v1S).(type) {
		case uint64:
			sMulSubScalarTo(vOut, v0, v1, v1S, q)
		}

	case []uint64:
		switch v1S := any(v1S).(type) {
		case []uint64:
			sMulSubTo(vOut, v0, v1, v1S, q)
		}
	}
}

// SMulLazy returns vOut = v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulLazy[T1 uint64 | []uint64](v0 []uint64, v1, v1S T1, q *num.Modulus) []uint64 {
	vOut := vecScalarOpOut(v0, v1)
	SMulLazyTo(vOut, v0, v1, v1S, q)
	return vOut
}

// SMulLazyTo computes vOut = v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulLazyTo[T1 uint64 | []uint64](vOut, v0 []uint64, v1, v1S T1, q *num.Modulus) {
	switch v1 := any(v1).(type) {
	case uint64:
		switch v1S := any(v1S).(type) {
		case uint64:
			sMulScalarLazyTo(vOut, v0, v1, v1S, q)
		}

	case []uint64:
		switch v1S := any(v1S).(type) {
		case []uint64:
			sMulLazyTo(vOut, v0, v1, v1S, q)
		}
	}
}

// SMulAddLazyTo computes vOut += v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulAddLazyTo[T1 uint64 | []uint64](vOut, v0 []uint64, v1, v1S T1, q *num.Modulus) {
	switch v1 := any(v1).(type) {
	case uint64:
		switch v1S := any(v1S).(type) {
		case uint64:
			sMulAddScalarLazyTo(vOut, v0, v1, v1S, q)
		}

	case []uint64:
		switch v1S := any(v1S).(type) {
		case []uint64:
			sMulAddLazyTo(vOut, v0, v1, v1S, q)
		}
	}
}

// SMulSubLazyTo computes vOut -= v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics when v0, v1 have different lengths, or q is nil.
func SMulSubLazyTo[T1 uint64 | []uint64](vOut, v0 []uint64, v1, v1S T1, q *num.Modulus) {
	switch v1 := any(v1).(type) {
	case uint64:
		switch v1S := any(v1S).(type) {
		case uint64:
			sMulSubScalarLazyTo(vOut, v0, v1, v1S, q)
		}

	case []uint64:
		switch v1S := any(v1S).(type) {
		case []uint64:
			sMulSubLazyTo(vOut, v0, v1, v1S, q)
		}
	}
}

// Reduce returns vOut = v mod q.
//
// Panics if q is nil.
func Reduce[T num.Integer](v []T, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ReduceTo(vOut, v, q)
	return vOut
}

// Reduce2Q returns vOut = v mod q assuming v is in [0, 2q).
//
// Panics if q is nil.
func Reduce2Q(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	Reduce2QTo(vOut, v, q)
	return vOut
}

// Reduce4Q returns vOut = v mod q assuming v is in [0, 4q).
//
// Panics if q is nil.
func Reduce4Q(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	Reduce4QTo(vOut, v, q)
	return vOut
}
