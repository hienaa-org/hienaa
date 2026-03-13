package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
)

// Add returns v0 + v1 mod q.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns v0 + v1.
func Add(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	AddTo(vOut, v0, v1, q)
	return vOut
}

// AddScalar returns v + c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v + c.
func AddScalar(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	AddScalarTo(vOut, v, c, q)
	return vOut
}

// Sub returns v0 - v1 mod q.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns v0 - v1.
func Sub(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SubTo(vOut, v0, v1, q)
	return vOut
}

// SubScalar returns v - c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v - c.
func SubScalar(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	SubScalarTo(vOut, v, c, q)
	return vOut
}

// Neg returns -v mod q.
// v must be in [0, q).
// If q is nil, then it returns -v.
func Neg(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	NegTo(vOut, v, q)
	return vOut
}

// NegTo computes -v mod q.
// v must be in [0, q).
// If q is nil, then it returns -v.
func NegTo(vOut, v []uint64, q *num.Modulus) {
	if q != nil {
		negTo(vOut, v, q)
		return
	}
	negWordTo(vOut, v)
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

// MulScalar returns v * c mod q using Shoup multiplication.
//
// Panics if q is nil.
func MulScalar(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MulScalarTo(vOut, v, c, q)
	return vOut
}

// MulScalarTo computes vOut = v * c mod q using Shoup multiplication.
// If q is nil, then it returns x0 * x1.
func MulScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		divHi, _ := q.Div()
		SMulScalarTo(vOut, v, c, modops.SForm(c, q.Value(), divHi), q)
		return
	}
	mulScalarWordTo(vOut, v, c)
}

// MulAddScalarTo computes vOut += v * c mod q using Shoup multiplication.
// If q is nil, then it returns vOut += v * c.
func MulAddScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		divHi, _ := q.Div()
		SMulAddScalarTo(vOut, v, c, modops.SForm(c, q.Value(), divHi), q)
		return
	}
	mulAddScalarWordTo(vOut, v, c)
}

// MulSubScalarTo computes vOut -= v * c mod q using Shoup multiplication.
// If q is nil, then it returns vOut -= v * c.
func MulSubScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		divHi, _ := q.Div()
		SMulSubScalarTo(vOut, v, c, modops.SForm(c, q.Value(), divHi), q)
		return
	}
	mulSubScalarWordTo(vOut, v, c)
}

// SMulScalar returns v * c mod q using Shoup multiplication.
//
// Panics if q is nil.
func SMulScalar(v []uint64, c, cS uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	SMulScalarTo(vOut, v, c, cS, q)
	return vOut
}

// MulScalarLazy returns v * c mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func MulScalarLazy(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MulScalarLazyTo(vOut, v, c, q)
	return vOut
}

// MulScalarLazyTo computes vOut = c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func MulScalarLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	divHi, _ := q.Div()
	SMulScalarLazyTo(vOut, v, c, modops.SForm(c, q.Value(), divHi), q)
}

// MulAddScalarLazyTo computes vOut += c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func MulAddScalarLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	divHi, _ := q.Div()
	SMulAddScalarLazyTo(vOut, v, c, modops.SForm(c, q.Value(), divHi), q)
}

// MulSubScalarLazyTo computes vOut -= c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func MulSubScalarLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	divHi, _ := q.Div()
	SMulSubScalarLazyTo(vOut, v, c, modops.SForm(c, q.Value(), divHi), q)
}

// SMulScalarLazy returns v * c mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func SMulScalarLazy(v []uint64, c, cS uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	SMulScalarLazyTo(vOut, v, c, cS, q)
	return vOut
}

// MMulScalar returns v * c mod q using Montgomery multiplication.
// When c is in Montgomery form, the output is the same form as v.
//
// Panics if q is nil.
func MMulScalar(vM []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(vM))
	MMulScalarTo(vOutM, vM, cM, q)
	return vOutM
}

// Mul returns v0 * v1 mod q using Barrett reduction.
func Mul(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	MulTo(vOut, v0, v1, q)
	return vOut
}

// MulTo computes vOut = v0 * v1 mod q using Barrett reduction.
// If q is nil, then it returns v0 * v1.
func MulTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulTo(vOut, v0, v1, q)
		return
	}
	mulWordTo(vOut, v0, v1)
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

// MulAddTo computes vOut += v0 * v1 mod q using Barrett reduction.
// If q is nil, then it returns vOut += v0 * v1.
func MulAddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulAddTo(vOut, v0, v1, q)
		return
	}
	mulAddWordTo(vOut, v0, v1)
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

// MulSubTo computes vOut -= v0 * v1 mod q using Barrett reduction.
// If q is nil, then it returns vOut -= v0 * v1.
func MulSubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulSubTo(vOut, v0, v1, q)
		return
	}
	mulSubWordTo(vOut, v0, v1)
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

// MulLazyTo computes vOut = v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func MulLazyTo(vOut, v0, v1 []uint64, q *num.Modulus) {
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

// MulAddLazyTo computes vOut += v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func MulAddLazyTo(vOut, v0, v1 []uint64, q *num.Modulus) {
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

// MulSubLazyTo computes vOut -= v0 * v1 mod q using Barrett reduction,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func MulSubLazyTo(vOut, v0, v1 []uint64, q *num.Modulus) {
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

// MMul returns v0 * v1 mod q in Montgomery form.
//
// Panics if q is nil.
func MMul(v0M, v1M []uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	MMulTo(vOutM, v0M, v1M, q)
	return vOutM
}

// MMulLazy returns v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func MMulLazy(v0M, v1M []uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	MMulLazyTo(vOutM, v0M, v1M, q)
	return vOutM
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
	divHi, _ := q.Div()

	M := (len(vOutS) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOutS))
	r := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(i)*L))

		wOut[0] = modops.SForm(w[0], qv, divHi)
		wOut[1] = modops.SForm(w[1], qv, divHi)
		wOut[2] = modops.SForm(w[2], qv, divHi)
		wOut[3] = modops.SForm(w[3], qv, divHi)

		wOut[4] = modops.SForm(w[4], qv, divHi)
		wOut[5] = modops.SForm(w[5], qv, divHi)
		wOut[6] = modops.SForm(w[6], qv, divHi)
		wOut[7] = modops.SForm(w[7], qv, divHi)
	}

	for i := M; i < len(vOutS); i++ {
		vOutS[i] = modops.SForm(v[i], qv, divHi)
	}
}

// SMul returns v0 * v1 mod q using Shoup multiplication.
//
// Panics if q is nil.
func SMul(v0, v1, v1S []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SMulTo(vOut, v0, v1, v1S, q)
	return vOut
}

// Reduce returns v mod q.
//
// Panics if q is nil.
func Reduce[T num.Integer](v []T, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ReduceTo(vOut, v, q)
	return vOut
}

// ReduceTo computes vOut = v mod q.
//
// Panics if q is nil.
func ReduceTo[T num.Integer](vOut []uint64, v []T, q *num.Modulus) {
	checkLength(len(vOut), len(v))

	qv := q.Value()
	divHi, _ := q.Div()

	M := (len(vOut) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))
	LT := unsafe.Sizeof(T(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vOut))
	r := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))
		w := (*[8]T)(unsafe.Add(r, uintptr(i)*LT))

		wOut[0] = modops.BMod(w[0], qv, divHi)
		wOut[1] = modops.BMod(w[1], qv, divHi)
		wOut[2] = modops.BMod(w[2], qv, divHi)
		wOut[3] = modops.BMod(w[3], qv, divHi)

		wOut[4] = modops.BMod(w[4], qv, divHi)
		wOut[5] = modops.BMod(w[5], qv, divHi)
		wOut[6] = modops.BMod(w[6], qv, divHi)
		wOut[7] = modops.BMod(w[7], qv, divHi)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.BMod(v[i], qv, divHi)
	}
}
