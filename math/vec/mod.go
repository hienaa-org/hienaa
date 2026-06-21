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

// MulScalar returns v * c mod q using Barrett or Float reduction.
// v and c must be in [0, q).
// If q is nil, then it returns v * c.
func MulScalar(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MulScalarTo(vOut, v, c, q)
	return vOut
}

// MulScalarTo computes vOut = v * c mod q using Barrett or Float reduction.
// v and c must be in [0, q).
// If q is nil, then it returns v * c.
func MulScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		mulScalarTo(vOut, v, c, q)
		return
	}
	mulScalarWordTo(vOut, v, c)
}

// MulAddScalarTo computes vOut += v * c mod q using Barrett or Float reduction.
// v and c must be in [0, q).
// If q is nil, then it returns vOut += v * c.
func MulAddScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		mulAddScalarTo(vOut, v, c, q)
		return
	}
	mulAddScalarWordTo(vOut, v, c)
}

// MulSubScalarTo computes vOut -= v * c mod q using Barrett or Float multiplication.
// v and c must be in [0, q).
// If q is nil, then it returns vOut -= v * c.
func MulSubScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		mulSubScalarTo(vOut, v, c, q)
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

// Mul returns v0 * v1 mod q using Barrett or Float reduction.
// v0 and v1 must be in [0, q).
func Mul(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	MulTo(vOut, v0, v1, q)
	return vOut
}

// MulTo computes vOut = v0 * v1 mod q using Barrett or Float reduction.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns v0 * v1.
func MulTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulTo(vOut, v0, v1, q)
		return
	}
	mulWordTo(vOut, v0, v1)
}

// MulAddTo computes vOut += v0 * v1 mod q using Barrett or Float reduction.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns vOut += v0 * v1.
func MulAddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulAddTo(vOut, v0, v1, q)
		return
	}
	mulAddWordTo(vOut, v0, v1)
}

// MulSubTo computes vOut -= v0 * v1 mod q using Barrett or Float reduction.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns vOut -= v0 * v1.
func MulSubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulSubTo(vOut, v0, v1, q)
		return
	}
	mulSubWordTo(vOut, v0, v1)
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
