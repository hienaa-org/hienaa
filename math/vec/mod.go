package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
)

// MulForm stores precomputed values for modular multiplication.
type MulForm struct {
	// Float is the float64 representation of the value.
	Float []float64
	// SForm is the Shoup Form of the value.
	SForm []uint64
}

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

// MulScalar returns v * c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v * c.
func MulScalar(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MulScalarTo(vOut, v, c, q)
	return vOut
}

// MulScalarTo computes vOut = v * c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v * c.
func MulScalarTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		mulScalarTo(vOut, v, c, q)
		return
	}
	mulScalarWordTo(vOut, v, c)
}

// MulAddScalarTo computes vOut += v * c mod q.
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

// FMulScalar returns v * c mod q using [num.MulForm] of c.
//
// Panics if q is nil.
func FMulScalar(v []uint64, c uint64, cM num.MulForm, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	FMulScalarTo(vOut, v, c, cM, q)
	return vOut
}

// Mul returns v0 * v1 mod q.
// v0 and v1 must be in [0, q).
func Mul(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	MulTo(vOut, v0, v1, q)
	return vOut
}

// MulTo computes vOut = v0 * v1 mod q.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns v0 * v1.
func MulTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulTo(vOut, v0, v1, q)
		return
	}
	mulWordTo(vOut, v0, v1)
}

// MulAddTo computes vOut += v0 * v1 mod q.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns vOut += v0 * v1.
func MulAddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulAddTo(vOut, v0, v1, q)
		return
	}
	mulAddWordTo(vOut, v0, v1)
}

// MulSubTo computes vOut -= v0 * v1 mod q.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns vOut -= v0 * v1.
func MulSubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		mulSubTo(vOut, v0, v1, q)
		return
	}
	mulSubWordTo(vOut, v0, v1)
}

// ToMulForm transforms v into [MulForm].
//
// Panics if q is nil
func ToMulForm(v []uint64, q *num.Modulus) MulForm {
	vF := make([]float64, len(v))
	vS := make([]uint64, len(v))

	qv := q.Value()

	M := (len(v) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rF := unsafe.Pointer(unsafe.SliceData(vF))
	rS := unsafe.Pointer(unsafe.SliceData(vS))
	r := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wF := (*[8]float64)(unsafe.Add(rF, uintptr(i)*L))
		wS := (*[8]uint64)(unsafe.Add(rS, uintptr(i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(i)*L))

		wF[0] = float64(w[0])
		wF[1] = float64(w[1])
		wF[2] = float64(w[2])
		wF[3] = float64(w[3])

		wF[4] = float64(w[4])
		wF[5] = float64(w[5])
		wF[6] = float64(w[6])
		wF[7] = float64(w[7])

		wS[0] = modops.SForm(w[0], qv)
		wS[1] = modops.SForm(w[1], qv)
		wS[2] = modops.SForm(w[2], qv)
		wS[3] = modops.SForm(w[3], qv)

		wS[4] = modops.SForm(w[4], qv)
		wS[5] = modops.SForm(w[5], qv)
		wS[6] = modops.SForm(w[6], qv)
		wS[7] = modops.SForm(w[7], qv)
	}

	for i := M; i < len(v); i++ {
		vF[i] = float64(v[i])
		vS[i] = modops.SForm(v[i], qv)
	}

	return MulForm{
		Float: vF,
		SForm: vS,
	}
}

// FMul returns v0 * v1 mod q using [MulForm] of v1.
//
// Panics if q is nil.
func FMul(v0, v1 []uint64, v1M MulForm, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	FMulTo(vOut, v0, v1, v1M, q)
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
