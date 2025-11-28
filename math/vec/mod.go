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

// ScalarAdd returns v + c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v + c.
func ScalarAdd(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ScalarAddTo(vOut, v, c, q)
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

// ScalarSub returns v - c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v - c.
func ScalarSub(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ScalarSubTo(vOut, v, c, q)
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

// negTo computes vOut = -v mod q.
func negTo(vOut, v []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.Neg(w[0], qv)
		wOut[1] = modops.Neg(w[1], qv)
		wOut[2] = modops.Neg(w[2], qv)
		wOut[3] = modops.Neg(w[3], qv)

		wOut[4] = modops.Neg(w[4], qv)
		wOut[5] = modops.Neg(w[5], qv)
		wOut[6] = modops.Neg(w[6], qv)
		wOut[7] = modops.Neg(w[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Neg(v[i], qv)
	}
}

// negWordTo computes vOut = -v.
func negWordTo(vOut, v []uint64) {
	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = -w[0]
		wOut[1] = -w[1]
		wOut[2] = -w[2]
		wOut[3] = -w[3]

		wOut[4] = -w[4]
		wOut[5] = -w[5]
		wOut[6] = -w[6]
		wOut[7] = -w[7]
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = -v[i]
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

// MFormTo transforms v to Montgomery form to vOutM.
//
// Panics if q is even or nil.
func MFormTo(vOutM, v []uint64, q *num.Modulus) {
	if q.Inv() == 0 {
		panic("MFormTo: modulus is even")
	}

	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.MForm(w[0], qv, divHi, divLo)
		wOut[1] = modops.MForm(w[1], qv, divHi, divLo)
		wOut[2] = modops.MForm(w[2], qv, divHi, divLo)
		wOut[3] = modops.MForm(w[3], qv, divHi, divLo)

		wOut[4] = modops.MForm(w[4], qv, divHi, divLo)
		wOut[5] = modops.MForm(w[5], qv, divHi, divLo)
		wOut[6] = modops.MForm(w[6], qv, divHi, divLo)
		wOut[7] = modops.MForm(w[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MForm(v[i], qv, divHi, divLo)
	}
}

// InvMForm transforms vM to Normal form.
//
// Panics if q is even or nil.
func InvMForm(vM []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(vM))
	InvMFormTo(vOut, vM, q)
	return vOut
}

// InvMFormTo computes vOut as vM in Normal form.
//
// Panics if q is even or nil.
func InvMFormTo(vOut, vM []uint64, q *num.Modulus) {
	if q.Inv() == 0 {
		panic("InvMFormTo: modulus is even")
	}

	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.InvMForm(w[0], qv, inv)
		wOut[1] = modops.InvMForm(w[1], qv, inv)
		wOut[2] = modops.InvMForm(w[2], qv, inv)
		wOut[3] = modops.InvMForm(w[3], qv, inv)

		wOut[4] = modops.InvMForm(w[4], qv, inv)
		wOut[5] = modops.InvMForm(w[5], qv, inv)
		wOut[6] = modops.InvMForm(w[6], qv, inv)
		wOut[7] = modops.InvMForm(w[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.InvMForm(vM[i], qv, inv)
	}
}

// ScalarMul returns c * v mod q using Shoup multiplication.
//
// Panics if q is nil.
func ScalarMul(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ScalarMulTo(vOut, v, c, q)
	return vOut
}

// ScalarMulTo computes vOut = c * v mod q using Shoup multiplication.
// If q is nil, then it returns x0 * x1.
func ScalarMulTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		scalarMulTo(vOut, v, c, q)
		return
	}
	scalarMulWordTo(vOut, v, c)
}

// ScalarMulTo computes vOut = c * v mod q using Shoup multiplication.
func scalarMulTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.SMul(w[0], c, cS, qv)
		wOut[1] = modops.SMul(w[1], c, cS, qv)
		wOut[2] = modops.SMul(w[2], c, cS, qv)
		wOut[3] = modops.SMul(w[3], c, cS, qv)

		wOut[4] = modops.SMul(w[4], c, cS, qv)
		wOut[5] = modops.SMul(w[5], c, cS, qv)
		wOut[6] = modops.SMul(w[6], c, cS, qv)
		wOut[7] = modops.SMul(w[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.SMul(v[i], c, cS, qv)
	}
}

// ScalarMulLazyTo computes vOut = c * v.
func scalarMulWordTo(vOut, v []uint64, c uint64) {
	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = w[0] * c
		wOut[1] = w[1] * c
		wOut[2] = w[2] * c
		wOut[3] = w[3] * c

		wOut[4] = w[4] * c
		wOut[5] = w[5] * c
		wOut[6] = w[6] * c
		wOut[7] = w[7] * c
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = v[i] * c
	}
}

// ScalarMulAddTo computes vOut += c * v mod q using Shoup multiplication.
// If q is nil, then it returns vOut += c * v.
func ScalarMulAddTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		scalarMulAddTo(vOut, v, c, q)
		return
	}
	scalarMulAddWordTo(vOut, v, c)
}

// scalarMulAddTo computes vOut += c * v mod q using Shoup multiplication.
func scalarMulAddTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.Add(wOut[0], modops.SMul(w0[0], c, cS, qv), qv)
		wOut[1] = modops.Add(wOut[1], modops.SMul(w0[1], c, cS, qv), qv)
		wOut[2] = modops.Add(wOut[2], modops.SMul(w0[2], c, cS, qv), qv)
		wOut[3] = modops.Add(wOut[3], modops.SMul(w0[3], c, cS, qv), qv)

		wOut[4] = modops.Add(wOut[4], modops.SMul(w0[4], c, cS, qv), qv)
		wOut[5] = modops.Add(wOut[5], modops.SMul(w0[5], c, cS, qv), qv)
		wOut[6] = modops.Add(wOut[6], modops.SMul(w0[6], c, cS, qv), qv)
		wOut[7] = modops.Add(wOut[7], modops.SMul(w0[7], c, cS, qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Add(vOut[i], modops.SMul(v[i], c, cS, qv), qv)
	}
}

// scalarMulAddWordTo computes vOut += c * v.
func scalarMulAddWordTo(vOut, v []uint64, c uint64) {
	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] += w[0] * c
		wOut[1] += w[1] * c
		wOut[2] += w[2] * c
		wOut[3] += w[3] * c

		wOut[4] += w[4] * c
		wOut[5] += w[5] * c
		wOut[6] += w[6] * c
		wOut[7] += w[7] * c
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += v[i] * c
	}
}

// ScalarMulSubTo computes vOut -= c * v mod q using Shoup multiplication.
// If q is nil, then it returns vOut -= c * v.
func ScalarMulSubTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		scalarMulSubTo(vOut, v, c, q)
		return
	}
	scalarMulSubWordTo(vOut, v, c)
}

// scalarMulSubTo computes vOut -= c * v mod q using Shoup multiplication.
func scalarMulSubTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.Sub(wOut[0], modops.SMul(w[0], c, cS, qv), qv)
		wOut[1] = modops.Sub(wOut[1], modops.SMul(w[1], c, cS, qv), qv)
		wOut[2] = modops.Sub(wOut[2], modops.SMul(w[2], c, cS, qv), qv)
		wOut[3] = modops.Sub(wOut[3], modops.SMul(w[3], c, cS, qv), qv)

		wOut[4] = modops.Sub(wOut[4], modops.SMul(w[4], c, cS, qv), qv)
		wOut[5] = modops.Sub(wOut[5], modops.SMul(w[5], c, cS, qv), qv)
		wOut[6] = modops.Sub(wOut[6], modops.SMul(w[6], c, cS, qv), qv)
		wOut[7] = modops.Sub(wOut[7], modops.SMul(w[7], c, cS, qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Sub(vOut[i], modops.SMul(v[i], c, cS, qv), qv)
	}
}

// scalarMulSubWordTo computes vOut -= c * v.
func scalarMulSubWordTo(vOut, v []uint64, c uint64) {
	M := (len(vOut) >> 3) << 3

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] -= w[0] * c
		wOut[1] -= w[1] * c
		wOut[2] -= w[2] * c
		wOut[3] -= w[3] * c

		wOut[4] -= w[4] * c
		wOut[5] -= w[5] * c
		wOut[6] -= w[6] * c
		wOut[7] -= w[7] * c
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] -= v[i] * c
	}
}

// ScalarMulLazy returns c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func ScalarMulLazy(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ScalarMulLazyTo(vOut, v, c, q)
	return vOut
}

// ScalarMulLazyTo computes vOut = c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func ScalarMulLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.SMulLazy(w[0], c, cS, qv)
		wOut[1] = modops.SMulLazy(w[1], c, cS, qv)
		wOut[2] = modops.SMulLazy(w[2], c, cS, qv)
		wOut[3] = modops.SMulLazy(w[3], c, cS, qv)

		wOut[4] = modops.SMulLazy(w[4], c, cS, qv)
		wOut[5] = modops.SMulLazy(w[5], c, cS, qv)
		wOut[6] = modops.SMulLazy(w[6], c, cS, qv)
		wOut[7] = modops.SMulLazy(w[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.SMulLazy(v[i], c, cS, qv)
	}
}

// ScalarMulAddLazyTo computes vOut += c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func ScalarMulAddLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] += modops.SMulLazy(w[0], c, cS, qv)
		wOut[1] += modops.SMulLazy(w[1], c, cS, qv)
		wOut[2] += modops.SMulLazy(w[2], c, cS, qv)
		wOut[3] += modops.SMulLazy(w[3], c, cS, qv)

		wOut[4] += modops.SMulLazy(w[4], c, cS, qv)
		wOut[5] += modops.SMulLazy(w[5], c, cS, qv)
		wOut[6] += modops.SMulLazy(w[6], c, cS, qv)
		wOut[7] += modops.SMulLazy(w[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.SMulLazy(v[i], c, cS, qv)
	}
}

// ScalarMulSubLazyTo computes vOut -= c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func ScalarMulSubLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()
	cNeg := qv - c
	cNegS := modops.SForm(cNeg, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] += modops.SMulLazy(w[0], cNeg, cNegS, qv)
		wOut[1] += modops.SMulLazy(w[1], cNeg, cNegS, qv)
		wOut[2] += modops.SMulLazy(w[2], cNeg, cNegS, qv)
		wOut[3] += modops.SMulLazy(w[3], cNeg, cNegS, qv)

		wOut[4] += modops.SMulLazy(w[4], cNeg, cNegS, qv)
		wOut[5] += modops.SMulLazy(w[5], cNeg, cNegS, qv)
		wOut[6] += modops.SMulLazy(w[6], cNeg, cNegS, qv)
		wOut[7] += modops.SMulLazy(w[7], cNeg, cNegS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.SMulLazy(v[i], cNeg, cNegS, qv)
	}
}

// ScalarMMul returns c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, the output is the same form as v.
//
// Panics if q is nil.
func ScalarMMul(vM []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(vM))
	ScalarMMulTo(vOutM, vM, cM, q)
	return vOutM
}

// ScalarMMulTo computes vOut = c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.MMul(w[0], cM, qv, inv)
		wOut[1] = modops.MMul(w[1], cM, qv, inv)
		wOut[2] = modops.MMul(w[2], cM, qv, inv)
		wOut[3] = modops.MMul(w[3], cM, qv, inv)

		wOut[4] = modops.MMul(w[4], cM, qv, inv)
		wOut[5] = modops.MMul(w[5], cM, qv, inv)
		wOut[6] = modops.MMul(w[6], cM, qv, inv)
		wOut[7] = modops.MMul(w[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MMul(vM[i], cM, qv, inv)
	}
}

// ScalarMMulAddTo computes vOut += c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulAddTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.Add(wOut[0], modops.MMul(w[0], cM, qv, inv), qv)
		wOut[1] = modops.Add(wOut[1], modops.MMul(w[1], cM, qv, inv), qv)
		wOut[2] = modops.Add(wOut[2], modops.MMul(w[2], cM, qv, inv), qv)
		wOut[3] = modops.Add(wOut[3], modops.MMul(w[3], cM, qv, inv), qv)

		wOut[4] = modops.Add(wOut[4], modops.MMul(w[4], cM, qv, inv), qv)
		wOut[5] = modops.Add(wOut[5], modops.MMul(w[5], cM, qv, inv), qv)
		wOut[6] = modops.Add(wOut[6], modops.MMul(w[6], cM, qv, inv), qv)
		wOut[7] = modops.Add(wOut[7], modops.MMul(w[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.Add(vOutM[i], modops.MMul(vM[i], cM, qv, inv), qv)
	}
}

// ScalarMMulSubTo computes vOut -= c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulSubTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.Sub(wOut[0], modops.MMul(w[0], cM, qv, inv), qv)
		wOut[1] = modops.Sub(wOut[1], modops.MMul(w[1], cM, qv, inv), qv)
		wOut[2] = modops.Sub(wOut[2], modops.MMul(w[2], cM, qv, inv), qv)
		wOut[3] = modops.Sub(wOut[3], modops.MMul(w[3], cM, qv, inv), qv)

		wOut[4] = modops.Sub(wOut[4], modops.MMul(w[4], cM, qv, inv), qv)
		wOut[5] = modops.Sub(wOut[5], modops.MMul(w[5], cM, qv, inv), qv)
		wOut[6] = modops.Sub(wOut[6], modops.MMul(w[6], cM, qv, inv), qv)
		wOut[7] = modops.Sub(wOut[7], modops.MMul(w[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.Sub(vOutM[i], modops.MMul(vM[i], cM, qv, inv), qv)
	}
}

// ScalarMMulLazy returns c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, the output is the same form as v.
//
// Panics if q is nil.
func ScalarMMulLazy(vM []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(vM))
	ScalarMMulLazyTo(vOutM, vM, cM, q)
	return vOutM
}

// ScalarMMulLazyTo computes vOut = c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.MMulLazy(w[0], cM, qv, inv)
		wOut[1] = modops.MMulLazy(w[1], cM, qv, inv)
		wOut[2] = modops.MMulLazy(w[2], cM, qv, inv)
		wOut[3] = modops.MMulLazy(w[3], cM, qv, inv)

		wOut[4] = modops.MMulLazy(w[4], cM, qv, inv)
		wOut[5] = modops.MMulLazy(w[5], cM, qv, inv)
		wOut[6] = modops.MMulLazy(w[6], cM, qv, inv)
		wOut[7] = modops.MMulLazy(w[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MMulLazy(vM[i], cM, qv, inv)
	}
}

// ScalarMMulAddLazyTo computes vOut += c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulAddLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] += modops.MMulLazy(w[0], cM, qv, inv)
		wOut[1] += modops.MMulLazy(w[1], cM, qv, inv)
		wOut[2] += modops.MMulLazy(w[2], cM, qv, inv)
		wOut[3] += modops.MMulLazy(w[3], cM, qv, inv)

		wOut[4] += modops.MMulLazy(w[4], cM, qv, inv)
		wOut[5] += modops.MMulLazy(w[5], cM, qv, inv)
		wOut[6] += modops.MMulLazy(w[6], cM, qv, inv)
		wOut[7] += modops.MMulLazy(w[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += modops.MMulLazy(vM[i], cM, qv, inv)
	}
}

// ScalarMMulSubLazyTo computes vOut -= c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulSubLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	cMNeg := qv - cM

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] += modops.MMulLazy(w[0], cMNeg, qv, inv)
		wOut[1] += modops.MMulLazy(w[1], cMNeg, qv, inv)
		wOut[2] += modops.MMulLazy(w[2], cMNeg, qv, inv)
		wOut[3] += modops.MMulLazy(w[3], cMNeg, qv, inv)

		wOut[4] += modops.MMulLazy(w[4], cMNeg, qv, inv)
		wOut[5] += modops.MMulLazy(w[5], cMNeg, qv, inv)
		wOut[6] += modops.MMulLazy(w[6], cMNeg, qv, inv)
		wOut[7] += modops.MMulLazy(w[7], cMNeg, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += modops.MMulLazy(vM[i], cMNeg, qv, inv)
	}
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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))

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
	M := (len(vOutS) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutS[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.SForm(w0[0], qv, divHi)
		wOut[1] = modops.SForm(w0[1], qv, divHi)
		wOut[2] = modops.SForm(w0[2], qv, divHi)
		wOut[3] = modops.SForm(w0[3], qv, divHi)

		wOut[4] = modops.SForm(w0[4], qv, divHi)
		wOut[5] = modops.SForm(w0[5], qv, divHi)
		wOut[6] = modops.SForm(w0[6], qv, divHi)
		wOut[7] = modops.SForm(w0[7], qv, divHi)
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

// SMulTo computes vOut = v0 * v1 mod q using Shoup multiplication.
//
// Panics if q is nil.
func SMulTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		w1S := (*[8]uint64)(unsafe.Pointer(&v1S[i]))

		wOut[0] = modops.SMul(w0[0], w1[0], w1S[0], qv)
		wOut[1] = modops.SMul(w0[1], w1[1], w1S[1], qv)
		wOut[2] = modops.SMul(w0[2], w1[2], w1S[2], qv)
		wOut[3] = modops.SMul(w0[3], w1[3], w1S[3], qv)

		wOut[4] = modops.SMul(w0[4], w1[4], w1S[4], qv)
		wOut[5] = modops.SMul(w0[5], w1[5], w1S[5], qv)
		wOut[6] = modops.SMul(w0[6], w1[6], w1S[6], qv)
		wOut[7] = modops.SMul(w0[7], w1[7], w1S[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.SMul(v0[i], v1[i], v1S[i], qv)
	}
}

// SMulAddTo computes vOut += v0 * v1 mod q using Shoup multiplication.
//
// Panics if q is nil.
func SMulAddTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		w1S := (*[8]uint64)(unsafe.Pointer(&v1S[i]))

		wOut[0] = modops.Add(wOut[0], modops.SMul(w0[0], w1[0], w1S[0], qv), qv)
		wOut[1] = modops.Add(wOut[1], modops.SMul(w0[1], w1[1], w1S[1], qv), qv)
		wOut[2] = modops.Add(wOut[2], modops.SMul(w0[2], w1[2], w1S[2], qv), qv)
		wOut[3] = modops.Add(wOut[3], modops.SMul(w0[3], w1[3], w1S[3], qv), qv)

		wOut[4] = modops.Add(wOut[4], modops.SMul(w0[4], w1[4], w1S[4], qv), qv)
		wOut[5] = modops.Add(wOut[5], modops.SMul(w0[5], w1[5], w1S[5], qv), qv)
		wOut[6] = modops.Add(wOut[6], modops.SMul(w0[6], w1[6], w1S[6], qv), qv)
		wOut[7] = modops.Add(wOut[7], modops.SMul(w0[7], w1[7], w1S[7], qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Add(vOut[i], modops.SMul(v0[i], v1[i], v1S[i], qv), qv)
	}
}

// SMulSubTo computes vOut -= v0 * v1 mod q using Shoup multiplication.
//
// Panics if q is nil.
func SMulSubTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		w1S := (*[8]uint64)(unsafe.Pointer(&v1S[i]))

		wOut[0] = modops.Sub(wOut[0], modops.SMul(w0[0], w1[0], w1S[0], qv), qv)
		wOut[1] = modops.Sub(wOut[1], modops.SMul(w0[1], w1[1], w1S[1], qv), qv)
		wOut[2] = modops.Sub(wOut[2], modops.SMul(w0[2], w1[2], w1S[2], qv), qv)
		wOut[3] = modops.Sub(wOut[3], modops.SMul(w0[3], w1[3], w1S[3], qv), qv)

		wOut[4] = modops.Sub(wOut[4], modops.SMul(w0[4], w1[4], w1S[4], qv), qv)
		wOut[5] = modops.Sub(wOut[5], modops.SMul(w0[5], w1[5], w1S[5], qv), qv)
		wOut[6] = modops.Sub(wOut[6], modops.SMul(w0[6], w1[6], w1S[6], qv), qv)
		wOut[7] = modops.Sub(wOut[7], modops.SMul(w0[7], w1[7], w1S[7], qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Sub(vOut[i], modops.SMul(v0[i], v1[i], v1S[i], qv), qv)
	}
}

// SMulLazyTo computes vOut = v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func SMulLazyTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		w1S := (*[8]uint64)(unsafe.Pointer(&v1S[i]))

		wOut[0] = modops.SMulLazy(w0[0], w1[0], w1S[0], qv)
		wOut[1] = modops.SMulLazy(w0[1], w1[1], w1S[1], qv)
		wOut[2] = modops.SMulLazy(w0[2], w1[2], w1S[2], qv)
		wOut[3] = modops.SMulLazy(w0[3], w1[3], w1S[3], qv)

		wOut[4] = modops.SMulLazy(w0[4], w1[4], w1S[4], qv)
		wOut[5] = modops.SMulLazy(w0[5], w1[5], w1S[5], qv)
		wOut[6] = modops.SMulLazy(w0[6], w1[6], w1S[6], qv)
		wOut[7] = modops.SMulLazy(w0[7], w1[7], w1S[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.SMulLazy(v0[i], v1[i], v1S[i], qv)
	}
}

// SMulAddLazyTo computes vOut += v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func SMulAddLazyTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		w1S := (*[8]uint64)(unsafe.Pointer(&v1S[i]))

		wOut[0] += modops.SMulLazy(w0[0], w1[0], w1S[0], qv)
		wOut[1] += modops.SMulLazy(w0[1], w1[1], w1S[1], qv)
		wOut[2] += modops.SMulLazy(w0[2], w1[2], w1S[2], qv)
		wOut[3] += modops.SMulLazy(w0[3], w1[3], w1S[3], qv)

		wOut[4] += modops.SMulLazy(w0[4], w1[4], w1S[4], qv)
		wOut[5] += modops.SMulLazy(w0[5], w1[5], w1S[5], qv)
		wOut[6] += modops.SMulLazy(w0[6], w1[6], w1S[6], qv)
		wOut[7] += modops.SMulLazy(w0[7], w1[7], w1S[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.SMulLazy(v0[i], v1[i], v1S[i], qv)
	}
}

// SMulSubLazyTo computes vOut -= v0 * v1 mod q using Shoup multiplication,
// but the result is in [0, 3q).
//
// Panics if q is nil.
func SMulSubLazyTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1[i]))
		w1S := (*[8]uint64)(unsafe.Pointer(&v1S[i]))

		wOut[0] += modops.SMulLazy(qv-w0[0], w1[0], w1S[0], qv)
		wOut[1] += modops.SMulLazy(qv-w0[1], w1[1], w1S[1], qv)
		wOut[2] += modops.SMulLazy(qv-w0[2], w1[2], w1S[2], qv)
		wOut[3] += modops.SMulLazy(qv-w0[3], w1[3], w1S[3], qv)

		wOut[4] += modops.SMulLazy(qv-w0[4], w1[4], w1S[4], qv)
		wOut[5] += modops.SMulLazy(qv-w0[5], w1[5], w1S[5], qv)
		wOut[6] += modops.SMulLazy(qv-w0[6], w1[6], w1S[6], qv)
		wOut[7] += modops.SMulLazy(qv-w0[7], w1[7], w1S[7], qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.SMulLazy(qv-v0[i], v1[i], v1S[i], qv)
	}
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
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]T)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.BMod(w0[0], qv, divHi)
		wOut[1] = modops.BMod(w0[1], qv, divHi)
		wOut[2] = modops.BMod(w0[2], qv, divHi)
		wOut[3] = modops.BMod(w0[3], qv, divHi)

		wOut[4] = modops.BMod(w0[4], qv, divHi)
		wOut[5] = modops.BMod(w0[5], qv, divHi)
		wOut[6] = modops.BMod(w0[6], qv, divHi)
		wOut[7] = modops.BMod(w0[7], qv, divHi)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.BMod(v[i], qv, divHi)
	}
}
