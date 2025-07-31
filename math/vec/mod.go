package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
)

// Add returns v0 + v1 mod q.
func Add(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	AddTo(vOut, v0, v1, q)
	return vOut
}

// AddLazy returns v0 + v1.
func AddLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	AddLazyTo(vOut, v0, v1)
	return vOut
}

// Sub returns v0 - v1 mod q.
func Sub(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SubTo(vOut, v0, v1, q)
	return vOut
}

// SubLazy returns v0 - v1.
func SubLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	SubLazyTo(vOut, v0, v1)
	return vOut
}

// MForm returns v in Montgomery form.
func MForm(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	MFormTo(vOut, v, q)
	return vOut
}

// MFormTo transforms v to Montgomery form to vOutM.
func MFormTo(vOutM, v []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.MForm(w0[0], qv, divHi, divLo)
		wOut[1] = modops.MForm(w0[1], qv, divHi, divLo)
		wOut[2] = modops.MForm(w0[2], qv, divHi, divLo)
		wOut[3] = modops.MForm(w0[3], qv, divHi, divLo)

		wOut[4] = modops.MForm(w0[4], qv, divHi, divLo)
		wOut[5] = modops.MForm(w0[5], qv, divHi, divLo)
		wOut[6] = modops.MForm(w0[6], qv, divHi, divLo)
		wOut[7] = modops.MForm(w0[7], qv, divHi, divLo)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MForm(v[i], qv, divHi, divLo)
	}
}

// InvMForm transforms vM to Normal form.
func InvMForm(vM []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(vM))
	InvMFormTo(vOut, vM, q)
	return vOut
}

// InvMFormTo computes vOut as vM in Normal form.
func InvMFormTo(vOut, vM []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.InvMForm(w0[0], qv, inv)
		wOut[1] = modops.InvMForm(w0[1], qv, inv)
		wOut[2] = modops.InvMForm(w0[2], qv, inv)
		wOut[3] = modops.InvMForm(w0[3], qv, inv)

		wOut[4] = modops.InvMForm(w0[4], qv, inv)
		wOut[5] = modops.InvMForm(w0[5], qv, inv)
		wOut[6] = modops.InvMForm(w0[6], qv, inv)
		wOut[7] = modops.InvMForm(w0[7], qv, inv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.InvMForm(vM[i], qv, inv)
	}
}

// ScalarMul returns c * v mod q using Shoup multiplication.
func ScalarMul(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ScalarMulTo(vOut, v, c, q)
	return vOut
}

// ScalarMulTo computes vOut = c * v mod q using Shoup multiplication.
func ScalarMulTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := modops.SForm(c, qv)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.SMul(w0[0], c, cS, qv)
		wOut[1] = modops.SMul(w0[1], c, cS, qv)
		wOut[2] = modops.SMul(w0[2], c, cS, qv)
		wOut[3] = modops.SMul(w0[3], c, cS, qv)

		wOut[4] = modops.SMul(w0[4], c, cS, qv)
		wOut[5] = modops.SMul(w0[5], c, cS, qv)
		wOut[6] = modops.SMul(w0[6], c, cS, qv)
		wOut[7] = modops.SMul(w0[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.SMul(v[i], c, cS, qv)
	}
}

// ScalarMulAddTo computes vOut += c * v mod q using Shoup multiplication.
func ScalarMulAddTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := modops.SForm(c, qv)

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

// ScalarMulSubTo computes vOut -= c * v mod q using Shoup multiplication.
func ScalarMulSubTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := modops.SForm(c, qv)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.Sub(wOut[0], modops.SMul(w0[0], c, cS, qv), qv)
		wOut[1] = modops.Sub(wOut[1], modops.SMul(w0[1], c, cS, qv), qv)
		wOut[2] = modops.Sub(wOut[2], modops.SMul(w0[2], c, cS, qv), qv)
		wOut[3] = modops.Sub(wOut[3], modops.SMul(w0[3], c, cS, qv), qv)

		wOut[4] = modops.Sub(wOut[4], modops.SMul(w0[4], c, cS, qv), qv)
		wOut[5] = modops.Sub(wOut[5], modops.SMul(w0[5], c, cS, qv), qv)
		wOut[6] = modops.Sub(wOut[6], modops.SMul(w0[6], c, cS, qv), qv)
		wOut[7] = modops.Sub(wOut[7], modops.SMul(w0[7], c, cS, qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Sub(vOut[i], modops.SMul(v[i], c, cS, qv), qv)
	}
}

// ScalarMulLazy returns c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
func ScalarMulLazy(v []uint64, c uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ScalarMulLazyTo(vOut, v, c, q)
	return vOut
}

// ScalarMulLazyTo computes vOut = c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
func ScalarMulLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := modops.SForm(c, qv)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.SMulLazy(w0[0], c, cS, qv)
		wOut[1] = modops.SMulLazy(w0[1], c, cS, qv)
		wOut[2] = modops.SMulLazy(w0[2], c, cS, qv)
		wOut[3] = modops.SMulLazy(w0[3], c, cS, qv)

		wOut[4] = modops.SMulLazy(w0[4], c, cS, qv)
		wOut[5] = modops.SMulLazy(w0[5], c, cS, qv)
		wOut[6] = modops.SMulLazy(w0[6], c, cS, qv)
		wOut[7] = modops.SMulLazy(w0[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.SMulLazy(v[i], c, cS, qv)
	}
}

// ScalarMulAddTo computes vOut += c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
func ScalarMulAddLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cS := modops.SForm(c, qv)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] += modops.SMulLazy(w0[0], c, cS, qv)
		wOut[1] += modops.SMulLazy(w0[1], c, cS, qv)
		wOut[2] += modops.SMulLazy(w0[2], c, cS, qv)
		wOut[3] += modops.SMulLazy(w0[3], c, cS, qv)

		wOut[4] += modops.SMulLazy(w0[4], c, cS, qv)
		wOut[5] += modops.SMulLazy(w0[5], c, cS, qv)
		wOut[6] += modops.SMulLazy(w0[6], c, cS, qv)
		wOut[7] += modops.SMulLazy(w0[7], c, cS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.SMulLazy(v[i], c, cS, qv)
	}
}

// ScalarMulSubLazyTo computes vOut -= c * v mod q using Shoup multiplication,
// but the result is in [0, 3q).
func ScalarMulSubLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	cNeg := qv - c
	cNegS := modops.SForm(cNeg, qv)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] += modops.SMulLazy(w0[0], cNeg, cNegS, qv)
		wOut[1] += modops.SMulLazy(w0[1], cNeg, cNegS, qv)
		wOut[2] += modops.SMulLazy(w0[2], cNeg, cNegS, qv)
		wOut[3] += modops.SMulLazy(w0[3], cNeg, cNegS, qv)

		wOut[4] += modops.SMulLazy(w0[4], cNeg, cNegS, qv)
		wOut[5] += modops.SMulLazy(w0[5], cNeg, cNegS, qv)
		wOut[6] += modops.SMulLazy(w0[6], cNeg, cNegS, qv)
		wOut[7] += modops.SMulLazy(w0[7], cNeg, cNegS, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += modops.SMulLazy(v[i], cNeg, cNegS, qv)
	}
}

// ScalarMMul returns c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, the output is the same form as v.
func ScalarMMul(vM []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(vM))
	ScalarMMulTo(vOutM, vM, cM, q)
	return vOutM
}

// ScalarMMulTo computes vOut = c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.MMul(w0[0], cM, qv, inv)
		wOut[1] = modops.MMul(w0[1], cM, qv, inv)
		wOut[2] = modops.MMul(w0[2], cM, qv, inv)
		wOut[3] = modops.MMul(w0[3], cM, qv, inv)

		wOut[4] = modops.MMul(w0[4], cM, qv, inv)
		wOut[5] = modops.MMul(w0[5], cM, qv, inv)
		wOut[6] = modops.MMul(w0[6], cM, qv, inv)
		wOut[7] = modops.MMul(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MMul(vM[i], cM, qv, inv)
	}
}

// ScalarMMulAddTo computes vOut += c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulAddTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.Add(wOut[0], modops.MMul(w0[0], cM, qv, inv), qv)
		wOut[1] = modops.Add(wOut[1], modops.MMul(w0[1], cM, qv, inv), qv)
		wOut[2] = modops.Add(wOut[2], modops.MMul(w0[2], cM, qv, inv), qv)
		wOut[3] = modops.Add(wOut[3], modops.MMul(w0[3], cM, qv, inv), qv)

		wOut[4] = modops.Add(wOut[4], modops.MMul(w0[4], cM, qv, inv), qv)
		wOut[5] = modops.Add(wOut[5], modops.MMul(w0[5], cM, qv, inv), qv)
		wOut[6] = modops.Add(wOut[6], modops.MMul(w0[6], cM, qv, inv), qv)
		wOut[7] = modops.Add(wOut[7], modops.MMul(w0[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.Add(vOutM[i], modops.MMul(vM[i], cM, qv, inv), qv)
	}
}

// ScalarMMulSubTo computes vOut -= c * v mod q using Montgomery multiplication.
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulSubTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.Sub(wOut[0], modops.MMul(w0[0], cM, qv, inv), qv)
		wOut[1] = modops.Sub(wOut[1], modops.MMul(w0[1], cM, qv, inv), qv)
		wOut[2] = modops.Sub(wOut[2], modops.MMul(w0[2], cM, qv, inv), qv)
		wOut[3] = modops.Sub(wOut[3], modops.MMul(w0[3], cM, qv, inv), qv)

		wOut[4] = modops.Sub(wOut[4], modops.MMul(w0[4], cM, qv, inv), qv)
		wOut[5] = modops.Sub(wOut[5], modops.MMul(w0[5], cM, qv, inv), qv)
		wOut[6] = modops.Sub(wOut[6], modops.MMul(w0[6], cM, qv, inv), qv)
		wOut[7] = modops.Sub(wOut[7], modops.MMul(w0[7], cM, qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.Sub(vOutM[i], modops.MMul(vM[i], cM, qv, inv), qv)
	}
}

// ScalarMMulLazy returns c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, the output is the same form as v.
func ScalarMMulLazy(vM []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(vM))
	ScalarMMulLazyTo(vOutM, vM, cM, q)
	return vOutM
}

// ScalarMMulLazyTo computes vOut = c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] = modops.MMulLazy(w0[0], cM, qv, inv)
		wOut[1] = modops.MMulLazy(w0[1], cM, qv, inv)
		wOut[2] = modops.MMulLazy(w0[2], cM, qv, inv)
		wOut[3] = modops.MMulLazy(w0[3], cM, qv, inv)

		wOut[4] = modops.MMulLazy(w0[4], cM, qv, inv)
		wOut[5] = modops.MMulLazy(w0[5], cM, qv, inv)
		wOut[6] = modops.MMulLazy(w0[6], cM, qv, inv)
		wOut[7] = modops.MMulLazy(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MMulLazy(vM[i], cM, qv, inv)
	}
}

// ScalarMMulAddLazyTo computes vOut += c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulAddLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] += modops.MMulLazy(w0[0], cM, qv, inv)
		wOut[1] += modops.MMulLazy(w0[1], cM, qv, inv)
		wOut[2] += modops.MMulLazy(w0[2], cM, qv, inv)
		wOut[3] += modops.MMulLazy(w0[3], cM, qv, inv)

		wOut[4] += modops.MMulLazy(w0[4], cM, qv, inv)
		wOut[5] += modops.MMulLazy(w0[5], cM, qv, inv)
		wOut[6] += modops.MMulLazy(w0[6], cM, qv, inv)
		wOut[7] += modops.MMulLazy(w0[7], cM, qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += modops.MMulLazy(vM[i], cM, qv, inv)
	}
}

// ScalarMMulSubLazyTo computes vOut -= c * v mod q using Montgomery multiplication,
// but the result is in [0, 3q).
// When c is in Motgomery form, vOut is the same form as v.
func ScalarMMulSubLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	cMNeg := qv - cM

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&vM[i]))

		wOut[0] += modops.MMulLazy(w0[0], cMNeg, qv, inv)
		wOut[1] += modops.MMulLazy(w0[1], cMNeg, qv, inv)
		wOut[2] += modops.MMulLazy(w0[2], cMNeg, qv, inv)
		wOut[3] += modops.MMulLazy(w0[3], cMNeg, qv, inv)

		wOut[4] += modops.MMulLazy(w0[4], cMNeg, qv, inv)
		wOut[5] += modops.MMulLazy(w0[5], cMNeg, qv, inv)
		wOut[6] += modops.MMulLazy(w0[6], cMNeg, qv, inv)
		wOut[7] += modops.MMulLazy(w0[7], cMNeg, qv, inv)
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
func MulTo(vOut, v0, v1 []uint64, q *num.Modulus) {
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
func MulAddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
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
func MulSubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
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
func MMul(v0M, v1M []uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	MMulTo(vOutM, v0M, v1M, q)
	return vOutM
}

// MMulTo computes vOut = v0 * v1 mod q in Montgomery form.
func MMulTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))

		wOut[0] = modops.MMul(w0[0], w1[0], qv, inv)
		wOut[1] = modops.MMul(w0[1], w1[1], qv, inv)
		wOut[2] = modops.MMul(w0[2], w1[2], qv, inv)
		wOut[3] = modops.MMul(w0[3], w1[3], qv, inv)

		wOut[4] = modops.MMul(w0[4], w1[4], qv, inv)
		wOut[5] = modops.MMul(w0[5], w1[5], qv, inv)
		wOut[6] = modops.MMul(w0[6], w1[6], qv, inv)
		wOut[7] = modops.MMul(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MMul(v0M[i], v1M[i], qv, inv)
	}
}

// MMulAddTo computes vOut += v0 * v1 mod q in Montgomery form.
func MMulAddTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))

		wOut[0] = modops.Add(wOut[0], modops.MMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = modops.Add(wOut[1], modops.MMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = modops.Add(wOut[2], modops.MMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = modops.Add(wOut[3], modops.MMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = modops.Add(wOut[4], modops.MMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = modops.Add(wOut[5], modops.MMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = modops.Add(wOut[6], modops.MMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = modops.Add(wOut[7], modops.MMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.Add(vOutM[i], modops.MMul(v0M[i], v1M[i], qv, inv), qv)
	}
}

// MMulSubTo computes vOut -= v0 * v1 mod q in Montgomery form.
func MMulSubTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))

		wOut[0] = modops.Sub(wOut[0], modops.MMul(w0[0], w1[0], qv, inv), qv)
		wOut[1] = modops.Sub(wOut[1], modops.MMul(w0[1], w1[1], qv, inv), qv)
		wOut[2] = modops.Sub(wOut[2], modops.MMul(w0[2], w1[2], qv, inv), qv)
		wOut[3] = modops.Sub(wOut[3], modops.MMul(w0[3], w1[3], qv, inv), qv)

		wOut[4] = modops.Sub(wOut[4], modops.MMul(w0[4], w1[4], qv, inv), qv)
		wOut[5] = modops.Sub(wOut[5], modops.MMul(w0[5], w1[5], qv, inv), qv)
		wOut[6] = modops.Sub(wOut[6], modops.MMul(w0[6], w1[6], qv, inv), qv)
		wOut[7] = modops.Sub(wOut[7], modops.MMul(w0[7], w1[7], qv, inv), qv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.Sub(vOutM[i], modops.MMul(v0M[i], v1M[i], qv, inv), qv)
	}
}

// MMulLazy returns v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(v0M, v1M []uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(v0M))
	MMulLazyTo(vOutM, v0M, v1M, q)
	return vOutM
}

// MMulLazyTo computes vOut = v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazyTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))

		wOut[0] = modops.MMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] = modops.MMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] = modops.MMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] = modops.MMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] = modops.MMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] = modops.MMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] = modops.MMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] = modops.MMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] = modops.MMulLazy(v0M[i], v1M[i], qv, inv)
	}
}

// MMulAddLazyTo computes vOut += v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulAddLazyTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))

		wOut[0] += modops.MMulLazy(w0[0], w1[0], qv, inv)
		wOut[1] += modops.MMulLazy(w0[1], w1[1], qv, inv)
		wOut[2] += modops.MMulLazy(w0[2], w1[2], qv, inv)
		wOut[3] += modops.MMulLazy(w0[3], w1[3], qv, inv)

		wOut[4] += modops.MMulLazy(w0[4], w1[4], qv, inv)
		wOut[5] += modops.MMulLazy(w0[5], w1[5], qv, inv)
		wOut[6] += modops.MMulLazy(w0[6], w1[6], qv, inv)
		wOut[7] += modops.MMulLazy(w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += modops.MMulLazy(v0M[i], v1M[i], qv, inv)
	}
}

// MMulSubLazyTo computes vOut -= v0 * v1 mod q in Montgomery form,
// but the result is in [0, 3q).
func MMulSubLazyTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	M := (len(vOutM) >> 3) << 3

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutM[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v0M[i]))
		w1 := (*[8]uint64)(unsafe.Pointer(&v1M[i]))

		wOut[0] += modops.MMulLazy(qv-w0[0], w1[0], qv, inv)
		wOut[1] += modops.MMulLazy(qv-w0[1], w1[1], qv, inv)
		wOut[2] += modops.MMulLazy(qv-w0[2], w1[2], qv, inv)
		wOut[3] += modops.MMulLazy(qv-w0[3], w1[3], qv, inv)

		wOut[4] += modops.MMulLazy(qv-w0[4], w1[4], qv, inv)
		wOut[5] += modops.MMulLazy(qv-w0[5], w1[5], qv, inv)
		wOut[6] += modops.MMulLazy(qv-w0[6], w1[6], qv, inv)
		wOut[7] += modops.MMulLazy(qv-w0[7], w1[7], qv, inv)
	}

	for i := M; i < len(vOutM); i++ {
		vOutM[i] += modops.MMulLazy(qv-v0M[i], v1M[i], qv, inv)
	}
}

// SForm returns v in Shoup form.
func SForm(v []uint64, q *num.Modulus) []uint64 {
	vOutS := make([]uint64, len(v))
	SFormTo(vOutS, v, q)
	return vOutS
}

// SFormTo transforms v to Shoup form to vOutS.
func SFormTo(vOutS, v []uint64, q *num.Modulus) {
	M := (len(vOutS) >> 3) << 3

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&vOutS[i]))
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))

		wOut[0] = modops.SForm(w0[0], qv)
		wOut[1] = modops.SForm(w0[1], qv)
		wOut[2] = modops.SForm(w0[2], qv)
		wOut[3] = modops.SForm(w0[3], qv)

		wOut[4] = modops.SForm(w0[4], qv)
		wOut[5] = modops.SForm(w0[5], qv)
		wOut[6] = modops.SForm(w0[6], qv)
		wOut[7] = modops.SForm(w0[7], qv)
	}

	for i := M; i < len(vOutS); i++ {
		vOutS[i] = modops.SForm(v[i], qv)
	}
}

// SMul returns v0 * v1 mod q using Shoup multiplication.
func SMul(v0, v1, v1S []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SMulTo(vOut, v0, v1, v1S, q)
	return vOut
}

// SMulTo computes vOut = v0 * v1 mod q using Shoup multiplication.
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
func Reduce(v []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v))
	ReduceTo(vOut, v, q)
	return vOut
}

// ReduceTo computes vOut = v mod q.
func ReduceTo(vOut, v []uint64, q *num.Modulus) {
	M := (len(vOut) >> 3) << 3

	qv := q.Value()
	divHi, _ := q.Div()

	for i := 0; i < M; i += 8 {
		w0 := (*[8]uint64)(unsafe.Pointer(&v[i]))
		wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i]))

		wOut[1] = modops.BMod64(w0[1], qv, divHi)
		wOut[0] = modops.BMod64(w0[0], qv, divHi)
		wOut[2] = modops.BMod64(w0[2], qv, divHi)
		wOut[3] = modops.BMod64(w0[3], qv, divHi)

		wOut[4] = modops.BMod64(w0[4], qv, divHi)
		wOut[5] = modops.BMod64(w0[5], qv, divHi)
		wOut[6] = modops.BMod64(w0[6], qv, divHi)
		wOut[7] = modops.BMod64(w0[7], qv, divHi)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.BMod64(v[i], qv, divHi)
	}
}
