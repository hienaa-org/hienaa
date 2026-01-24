//go:build amd64 && !purego

package vec

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
	"golang.org/x/sys/cpu"
)

// AddTo computes vOut = v0 + v1 mod q.
// x0 and x1 must be in [0, q).
// If q is nil, then it returns x0 + x1.
func AddTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		addTo(vOut, v0, v1, q)
		return
	}
	addWordTo(vOut, v0, v1)
}

// addTo computes vOut = v0 + v1 mod q.
func addTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		addToAVX512(vOut, v0, v1, q.Value())
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		addToAVX2(vOut, v0, v1, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// addWordTo computes vOut = v0 + v1.
func addWordTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		addWordToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		addWordToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// ScalarAddTo computes vOut = v + c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v + c.
func ScalarAddTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		scalarAddTo(vOut, v, c, q)
		return
	}
	scalarAddWordTo(vOut, v, c)
}

// scalarAddTo computes vOut = v + c mod q.
func scalarAddTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		scalarAddToAVX512(vOut, v, c, q.Value())
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarAddToAVX2(vOut, v, c, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] = modops.Add(w[0], c, qv)
		wOut[1] = modops.Add(w[1], c, qv)
		wOut[2] = modops.Add(w[2], c, qv)
		wOut[3] = modops.Add(w[3], c, qv)

		wOut[4] = modops.Add(w[4], c, qv)
		wOut[5] = modops.Add(w[5], c, qv)
		wOut[6] = modops.Add(w[6], c, qv)
		wOut[7] = modops.Add(w[7], c, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Add(v[i], c, qv)
	}
}

// scalarAddWordTo computes vOut = v + c.
func scalarAddWordTo(vOut, v []uint64, c uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		scalarAddWordToAVX512(vOut, v, c)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarAddWordToAVX2(vOut, v, c)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] = w[0] + c
		wOut[1] = w[1] + c
		wOut[2] = w[2] + c
		wOut[3] = w[3] + c

		wOut[4] = w[4] + c
		wOut[5] = w[5] + c
		wOut[6] = w[6] + c
		wOut[7] = w[7] + c
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = v[i] + c
	}
}

// Sub returns v0 - v1 mod q.
// v0 and v1 must be in [0, q).
// If q is nil, then it returns v0 - v1.
func SubTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	if q != nil {
		subTo(vOut, v0, v1, q)
		return
	}
	subWordTo(vOut, v0, v1)
}

// subTo computes vOut = v0 - v1 mod q.
func subTo(vOut, v0, v1 []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		subToAVX512(vOut, v0, v1, q.Value())
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		subToAVX2(vOut, v0, v1, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// subWordTo computes vOut = v0 - v1.
func subWordTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		subWordToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		subWordToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// ScalarSub returns v - c mod q.
// v and c must be in [0, q).
// If q is nil, then it returns v - c.
func ScalarSubTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	if q != nil {
		scalarSubTo(vOut, v, c, q)
		return
	}
	scalarSubWordTo(vOut, v, c)
}

// scalarSubTo computes vOut = v - c mod q.
func scalarSubTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		scalarSubToAVX512(vOut, v, c, q.Value())
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarSubToAVX2(vOut, v, c, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] = modops.Sub(w[0], c, qv)
		wOut[1] = modops.Sub(w[1], c, qv)
		wOut[2] = modops.Sub(w[2], c, qv)
		wOut[3] = modops.Sub(w[3], c, qv)

		wOut[4] = modops.Sub(w[4], c, qv)
		wOut[5] = modops.Sub(w[5], c, qv)
		wOut[6] = modops.Sub(w[6], c, qv)
		wOut[7] = modops.Sub(w[7], c, qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Sub(v[i], c, qv)
	}
}

// scalarSubWordTo computes vOut = v - c.
func scalarSubWordTo(vOut, v []uint64, c uint64) {
	switch {
	case cpu.X86.HasAVX512F:
		scalarSubWordToAVX512(vOut, v, c)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarSubWordToAVX2(vOut, v, c)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] = w[0] - c
		wOut[1] = w[1] - c
		wOut[2] = w[2] - c
		wOut[3] = w[3] - c

		wOut[4] = w[4] - c
		wOut[5] = w[5] - c
		wOut[6] = w[6] - c
		wOut[7] = w[7] - c
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = v[i] - c
	}
}

// negTo computes vOut = -v mod q.
func negTo(vOut, v []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512F:
		negToAVX512(vOut, v, q.Value())
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		negToAVX2(vOut, v, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512F:
		negWordToAVX512(vOut, v)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		negWordToAVX2(vOut, v)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// MFormTo transforms v to Montgomery form to vOutM.
//
// Panics if q is even or nil.
func MFormTo(vOutM, v []uint64, q *num.Modulus) {
	if q.Inv() == 0 {
		panic("modulus must be odd")
	}

	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		divHi, divLo := q.Div()
		mFormToAVX512(vOutM, v, q.Value(), divHi, divLo)
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, divLo := q.Div()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// InvMFormTo computes vOut as vM in Normal form.
//
// Panics if q is even or nil.
func InvMFormTo(vOut, vM []uint64, q *num.Modulus) {
	if q.Inv() == 0 {
		panic("modulus must be odd")
	}

	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		invMFormToAVX512(vOut, vM, q.Value(), q.Inv())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// scalarMulTo computes vOut = c * v mod q using Shoup multiplication.
func scalarMulTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		divHi, _ := q.Div()
		scalarMulToAVX512(vOut, v, c, modops.SForm(c, qv, divHi), qv)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// scalarMulAddTo computes vOut += c * v mod q using Shoup multiplication.
func scalarMulAddTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		divHi, _ := q.Div()
		scalarMulAddToAVX512(vOut, v, c, modops.SForm(c, qv, divHi), qv)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] = modops.Add(wOut[0], modops.SMul(w[0], c, cS, qv), qv)
		wOut[1] = modops.Add(wOut[1], modops.SMul(w[1], c, cS, qv), qv)
		wOut[2] = modops.Add(wOut[2], modops.SMul(w[2], c, cS, qv), qv)
		wOut[3] = modops.Add(wOut[3], modops.SMul(w[3], c, cS, qv), qv)

		wOut[4] = modops.Add(wOut[4], modops.SMul(w[4], c, cS, qv), qv)
		wOut[5] = modops.Add(wOut[5], modops.SMul(w[5], c, cS, qv), qv)
		wOut[6] = modops.Add(wOut[6], modops.SMul(w[6], c, cS, qv), qv)
		wOut[7] = modops.Add(wOut[7], modops.SMul(w[7], c, cS, qv), qv)
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = modops.Add(vOut[i], modops.SMul(v[i], c, cS, qv), qv)
	}
}

// scalarMulSubTo computes vOut -= c * v mod q using Shoup multiplication.
func scalarMulSubTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		divHi, _ := q.Div()
		scalarMulSubToAVX512(vOut, v, c, modops.SForm(c, qv, divHi), qv)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// scalarMulWordTo computes vOut = c * v.
func scalarMulWordTo(vOut, v []uint64, c uint64) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMulWordToAVX512(vOut, v, c)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarMulWordToAVX2(vOut, v, c)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// scalarMulAddWordTo computes vOut += c * v.
func scalarMulAddWordTo(vOut, v []uint64, c uint64) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMulAddWordToAVX512(vOut, v, c)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarMulAddWordToAVX2(vOut, v, c)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// scalarMulSubWordTo computes vOut -= c * v.
func scalarMulSubWordTo(vOut, v []uint64, c uint64) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMulSubWordToAVX512(vOut, v, c)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		scalarMulSubWordToAVX2(vOut, v, c)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// ScalarMulLazyTo computes vOut = c * v mod q using Shoup multiplication,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func ScalarMulLazyTo(vOut, v []uint64, c uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		divHi, _ := q.Div()
		scalarMulLazyToAVX512(vOut, v, c, modops.SForm(c, qv, divHi), qv)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		divHi, _ := q.Div()
		scalarMulAddLazyToAVX512(vOut, v, c, modops.SForm(c, qv, divHi), qv)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, _ := q.Div()
	cS := modops.SForm(c, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		divHi, _ := q.Div()
		cNeg := modops.Neg(c, qv)
		cNegS := modops.SForm(cNeg, qv, divHi)
		scalarMulSubLazyToAVX512(vOut, v, cNeg, cNegS, qv)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr := unsafe.Pointer(unsafe.SliceData(v))

	qv := q.Value()
	divHi, _ := q.Div()
	cNeg := modops.Neg(c, qv)
	cNegS := modops.SForm(cNeg, qv, divHi)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// ScalarMMulTo computes vOut = c * v mod q using Montgomery multiplication.
// When c is in Montgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMMulToAVX512(vOutM, vM, cM, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
// When c is in Montgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulAddTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMMulAddToAVX512(vOutM, vM, cM, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
// When c is in Montgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulSubTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMMulSubToAVX512(vOutM, vM, cM, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
// When c is in Montgomery form, the output is the same form as v.
//
// Panics if q is nil.
func ScalarMMulLazy(vM []uint64, cM uint64, q *num.Modulus) []uint64 {
	vOutM := make([]uint64, len(vM))
	ScalarMMulLazyTo(vOutM, vM, cM, q)
	return vOutM
}

// ScalarMMulLazyTo computes vOut = c * v mod q using Montgomery multiplication,
// but the result is in [0, 2q).
// When c is in Montgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMMulLazyToAVX512(vOutM, vM, cM, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
// When c is in Montgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulAddLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		scalarMMulAddLazyToAVX512(vOutM, vM, cM, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
// When c is in Montgomery form, vOut is the same form as v.
//
// Panics if q is nil.
func ScalarMMulSubLazyTo(vOutM, vM []uint64, cM uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		qv := q.Value()
		scalarMMulSubLazyToAVX512(vOutM, vM, modops.Neg(cM, qv), qv, q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr := unsafe.Pointer(unsafe.SliceData(vM))

	qv := q.Value()
	inv := q.Inv()

	cMNeg := modops.Neg(cM, qv)

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w := (*[8]uint64)(unsafe.Add(ptr, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// mulWordTo computes vOut = v0 * v1.
func mulWordTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mulWordToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		mulWordToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] = w0[0] * w1[0]
		wOut[1] = w0[1] * w1[1]
		wOut[2] = w0[2] * w1[2]
		wOut[3] = w0[3] * w1[3]

		wOut[4] = w0[4] * w1[4]
		wOut[5] = w0[5] * w1[5]
		wOut[6] = w0[6] * w1[6]
		wOut[7] = w0[7] * w1[7]
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] = v0[i] * v1[i]
	}
}

// mulAddWordTo computes vOut += v0 * v1.
func mulAddWordTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mulAddWordToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		mulAddWordToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] += w0[0] * w1[0]
		wOut[1] += w0[1] * w1[1]
		wOut[2] += w0[2] * w1[2]
		wOut[3] += w0[3] * w1[3]

		wOut[4] += w0[4] * w1[4]
		wOut[5] += w0[5] * w1[5]
		wOut[6] += w0[6] * w1[6]
		wOut[7] += w0[7] * w1[7]
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] += v0[i] * v1[i]
	}
}

// mulSubWordTo computes vOut -= v0 * v1.
func mulSubWordTo(vOut, v0, v1 []uint64) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mulSubWordToAVX512(vOut, v0, v1)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2:
		mulSubWordToAVX2(vOut, v0, v1)
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

		wOut[0] -= w0[0] * w1[0]
		wOut[1] -= w0[1] * w1[1]
		wOut[2] -= w0[2] * w1[2]
		wOut[3] -= w0[3] * w1[3]

		wOut[4] -= w0[4] * w1[4]
		wOut[5] -= w0[5] * w1[5]
		wOut[6] -= w0[6] * w1[6]
		wOut[7] -= w0[7] * w1[7]
	}

	for i := M; i < len(vOut); i++ {
		vOut[i] -= v0[i] * v1[i]
	}
}

// MMulTo computes vOut = v0 * v1 mod q in Montgomery form.
//
// Panics if q is nil.
func MMulTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mMulToAVX512(vOutM, v0M, v1M, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0M))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1M))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
//
// Panics if q is nil.
func MMulAddTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mMulAddToAVX512(vOutM, v0M, v1M, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0M))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1M))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
//
// Panics if q is nil.
func MMulSubTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mMulSubToAVX512(vOutM, v0M, v1M, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0M))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1M))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// MMulLazyTo computes vOut = v0 * v1 mod q in Montgomery form,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func MMulLazyTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mMulLazyToAVX512(vOutM, v0M, v1M, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0M))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1M))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
//
// Panics if q is nil.
func MMulAddLazyTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mMulAddLazyToAVX512(vOutM, v0M, v1M, q.Value(), q.Inv())
		return
	}
	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0M))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1M))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
//
// Panics if q is nil.
func MMulSubLazyTo(vOutM, v0M, v1M []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		mMulSubLazyToAVX512(vOutM, v0M, v1M, q.Value(), q.Inv())
		return
	}

	M := (len(vOutM) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOutM))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0M))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1M))

	qv := q.Value()
	inv := q.Inv()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))

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

// SMulTo computes vOut = v0 * v1 mod q using Shoup multiplication.
//
// Panics if q is nil.
func SMulTo(vOut, v0, v1, v1S []uint64, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		sMulToAVX512(vOut, v0, v1, v1S, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))
	ptr1S := unsafe.Pointer(unsafe.SliceData(v1S))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1S := (*[8]uint64)(unsafe.Add(ptr1S, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		sMulAddToAVX512(vOut, v0, v1, v1S, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))
	ptr1S := unsafe.Pointer(unsafe.SliceData(v1S))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1S := (*[8]uint64)(unsafe.Add(ptr1S, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		sMulSubToAVX512(vOut, v0, v1, v1S, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))
	ptr1S := unsafe.Pointer(unsafe.SliceData(v1S))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1S := (*[8]uint64)(unsafe.Add(ptr1S, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		sMulLazyToAVX512(vOut, v0, v1, v1S, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))
	ptr1S := unsafe.Pointer(unsafe.SliceData(v1S))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1S := (*[8]uint64)(unsafe.Add(ptr1S, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		sMulAddLazyToAVX512(vOut, v0, v1, v1S, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))
	ptr1S := unsafe.Pointer(unsafe.SliceData(v1S))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1S := (*[8]uint64)(unsafe.Add(ptr1S, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
	switch {
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasBMI2:
		sMulSubLazyToAVX512(vOut, v0, v1, v1S, q.Value())
		return
	}

	M := (len(vOut) >> 3) << 3
	ptrOut := unsafe.Pointer(unsafe.SliceData(vOut))
	ptr0 := unsafe.Pointer(unsafe.SliceData(v0))
	ptr1 := unsafe.Pointer(unsafe.SliceData(v1))
	ptr1S := unsafe.Pointer(unsafe.SliceData(v1S))

	qv := q.Value()

	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(ptrOut, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w0 := (*[8]uint64)(unsafe.Add(ptr0, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1 := (*[8]uint64)(unsafe.Add(ptr1, uintptr(i)*unsafe.Sizeof(uint64(0))))
		w1S := (*[8]uint64)(unsafe.Add(ptr1S, uintptr(i)*unsafe.Sizeof(uint64(0))))

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
