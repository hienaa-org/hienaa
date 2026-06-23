//go:build amd64 && !purego

package dft

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/modops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
	"golang.org/x/sys/cpu"
)

// fwdNTTInPlacePow2Unroll computes the NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= [nttUnrollBound].
func fwdNTTInPlacePow2Unroll(coeffs, tw []uint64, twM vec.MulForm, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512IFMA && cpu.X86.HasAVX512F && cpu.X86.HasAVX512VL:
		fwdNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, tw, twM.SForm, q.Value(), 0, uint64(len(coeffs)), 1)
		return
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasAVX512VL && cpu.X86.HasBMI2:
		fwdNTTInPlacePow2UnrollAVX512(coeffs, twM.Float, q.Float(), q.FloatInv())
		return
	case cpu.X86.HasAVX2 && cpu.X86.HasAVX && cpu.X86.HasFMA:
		fwdNTTInPlacePow2UnrollAVX2(coeffs, twM.Float, q.Float(), q.FloatInv())
		return
	}

	qv := q.Value()
	twoQv := qv << 1

	N := len(coeffs)
	L := unsafe.Sizeof(uint64(0))

	r := unsafe.Pointer(unsafe.SliceData(tw))
	rS := unsafe.Pointer(unsafe.SliceData(twM.SForm))
	v := unsafe.Pointer(unsafe.SliceData(coeffs))

	t := N / 2
	w := *(*uint64)(unsafe.Add(r, 1*L))
	wS := *(*uint64)(unsafe.Add(rS, 1*L))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
		c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

		c0[0], c1[0] = fwdButterflyPow2(c0[0], c1[0], w, wS, qv, twoQv)
		c0[1], c1[1] = fwdButterflyPow2(c0[1], c1[1], w, wS, qv, twoQv)
		c0[2], c1[2] = fwdButterflyPow2(c0[2], c1[2], w, wS, qv, twoQv)
		c0[3], c1[3] = fwdButterflyPow2(c0[3], c1[3], w, wS, qv, twoQv)

		c0[4], c1[4] = fwdButterflyPow2(c0[4], c1[4], w, wS, qv, twoQv)
		c0[5], c1[5] = fwdButterflyPow2(c0[5], c1[5], w, wS, qv, twoQv)
		c0[6], c1[6] = fwdButterflyPow2(c0[6], c1[6], w, wS, qv, twoQv)
		c0[7], c1[7] = fwdButterflyPow2(c0[7], c1[7], w, wS, qv, twoQv)
	}

	for m := 2; m <= N/16; m <<= 1 {
		t >>= 1
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t

			w := *(*uint64)(unsafe.Add(r, uintptr(m+i)*L))
			wS := *(*uint64)(unsafe.Add(rS, uintptr(m+i)*L))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
				c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

				c0[0], c1[0] = fwdButterflyPow2(c0[0], c1[0], w, wS, qv, twoQv)
				c0[1], c1[1] = fwdButterflyPow2(c0[1], c1[1], w, wS, qv, twoQv)
				c0[2], c1[2] = fwdButterflyPow2(c0[2], c1[2], w, wS, qv, twoQv)
				c0[3], c1[3] = fwdButterflyPow2(c0[3], c1[3], w, wS, qv, twoQv)

				c0[4], c1[4] = fwdButterflyPow2(c0[4], c1[4], w, wS, qv, twoQv)
				c0[5], c1[5] = fwdButterflyPow2(c0[5], c1[5], w, wS, qv, twoQv)
				c0[6], c1[6] = fwdButterflyPow2(c0[6], c1[6], w, wS, qv, twoQv)
				c0[7], c1[7] = fwdButterflyPow2(c0[7], c1[7], w, wS, qv, twoQv)
			}
		}
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(8*i)*L))
		w := *(*uint64)(unsafe.Add(r, uintptr(i+N/8)*L))
		wS := *(*uint64)(unsafe.Add(rS, uintptr(i+N/8)*L))

		c[0], c[4] = fwdButterflyPow2(c[0], c[4], w, wS, qv, twoQv)
		c[1], c[5] = fwdButterflyPow2(c[1], c[5], w, wS, qv, twoQv)
		c[2], c[6] = fwdButterflyPow2(c[2], c[6], w, wS, qv, twoQv)
		c[3], c[7] = fwdButterflyPow2(c[3], c[7], w, wS, qv, twoQv)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(4*i)*L))
		w := (*[2]uint64)(unsafe.Add(r, uintptr(i+N/4)*L))
		wS := (*[2]uint64)(unsafe.Add(rS, uintptr(i+N/4)*L))

		c[0], c[2] = fwdButterflyPow2(c[0], c[2], w[0], wS[0], qv, twoQv)
		c[1], c[3] = fwdButterflyPow2(c[1], c[3], w[0], wS[0], qv, twoQv)

		c[4], c[6] = fwdButterflyPow2(c[4], c[6], w[1], wS[1], qv, twoQv)
		c[5], c[7] = fwdButterflyPow2(c[5], c[7], w[1], wS[1], qv, twoQv)
	}

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(2*i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(i+N/2)*L))
		wS := (*[8]uint64)(unsafe.Add(rS, uintptr(i+N/2)*L))

		c[0], c[1] = fwdButterflyPow2(c[0], c[1], w[0], wS[0], qv, twoQv)
		c[2], c[3] = fwdButterflyPow2(c[2], c[3], w[1], wS[1], qv, twoQv)
		c[4], c[5] = fwdButterflyPow2(c[4], c[5], w[2], wS[2], qv, twoQv)
		c[6], c[7] = fwdButterflyPow2(c[6], c[7], w[3], wS[3], qv, twoQv)
	}

	for i := 0; i < N; i += 8 {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(i)*L))

		c[0] = modops.Reduce4Q(c[0], qv, twoQv)
		c[1] = modops.Reduce4Q(c[1], qv, twoQv)
		c[2] = modops.Reduce4Q(c[2], qv, twoQv)
		c[3] = modops.Reduce4Q(c[3], qv, twoQv)

		c[4] = modops.Reduce4Q(c[4], qv, twoQv)
		c[5] = modops.Reduce4Q(c[5], qv, twoQv)
		c[6] = modops.Reduce4Q(c[6], qv, twoQv)
		c[7] = modops.Reduce4Q(c[7], qv, twoQv)
	}
}

func fwdNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, tw, twS []uint64, q uint64, idx, N, l uint64) {
	if N < nttRecurseBound {
		fwdNTTInPlacePow2UnrollAVX512IFMA(coeffs, tw, twS, q, idx, N, l)
		return
	}

	L := unsafe.Sizeof(uint64(0))

	r := unsafe.Pointer(unsafe.SliceData(tw))
	rS := unsafe.Pointer(unsafe.SliceData(twS))

	w := *(*uint64)(unsafe.Add(r, uintptr(l)*L))
	wS := *(*uint64)(unsafe.Add(rS, uintptr(l)*L))

	t := N >> 1

	fwdNTTInPlacePow2StrideUnrollAVX512IFMA(coeffs, w, wS, q, idx, t)

	fwdNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, tw, twS, q, idx, t, l<<1)
	fwdNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, tw, twS, q, idx+t, t, (l<<1)|1)
}

// invNTTInPlacePow2Unroll computes the Inverse NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= [nttUnrollBound].
func invNTTInPlacePow2Unroll(coeffs, twInv []uint64, twInvM vec.MulForm, q *num.Modulus) {
	switch {
	case cpu.X86.HasAVX512IFMA && cpu.X86.HasAVX512F && cpu.X86.HasAVX512VL:
		invNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, twInv, twInvM.SForm, q.Value(), 0, uint64(len(coeffs)), 1)
		return
	case cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasAVX512VL:
		invNTTInPlacePow2UnrollAVX512(coeffs, twInvM.Float, q.Float(), q.FloatInv())
		return
	case cpu.X86.HasAVX2 && cpu.X86.HasAVX && cpu.X86.HasFMA:
		invNTTInPlacePow2UnrollAVX2(coeffs, twInvM.Float, q.Float(), q.FloatInv())
		return
	}

	qv := q.Value()
	twoQv := qv << 1

	N := len(coeffs)
	L := unsafe.Sizeof(uint64(0))

	r := unsafe.Pointer(unsafe.SliceData(twInv))
	rS := unsafe.Pointer(unsafe.SliceData(twInvM.SForm))
	v := unsafe.Pointer(unsafe.SliceData(coeffs))

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(2*i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(i+N/2)*L))
		wS := (*[8]uint64)(unsafe.Add(rS, uintptr(i+N/2)*L))

		c[0], c[1] = invButterflyPow2(c[0], c[1], w[0], wS[0], qv, twoQv)
		c[2], c[3] = invButterflyPow2(c[2], c[3], w[1], wS[1], qv, twoQv)
		c[4], c[5] = invButterflyPow2(c[4], c[5], w[2], wS[2], qv, twoQv)
		c[6], c[7] = invButterflyPow2(c[6], c[7], w[3], wS[3], qv, twoQv)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(4*i)*L))
		w := (*[2]uint64)(unsafe.Add(r, uintptr(i+N/4)*L))
		wS := (*[2]uint64)(unsafe.Add(rS, uintptr(i+N/4)*L))

		c[0], c[2] = invButterflyPow2(c[0], c[2], w[0], wS[0], qv, twoQv)
		c[1], c[3] = invButterflyPow2(c[1], c[3], w[0], wS[0], qv, twoQv)

		c[4], c[6] = invButterflyPow2(c[4], c[6], w[1], wS[1], qv, twoQv)
		c[5], c[7] = invButterflyPow2(c[5], c[7], w[1], wS[1], qv, twoQv)
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(8*i)*L))
		w := *(*uint64)(unsafe.Add(r, uintptr(i+N/8)*L))
		wS := *(*uint64)(unsafe.Add(rS, uintptr(i+N/8)*L))

		c[0], c[4] = invButterflyPow2(c[0], c[4], w, wS, qv, twoQv)
		c[1], c[5] = invButterflyPow2(c[1], c[5], w, wS, qv, twoQv)
		c[2], c[6] = invButterflyPow2(c[2], c[6], w, wS, qv, twoQv)
		c[3], c[7] = invButterflyPow2(c[3], c[7], w, wS, qv, twoQv)
	}

	t := 8
	for m := N / 16; m >= 2; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t

			w := *(*uint64)(unsafe.Add(r, uintptr(m+i)*L))
			wS := *(*uint64)(unsafe.Add(rS, uintptr(m+i)*L))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
				c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

				c0[0], c1[0] = invButterflyPow2(c0[0], c1[0], w, wS, qv, twoQv)
				c0[1], c1[1] = invButterflyPow2(c0[1], c1[1], w, wS, qv, twoQv)
				c0[2], c1[2] = invButterflyPow2(c0[2], c1[2], w, wS, qv, twoQv)
				c0[3], c1[3] = invButterflyPow2(c0[3], c1[3], w, wS, qv, twoQv)

				c0[4], c1[4] = invButterflyPow2(c0[4], c1[4], w, wS, qv, twoQv)
				c0[5], c1[5] = invButterflyPow2(c0[5], c1[5], w, wS, qv, twoQv)
				c0[6], c1[6] = invButterflyPow2(c0[6], c1[6], w, wS, qv, twoQv)
				c0[7], c1[7] = invButterflyPow2(c0[7], c1[7], w, wS, qv, twoQv)
			}
		}
		t <<= 1
	}

	w := *(*uint64)(unsafe.Add(r, 1*L))
	wS := *(*uint64)(unsafe.Add(rS, 1*L))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
		c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

		c0[0], c1[0] = invButterflyPow2(c0[0], c1[0], w, wS, qv, twoQv)
		c0[1], c1[1] = invButterflyPow2(c0[1], c1[1], w, wS, qv, twoQv)
		c0[2], c1[2] = invButterflyPow2(c0[2], c1[2], w, wS, qv, twoQv)
		c0[3], c1[3] = invButterflyPow2(c0[3], c1[3], w, wS, qv, twoQv)

		c0[4], c1[4] = invButterflyPow2(c0[4], c1[4], w, wS, qv, twoQv)
		c0[5], c1[5] = invButterflyPow2(c0[5], c1[5], w, wS, qv, twoQv)
		c0[6], c1[6] = invButterflyPow2(c0[6], c1[6], w, wS, qv, twoQv)
		c0[7], c1[7] = invButterflyPow2(c0[7], c1[7], w, wS, qv, twoQv)
	}

	for i := 0; i < N; i += 8 {
		c := (*[8]uint64)(unsafe.Add(v, uintptr(i)*L))

		c[0] = modops.Reduce4Q(c[0], qv, twoQv)
		c[1] = modops.Reduce4Q(c[1], qv, twoQv)
		c[2] = modops.Reduce4Q(c[2], qv, twoQv)
		c[3] = modops.Reduce4Q(c[3], qv, twoQv)

		c[4] = modops.Reduce4Q(c[4], qv, twoQv)
		c[5] = modops.Reduce4Q(c[5], qv, twoQv)
		c[6] = modops.Reduce4Q(c[6], qv, twoQv)
		c[7] = modops.Reduce4Q(c[7], qv, twoQv)
	}
}

func invNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, twInv, twInvS []uint64, q uint64, idx, N, l uint64) {
	if N < nttRecurseBound {
		invNTTInPlacePow2UnrollAVX512IFMA(coeffs, twInv, twInvS, q, idx, N, l)
		return
	}

	L := unsafe.Sizeof(uint64(0))

	r := unsafe.Pointer(unsafe.SliceData(twInv))
	rS := unsafe.Pointer(unsafe.SliceData(twInvS))

	w := *(*uint64)(unsafe.Add(r, uintptr(l)*L))
	wS := *(*uint64)(unsafe.Add(rS, uintptr(l)*L))

	t := N >> 1

	invNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, twInv, twInvS, q, idx, t, l<<1)
	invNTTInPlacePow2UnrollRecurseAVX512IFMA(coeffs, twInv, twInvS, q, idx+t, t, (l<<1)|1)

	invNTTInPlacePow2StrideUnrollAVX512IFMA(coeffs, w, wS, q, idx, t)
}
