//go:build amd64 && !purego

package dft

import (
	"unsafe"

	"golang.org/x/sys/cpu"
)

// fwdNTTInPlacePow2Unroll computes the NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= 16.
func fwdNTTInPlacePow2Unroll(coeffs, tw, twS []uint64, q uint64) {
	switch {
	case cpu.X86.HasAVX2 && cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasAVX512VL && cpu.X86.HasBMI2:
		fwdNTTInPlacePow2UnrollAVX512(coeffs, tw, twS, q)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2 && cpu.X86.HasBMI2:
		fwdNTTInPlacePow2UnrollAVX2(coeffs, tw, twS, q)
		return
	}
	N := len(coeffs)
	twoQ := q << 1
	twPtr := unsafe.Pointer(unsafe.SliceData(tw))
	twSPtr := unsafe.Pointer(unsafe.SliceData(twS))
	coeffsPtr := unsafe.Pointer(unsafe.SliceData(coeffs))
	var w, wS uint64

	t := N / 2
	w = *(*uint64)(unsafe.Add(twPtr, 1*unsafe.Sizeof(uint64(0))))
	wS = *(*uint64)(unsafe.Add(twSPtr, 1*unsafe.Sizeof(uint64(0))))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))
		c1 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j+t)*unsafe.Sizeof(uint64(0))))

		c0[0], c1[0] = butterflyPow2(c0[0], c1[0], w, wS, q, twoQ)
		c0[1], c1[1] = butterflyPow2(c0[1], c1[1], w, wS, q, twoQ)
		c0[2], c1[2] = butterflyPow2(c0[2], c1[2], w, wS, q, twoQ)
		c0[3], c1[3] = butterflyPow2(c0[3], c1[3], w, wS, q, twoQ)

		c0[4], c1[4] = butterflyPow2(c0[4], c1[4], w, wS, q, twoQ)
		c0[5], c1[5] = butterflyPow2(c0[5], c1[5], w, wS, q, twoQ)
		c0[6], c1[6] = butterflyPow2(c0[6], c1[6], w, wS, q, twoQ)
		c0[7], c1[7] = butterflyPow2(c0[7], c1[7], w, wS, q, twoQ)
	}

	for m := 2; m <= N/16; m <<= 1 {
		t >>= 1
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t

			w = *(*uint64)(unsafe.Add(twPtr, uintptr(m+i)*unsafe.Sizeof(uint64(0))))
			wS = *(*uint64)(unsafe.Add(twSPtr, uintptr(m+i)*unsafe.Sizeof(uint64(0))))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))
				c1 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j+t)*unsafe.Sizeof(uint64(0))))

				c0[0], c1[0] = butterflyPow2(c0[0], c1[0], w, wS, q, twoQ)
				c0[1], c1[1] = butterflyPow2(c0[1], c1[1], w, wS, q, twoQ)
				c0[2], c1[2] = butterflyPow2(c0[2], c1[2], w, wS, q, twoQ)
				c0[3], c1[3] = butterflyPow2(c0[3], c1[3], w, wS, q, twoQ)

				c0[4], c1[4] = butterflyPow2(c0[4], c1[4], w, wS, q, twoQ)
				c0[5], c1[5] = butterflyPow2(c0[5], c1[5], w, wS, q, twoQ)
				c0[6], c1[6] = butterflyPow2(c0[6], c1[6], w, wS, q, twoQ)
				c0[7], c1[7] = butterflyPow2(c0[7], c1[7], w, wS, q, twoQ)
			}
		}
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		w = *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/8)*unsafe.Sizeof(uint64(0))))
		wS = *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/8)*unsafe.Sizeof(uint64(0))))

		j := 8 * i

		c := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))

		c[0], c[4] = butterflyPow2(c[0], c[4], w, wS, q, twoQ)
		c[1], c[5] = butterflyPow2(c[1], c[5], w, wS, q, twoQ)
		c[2], c[6] = butterflyPow2(c[2], c[6], w, wS, q, twoQ)
		c[3], c[7] = butterflyPow2(c[3], c[7], w, wS, q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		j := 4 * i

		c := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))

		w = *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/4)*unsafe.Sizeof(uint64(0))))
		wS = *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/4)*unsafe.Sizeof(uint64(0))))
		c[0], c[2] = butterflyPow2(c[0], c[2], w, wS, q, twoQ)
		c[1], c[3] = butterflyPow2(c[1], c[3], w, wS, q, twoQ)

		w = *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/4+1)*unsafe.Sizeof(uint64(0))))
		wS = *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/4+1)*unsafe.Sizeof(uint64(0))))
		c[4], c[6] = butterflyPow2(c[4], c[6], w, wS, q, twoQ)
		c[5], c[7] = butterflyPow2(c[5], c[7], w, wS, q, twoQ)
	}

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		j := 2 * i

		c := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))

		c[0], c[1] = butterflyPow2(c[0], c[1], *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/2+0)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/2+0)*unsafe.Sizeof(uint64(0)))), q, twoQ)
		c[2], c[3] = butterflyPow2(c[2], c[3], *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/2+1)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/2+1)*unsafe.Sizeof(uint64(0)))), q, twoQ)
		c[4], c[5] = butterflyPow2(c[4], c[5], *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/2+2)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/2+2)*unsafe.Sizeof(uint64(0)))), q, twoQ)
		c[6], c[7] = butterflyPow2(c[6], c[7], *(*uint64)(unsafe.Add(twPtr, uintptr(i+N/2+3)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twSPtr, uintptr(i+N/2+3)*unsafe.Sizeof(uint64(0)))), q, twoQ)
	}
}

// invNTTInPlacePow2Unroll computes the Inverse NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= 32.
func invNTTInPlacePow2Unroll(coeffs, twInv, twInvS []uint64, q uint64) {
	switch {
	case cpu.X86.HasAVX2 && cpu.X86.HasAVX512DQ && cpu.X86.HasAVX512F && cpu.X86.HasAVX512VL && cpu.X86.HasBMI2:
		invNTTInPlacePow2UnrollAVX512(coeffs, twInv, twInvS, q)
		return
	case cpu.X86.HasAVX && cpu.X86.HasAVX2 && cpu.X86.HasBMI2:
		invNTTInPlacePow2UnrollAVX2(coeffs, twInv, twInvS, q)
		return
	}

	N := len(coeffs)
	twoQ := q << 1
	twInvPtr := unsafe.Pointer(unsafe.SliceData(twInv))
	twInvSPtr := unsafe.Pointer(unsafe.SliceData(twInvS))
	coeffsPtr := unsafe.Pointer(unsafe.SliceData(coeffs))
	var w, wS uint64

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		j := 2 * i

		c := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))

		c[0], c[1] = invButterflyPow2(c[0], c[1], *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/2+0)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/2+0)*unsafe.Sizeof(uint64(0)))), q, twoQ)
		c[2], c[3] = invButterflyPow2(c[2], c[3], *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/2+1)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/2+1)*unsafe.Sizeof(uint64(0)))), q, twoQ)
		c[4], c[5] = invButterflyPow2(c[4], c[5], *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/2+2)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/2+2)*unsafe.Sizeof(uint64(0)))), q, twoQ)
		c[6], c[7] = invButterflyPow2(c[6], c[7], *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/2+3)*unsafe.Sizeof(uint64(0)))), *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/2+3)*unsafe.Sizeof(uint64(0)))), q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		j := 4 * i

		c := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))

		w = *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/4)*unsafe.Sizeof(uint64(0))))
		wS = *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/4)*unsafe.Sizeof(uint64(0))))
		c[0], c[2] = invButterflyPow2(c[0], c[2], w, wS, q, twoQ)
		c[1], c[3] = invButterflyPow2(c[1], c[3], w, wS, q, twoQ)

		w = *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/4+1)*unsafe.Sizeof(uint64(0))))
		wS = *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/4+1)*unsafe.Sizeof(uint64(0))))
		c[4], c[6] = invButterflyPow2(c[4], c[6], w, wS, q, twoQ)
		c[5], c[7] = invButterflyPow2(c[5], c[7], w, wS, q, twoQ)
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		w = *(*uint64)(unsafe.Add(twInvPtr, uintptr(i+N/8)*unsafe.Sizeof(uint64(0))))
		wS = *(*uint64)(unsafe.Add(twInvSPtr, uintptr(i+N/8)*unsafe.Sizeof(uint64(0))))

		j := 8 * i

		c := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))

		c[0], c[4] = invButterflyPow2(c[0], c[4], w, wS, q, twoQ)
		c[1], c[5] = invButterflyPow2(c[1], c[5], w, wS, q, twoQ)
		c[2], c[6] = invButterflyPow2(c[2], c[6], w, wS, q, twoQ)
		c[3], c[7] = invButterflyPow2(c[3], c[7], w, wS, q, twoQ)
	}

	t := 8
	for m := N / 16; m >= 2; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t

			w = *(*uint64)(unsafe.Add(twInvPtr, uintptr(m+i)*unsafe.Sizeof(uint64(0))))
			wS = *(*uint64)(unsafe.Add(twInvSPtr, uintptr(m+i)*unsafe.Sizeof(uint64(0))))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))
				c1 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j+t)*unsafe.Sizeof(uint64(0))))

				c0[0], c1[0] = invButterflyPow2(c0[0], c1[0], w, wS, q, twoQ)
				c0[1], c1[1] = invButterflyPow2(c0[1], c1[1], w, wS, q, twoQ)
				c0[2], c1[2] = invButterflyPow2(c0[2], c1[2], w, wS, q, twoQ)
				c0[3], c1[3] = invButterflyPow2(c0[3], c1[3], w, wS, q, twoQ)

				c0[4], c1[4] = invButterflyPow2(c0[4], c1[4], w, wS, q, twoQ)
				c0[5], c1[5] = invButterflyPow2(c0[5], c1[5], w, wS, q, twoQ)
				c0[6], c1[6] = invButterflyPow2(c0[6], c1[6], w, wS, q, twoQ)
				c0[7], c1[7] = invButterflyPow2(c0[7], c1[7], w, wS, q, twoQ)
			}
		}
		t <<= 1
	}

	w = *(*uint64)(unsafe.Add(twInvPtr, 1*unsafe.Sizeof(uint64(0))))
	wS = *(*uint64)(unsafe.Add(twInvSPtr, 1*unsafe.Sizeof(uint64(0))))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j)*unsafe.Sizeof(uint64(0))))
		c1 := (*[8]uint64)(unsafe.Add(coeffsPtr, uintptr(j+t)*unsafe.Sizeof(uint64(0))))

		c0[0], c1[0] = invButterflyPow2(c0[0], c1[0], w, wS, q, twoQ)
		c0[1], c1[1] = invButterflyPow2(c0[1], c1[1], w, wS, q, twoQ)
		c0[2], c1[2] = invButterflyPow2(c0[2], c1[2], w, wS, q, twoQ)
		c0[3], c1[3] = invButterflyPow2(c0[3], c1[3], w, wS, q, twoQ)

		c0[4], c1[4] = invButterflyPow2(c0[4], c1[4], w, wS, q, twoQ)
		c0[5], c1[5] = invButterflyPow2(c0[5], c1[5], w, wS, q, twoQ)
		c0[6], c1[6] = invButterflyPow2(c0[6], c1[6], w, wS, q, twoQ)
		c0[7], c1[7] = invButterflyPow2(c0[7], c1[7], w, wS, q, twoQ)
	}
}
