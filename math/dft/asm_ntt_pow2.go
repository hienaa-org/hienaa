//go:build !(amd64 && !purego)

package dft

import (
	"unsafe"
)

// fwdNTTInPlacePow2Unroll computes the NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= 32.
func fwdNTTInPlacePow2Unroll(coeffs, tw, twS []uint64, q uint64) {
	twoQ := q << 1

	N := len(coeffs)
	L := unsafe.Sizeof(uint64(0))

	r := unsafe.Pointer(unsafe.SliceData(tw))
	rS := unsafe.Pointer(unsafe.SliceData(twS))
	v := unsafe.Pointer(unsafe.SliceData(coeffs))
	var w, wS uint64

	t := N / 2
	w = *(*uint64)(unsafe.Add(r, 1*L))
	wS = *(*uint64)(unsafe.Add(rS, 1*L))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
		c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

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

			w = *(*uint64)(unsafe.Add(r, uintptr(m+i)*L))
			wS = *(*uint64)(unsafe.Add(rS, uintptr(m+i)*L))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
				c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

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
		w = *(*uint64)(unsafe.Add(r, uintptr(i+N/8)*L))
		wS = *(*uint64)(unsafe.Add(rS, uintptr(i+N/8)*L))

		j := 8 * i

		c := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))

		c[0], c[4] = butterflyPow2(c[0], c[4], w, wS, q, twoQ)
		c[1], c[5] = butterflyPow2(c[1], c[5], w, wS, q, twoQ)
		c[2], c[6] = butterflyPow2(c[2], c[6], w, wS, q, twoQ)
		c[3], c[7] = butterflyPow2(c[3], c[7], w, wS, q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		j := 4 * i

		c := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))

		w = *(*uint64)(unsafe.Add(r, uintptr(i+N/4)*L))
		wS = *(*uint64)(unsafe.Add(rS, uintptr(i+N/4)*L))
		c[0], c[2] = butterflyPow2(c[0], c[2], w, wS, q, twoQ)
		c[1], c[3] = butterflyPow2(c[1], c[3], w, wS, q, twoQ)

		w = *(*uint64)(unsafe.Add(r, uintptr(i+N/4+1)*L))
		wS = *(*uint64)(unsafe.Add(rS, uintptr(i+N/4+1)*L))
		c[4], c[6] = butterflyPow2(c[4], c[6], w, wS, q, twoQ)
		c[5], c[7] = butterflyPow2(c[5], c[7], w, wS, q, twoQ)
	}

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		j := 2 * i

		c := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))

		c[0], c[1] = butterflyPow2(c[0], c[1], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+0)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+0)*L)), q, twoQ)
		c[2], c[3] = butterflyPow2(c[2], c[3], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+1)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+1)*L)), q, twoQ)
		c[4], c[5] = butterflyPow2(c[4], c[5], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+2)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+2)*L)), q, twoQ)
		c[6], c[7] = butterflyPow2(c[6], c[7], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+3)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+3)*L)), q, twoQ)
	}
}

// invNTTInPlacePow2Unroll computes the Inverse NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= 32.
func invNTTInPlacePow2Unroll(coeffs, twInv, twInvS []uint64, q uint64) {
	twoQ := q << 1

	N := len(coeffs)
	L := unsafe.Sizeof(uint64(0))

	r := unsafe.Pointer(unsafe.SliceData(twInv))
	rS := unsafe.Pointer(unsafe.SliceData(twInvS))
	v := unsafe.Pointer(unsafe.SliceData(coeffs))
	var w, wS uint64

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		j := 2 * i

		c := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))

		c[0], c[1] = invButterflyPow2(c[0], c[1], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+0)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+0)*L)), q, twoQ)
		c[2], c[3] = invButterflyPow2(c[2], c[3], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+1)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+1)*L)), q, twoQ)
		c[4], c[5] = invButterflyPow2(c[4], c[5], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+2)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+2)*L)), q, twoQ)
		c[6], c[7] = invButterflyPow2(c[6], c[7], *(*uint64)(unsafe.Add(r, uintptr(i+N/2+3)*L)), *(*uint64)(unsafe.Add(rS, uintptr(i+N/2+3)*L)), q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		j := 4 * i

		c := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))

		w = *(*uint64)(unsafe.Add(r, uintptr(i+N/4)*L))
		wS = *(*uint64)(unsafe.Add(rS, uintptr(i+N/4)*L))
		c[0], c[2] = invButterflyPow2(c[0], c[2], w, wS, q, twoQ)
		c[1], c[3] = invButterflyPow2(c[1], c[3], w, wS, q, twoQ)

		w = *(*uint64)(unsafe.Add(r, uintptr(i+N/4+1)*L))
		wS = *(*uint64)(unsafe.Add(rS, uintptr(i+N/4+1)*L))
		c[4], c[6] = invButterflyPow2(c[4], c[6], w, wS, q, twoQ)
		c[5], c[7] = invButterflyPow2(c[5], c[7], w, wS, q, twoQ)
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		w = *(*uint64)(unsafe.Add(r, uintptr(i+N/8)*L))
		wS = *(*uint64)(unsafe.Add(rS, uintptr(i+N/8)*L))

		j := 8 * i

		c := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))

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

			w = *(*uint64)(unsafe.Add(r, uintptr(m+i)*L))
			wS = *(*uint64)(unsafe.Add(rS, uintptr(m+i)*L))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
				c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

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

	w = *(*uint64)(unsafe.Add(r, 1*L))
	wS = *(*uint64)(unsafe.Add(rS, 1*L))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Add(v, uintptr(j)*L))
		c1 := (*[8]uint64)(unsafe.Add(v, uintptr(j+t)*L))

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
