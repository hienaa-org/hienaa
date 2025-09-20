//go:build !(amd64 && !purego)

package dft

import (
	"unsafe"
)

// nttInPlacePow2Unroll computes the NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= 32.
func nttInPlacePow2Unroll(coeffs, tw, twS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1
	var w, wS uint64

	t := N / 2
	w, wS = tw[1], twS[1]
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))
		c1 := (*[8]uint64)(unsafe.Pointer(&coeffs[j+t]))

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

			w, wS = tw[m+i], twS[m+i]

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))
				c1 := (*[8]uint64)(unsafe.Pointer(&coeffs[j+t]))

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
		w, wS = tw[i+N/8], twS[i+N/8]

		j := 8 * i

		c := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))

		c[0], c[4] = butterflyPow2(c[0], c[4], w, wS, q, twoQ)
		c[1], c[5] = butterflyPow2(c[1], c[5], w, wS, q, twoQ)
		c[2], c[6] = butterflyPow2(c[2], c[6], w, wS, q, twoQ)
		c[3], c[7] = butterflyPow2(c[3], c[7], w, wS, q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		j := 4 * i

		c := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))

		w, wS = tw[i+N/4], twS[i+N/4]
		c[0], c[2] = butterflyPow2(c[0], c[2], w, wS, q, twoQ)
		c[1], c[3] = butterflyPow2(c[1], c[3], w, wS, q, twoQ)

		w, wS = tw[i+N/4+1], twS[i+N/4+1]
		c[4], c[6] = butterflyPow2(c[4], c[6], w, wS, q, twoQ)
		c[5], c[7] = butterflyPow2(c[5], c[7], w, wS, q, twoQ)
	}

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		j := 2 * i

		c := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))

		c[0], c[1] = butterflyPow2(c[0], c[1], tw[i+N/2+0], twS[i+N/2+0], q, twoQ)
		c[2], c[3] = butterflyPow2(c[2], c[3], tw[i+N/2+1], twS[i+N/2+1], q, twoQ)
		c[4], c[5] = butterflyPow2(c[4], c[5], tw[i+N/2+2], twS[i+N/2+2], q, twoQ)
		c[6], c[7] = butterflyPow2(c[6], c[7], tw[i+N/2+3], twS[i+N/2+3], q, twoQ)
	}
}

// inttInPlacePow2Unroll computes the Inverse NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) >= 32.
func inttInPlacePow2Unroll(coeffs, twInv, twInvS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1
	var w, wS uint64

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		j := 2 * i

		c := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))

		c[0], c[1] = invButterflyPow2(c[0], c[1], twInv[i+N/2+0], twInvS[i+N/2+0], q, twoQ)
		c[2], c[3] = invButterflyPow2(c[2], c[3], twInv[i+N/2+1], twInvS[i+N/2+1], q, twoQ)
		c[4], c[5] = invButterflyPow2(c[4], c[5], twInv[i+N/2+2], twInvS[i+N/2+2], q, twoQ)
		c[6], c[7] = invButterflyPow2(c[6], c[7], twInv[i+N/2+3], twInvS[i+N/2+3], q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		j := 4 * i

		c := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))

		w, wS = twInv[i+N/4], twInvS[i+N/4]
		c[0], c[2] = invButterflyPow2(c[0], c[2], w, wS, q, twoQ)
		c[1], c[3] = invButterflyPow2(c[1], c[3], w, wS, q, twoQ)

		w, wS = twInv[i+N/4+1], twInvS[i+N/4+1]
		c[4], c[6] = invButterflyPow2(c[4], c[6], w, wS, q, twoQ)
		c[5], c[7] = invButterflyPow2(c[5], c[7], w, wS, q, twoQ)
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		w, wS = twInv[i+N/8], twInvS[i+N/8]

		j := 8 * i

		c := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))

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

			w, wS = twInv[m+i], twInvS[m+i]

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))
				c1 := (*[8]uint64)(unsafe.Pointer(&coeffs[j+t]))

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

	w, wS = twInv[1], twInvS[1]
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]uint64)(unsafe.Pointer(&coeffs[j]))
		c1 := (*[8]uint64)(unsafe.Pointer(&coeffs[j+t]))

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
