package rns

import (
	"math/bits"
)

// nttInPlacePow3 computes the NTT transform in-place for power-of-three length coefficients.
func nttInPlacePow3(skip int, coeffs, tw, twS, root, rootS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	var quo uint64

	z1, z2 := root[1], root[2]
	z1S, z2S := rootS[1], rootS[2]

	t := N
	for m := 1; m <= (N/skip)/3; m *= 3 {
		t /= 3
		for i := 0; i < m; i++ {
			j1 := i * t * 3
			j2 := j1 + t

			w1, w2 := tw[m+i], tw[2*m+i]
			w1S, w2S := twS[m+i], twS[2*m+i]

			for j := j1; j < j2; j++ {
				u0, u1, u2 := coeffs[j], coeffs[j+t], coeffs[j+2*t]

				if u0 >= twoQ {
					u0 -= twoQ
				}

				quo, _ = bits.Mul64(u1, w1S)
				u1 = u1*w1 - quo*q

				quo, _ = bits.Mul64(u2, w2S)
				u2 = u2*w2 - quo*q

				x0 := u0 + u1
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u2

				u2 = twoQ - u2
				r0 := u0 + u2
				r1 := u1 + u2

				quo, _ = bits.Mul64(r1, z1S)
				x1 := r0 + r1*z1 - quo*q

				quo, _ = bits.Mul64(r1, z2S)
				x2 := r0 + r1*z2 - quo*q

				coeffs[j], coeffs[j+t], coeffs[j+2*t] = x0, x1, x2
			}
		}
	}
}

// inttInPlacePow3 computes the inverse NTT transform in-place for power-of-three length coefficients.
func inttInPlacePow3(skip int, coeffs, twInv, twInvS, root, rootS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	var quo uint64

	z1, z2 := root[1], root[2]
	z1S, z2S := rootS[1], rootS[2]

	t := skip
	for m := (N / skip) / 3; m >= 1; m /= 3 {
		for i := 0; i < m; i++ {
			j1 := i * t * 3
			j2 := j1 + t

			w1, w2 := twInv[m+i], twInv[2*m+i]
			w1S, w2S := twInvS[m+i], twInvS[2*m+i]

			for j := j1; j < j2; j++ {
				u0, u1, u2 := coeffs[j], coeffs[j+t], coeffs[j+2*t]

				x0 := u0 + u1
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u2
				if x0 >= twoQ {
					x0 -= twoQ
				}

				u2 = twoQ - u2
				r0 := u0 + u2
				r1 := u1 + u2

				quo, _ = bits.Mul64(r1, z2S)
				x1 := r0 + r1*z2 - quo*q

				quo, _ = bits.Mul64(r1, z1S)
				x2 := r0 + r1*z1 - quo*q

				quo, _ = bits.Mul64(x1, w1S)
				x1 = x1*w1 - quo*q

				quo, _ = bits.Mul64(x2, w2S)
				x2 = x2*w2 - quo*q

				coeffs[j], coeffs[j+t], coeffs[j+2*t] = x0, x1, x2
			}
		}
		t *= 3
	}
}
