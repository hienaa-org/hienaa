package dft

import (
	"math/bits"
)

// nttInPlacePow5 computes the NTT transform in-place for power-of-five length coefficients.
func nttInPlacePow5(skip int, coeffs, tw, twS, root, rootS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	var quo uint64

	z1, z2, z3, z4 := root[1], root[2], root[3], root[4]
	z1S, z2S, z3S, z4S := rootS[1], rootS[2], rootS[3], rootS[4]

	t := N
	for m := 1; m <= (N/skip)/5; m *= 5 {
		t /= 5
		for i := 0; i < m; i++ {
			j1 := i * t * 5
			j2 := j1 + t

			w1, w2, w3, w4 := tw[m+i], tw[2*m+i], tw[3*m+i], tw[4*m+i]
			w1S, w2S, w3S, w4S := twS[m+i], twS[2*m+i], twS[3*m+i], twS[4*m+i]

			for j := j1; j < j2; j++ {
				u0, u1, u2, u3, u4 := coeffs[j], coeffs[j+t], coeffs[j+2*t], coeffs[j+3*t], coeffs[j+4*t]

				if u0 >= twoQ {
					u0 -= twoQ
				}

				quo, _ = bits.Mul64(u1, w1S)
				u1 = u1*w1 - quo*q

				quo, _ = bits.Mul64(u2, w2S)
				u2 = u2*w2 - quo*q

				quo, _ = bits.Mul64(u3, w3S)
				u3 = u3*w3 - quo*q

				quo, _ = bits.Mul64(u4, w4S)
				u4 = u4*w4 - quo*q

				x0 := u0 + u1
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u2
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u3
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u4

				u4 = twoQ - u4
				r0 := u0 + u4
				r1 := u1 + u4
				r2 := u2 + u4
				r3 := u3 + u4

				quo, _ = bits.Mul64(r1, z1S)
				x1 := r0 + r1*z1 - quo*q
				if x1 >= twoQ {
					x1 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z2S)
				x1 += r2*z2 - quo*q
				if x1 >= twoQ {
					x1 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z3S)
				x1 += r3*z3 - quo*q

				quo, _ = bits.Mul64(r1, z2S)
				x2 := r0 + r1*z2 - quo*q
				if x2 >= twoQ {
					x2 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z4S)
				x2 += r2*z4 - quo*q
				if x2 >= twoQ {
					x2 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z1S)
				x2 += r3*z1 - quo*q

				quo, _ = bits.Mul64(r1, z3S)
				x3 := r0 + r1*z3 - quo*q
				if x3 >= twoQ {
					x3 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z1S)
				x3 += r2*z1 - quo*q
				if x3 >= twoQ {
					x3 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z4S)
				x3 += r3*z4 - quo*q

				quo, _ = bits.Mul64(r1, z4S)
				x4 := r0 + r1*z4 - quo*q
				if x4 >= twoQ {
					x4 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z3S)
				x4 += r2*z3 - quo*q
				if x4 >= twoQ {
					x4 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z2S)
				x4 += r3*z2 - quo*q

				coeffs[j], coeffs[j+t], coeffs[j+2*t], coeffs[j+3*t], coeffs[j+4*t] = x0, x1, x2, x3, x4
			}
		}
	}
}

// inttInPlacePow5 computes the inverse NTT transform in-place for power-of-five length coefficients.
func inttInPlacePow5(skip int, coeffs, twInv, twInvS, root, rootS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	var quo uint64

	z1, z2, z3, z4 := root[1], root[2], root[3], root[4]
	z1S, z2S, z3S, z4S := rootS[1], rootS[2], rootS[3], rootS[4]

	t := skip
	for m := (N / skip) / 5; m >= 1; m /= 5 {
		for i := 0; i < m; i++ {
			j1 := i * t * 5
			j2 := j1 + t

			w1, w2, w3, w4 := twInv[m+i], twInv[2*m+i], twInv[3*m+i], twInv[4*m+i]
			w1S, w2S, w3S, w4S := twInvS[m+i], twInvS[2*m+i], twInvS[3*m+i], twInvS[4*m+i]

			for j := j1; j < j2; j++ {
				u0, u1, u2, u3, u4 := coeffs[j], coeffs[j+t], coeffs[j+2*t], coeffs[j+3*t], coeffs[j+4*t]

				x0 := u0 + u1
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u2
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u3
				if x0 >= twoQ {
					x0 -= twoQ
				}
				x0 += u4
				if x0 >= twoQ {
					x0 -= twoQ
				}

				u4 = twoQ - u4
				r0 := u0 + u4
				r1 := u1 + u4
				r2 := u2 + u4
				r3 := u3 + u4

				quo, _ = bits.Mul64(r1, z4S)
				x1 := r0 + r1*z4 - quo*q
				if x1 >= twoQ {
					x1 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z3S)
				x1 += r2*z3 - quo*q
				if x1 >= twoQ {
					x1 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z2S)
				x1 += r3*z2 - quo*q

				quo, _ = bits.Mul64(r1, z3S)
				x2 := r0 + r1*z3 - quo*q
				if x2 >= twoQ {
					x2 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z1S)
				x2 += r2*z1 - quo*q
				if x2 >= twoQ {
					x2 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z4S)
				x2 += r3*z4 - quo*q

				quo, _ = bits.Mul64(r1, z2S)
				x3 := r0 + r1*z2 - quo*q
				if x3 >= twoQ {
					x3 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z4S)
				x3 += r2*z4 - quo*q
				if x3 >= twoQ {
					x3 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z1S)
				x3 += r3*z1 - quo*q

				quo, _ = bits.Mul64(r1, z1S)
				x4 := r0 + r1*z1 - quo*q
				if x4 >= twoQ {
					x4 -= twoQ
				}
				quo, _ = bits.Mul64(r2, z2S)
				x4 += r2*z2 - quo*q
				if x4 >= twoQ {
					x4 -= twoQ
				}
				quo, _ = bits.Mul64(r3, z3S)
				x4 += r3*z3 - quo*q

				quo, _ = bits.Mul64(x1, w1S)
				x1 = x1*w1 - quo*q

				quo, _ = bits.Mul64(x2, w2S)
				x2 = x2*w2 - quo*q

				quo, _ = bits.Mul64(x3, w3S)
				x3 = x3*w3 - quo*q

				quo, _ = bits.Mul64(x4, w4S)
				x4 = x4*w4 - quo*q

				coeffs[j], coeffs[j+t], coeffs[j+2*t], coeffs[j+3*t], coeffs[j+4*t] = x0, x1, x2, x3, x4
			}
		}
		t *= 5
	}
}
