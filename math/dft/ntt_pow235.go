package dft

import (
	"math/bits"
)

// fwdNTTInPlacePow2 computes the NTT transform in-place for power-of-two length coefficients.
func fwdNTTInPlacePow2(coeffs, tw, twS []uint64, q uint64) {
	if len(coeffs) < 32 {
		fwdNTTInPlacePow2Ref(coeffs, tw, twS, q)
		return
	}
	fwdNTTInPlacePow2Unroll(coeffs, tw, twS, q)
}

// fwdButterflyPow2 returns the Harvey butterfly.
func fwdButterflyPow2(u, v, w, wS, q, twoQ uint64) (uint64, uint64) {
	quo, _ := bits.Mul64(v, wS)
	t := v*w - quo*q
	if u >= twoQ {
		u -= twoQ
	}
	return u + t, u - t + twoQ
}

// fwdNTTInPlacePow2Ref computes the NTT transform in-place for power-of-two length coefficients.
func fwdNTTInPlacePow2Ref(coeffs, tw, twS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	t := N
	for m := 1; m <= N/2; m <<= 1 {
		t >>= 1
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w, wS := tw[m+i], twS[m+i]
			for j := j1; j < j2; j++ {
				coeffs[j], coeffs[j+t] = fwdButterflyPow2(coeffs[j], coeffs[j+t], w, wS, q, twoQ)
			}
		}
	}
}

// invNTTInPlacePow2 computes the inverse NTT transform in-place for power-of-two length coefficients.
func invNTTInPlacePow2(coeffs, twInv, twInvS []uint64, q uint64) {
	if len(coeffs) < 32 {
		invNTTInPlacePow2Ref(coeffs, twInv, twInvS, q)
		return
	}
	invNTTInPlacePow2Unroll(coeffs, twInv, twInvS, q)
}

// invButterflyPow2 returns the inverse Harvey butterfly.
func invButterflyPow2(u, v, w, wS, q, twoQ uint64) (uint64, uint64) {
	u, v = u+v, u-v+twoQ
	if u >= twoQ {
		u -= twoQ
	}
	quo, _ := bits.Mul64(v, wS)
	return u, v*w - quo*q
}

// invNTTInPlacePow2Ref computes the inverse NTT transform in-place for power-of-two length coefficients.
func invNTTInPlacePow2Ref(coeffs, twInv, twInvS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	t := 1
	for m := N / 2; m >= 1; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w, wS := twInv[m+i], twInvS[m+i]
			for j := j1; j < j2; j++ {
				coeffs[j], coeffs[j+t] = invButterflyPow2(coeffs[j], coeffs[j+t], w, wS, q, twoQ)
			}
		}
		t <<= 1
	}
}

// fwdNTTInPlacePow3 computes the NTT transform in-place for power-of-three length coefficients.
func fwdNTTInPlacePow3(skip int, coeffs, tw, twS, root, rootS []uint64, q uint64) {
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

				if r0 >= twoQ {
					r0 -= twoQ
				}

				quo, _ = bits.Mul64(r1, z1S)
				x1 := r0 + r1*z1 - quo*q

				quo, _ = bits.Mul64(r1, z2S)
				x2 := r0 + r1*z2 - quo*q

				coeffs[j], coeffs[j+t], coeffs[j+2*t] = x0, x1, x2
			}
		}
	}
}

// invNTTInPlacePow3 computes the inverse NTT transform in-place for power-of-three length coefficients.
func invNTTInPlacePow3(skip int, coeffs, twInv, twInvS, root, rootS []uint64, q uint64) {
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

				if r0 >= twoQ {
					r0 -= twoQ
				}

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

// fwdNTTInPlacePow5 computes the NTT transform in-place for power-of-five length coefficients.
func fwdNTTInPlacePow5(skip int, coeffs, tw, twS, root, rootS []uint64, q uint64) {
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

				if r0 >= twoQ {
					r0 -= twoQ
				}

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

// invNTTInPlacePow5 computes the inverse NTT transform in-place for power-of-five length coefficients.
func invNTTInPlacePow5(skip int, coeffs, twInv, twInvS, root, rootS []uint64, q uint64) {
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

				if r0 >= twoQ {
					r0 -= twoQ
				}

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
