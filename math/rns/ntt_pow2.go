package rns

import (
	"math/bits"
)

// nttInPlacePow2 computes the NTT transform in-place for power-of-two length coefficients.
func nttInPlacePow2(coeffs, tw, twS []uint64, q uint64) {
	if len(coeffs) < 16 {
		nttInPlacePow2Ref(coeffs, tw, twS, q)
		return
	}
	nttInPlacePow2Deg16(coeffs, tw, twS, q)
}

// butterflyNoCmp returns the Harvey butterfly without reduction.
func butterflyNoCmp(u, v, w, wS, q, twoQ uint64) (uint64, uint64) {
	quo, _ := bits.Mul64(v, wS)
	t := v*w - quo*q
	return u + t, u - t + twoQ
}

// butterfly returns the Harvey butterfly.
func butterfly(u, v, w, wS, q, twoQ uint64) (uint64, uint64) {
	if u >= twoQ {
		u -= twoQ
	}
	quo, _ := bits.Mul64(v, wS)
	t := v*w - quo*q
	return u + t, u - t + twoQ
}

// nttInPlacePow2Ref computes the NTT transform in-place for power-of-two length coefficients.
func nttInPlacePow2Ref(coeffs, tw, twS []uint64, q uint64) {
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
				coeffs[j], coeffs[j+t] = butterfly(coeffs[j], coeffs[j+t], w, wS, q, twoQ)
			}
		}
	}
}

// invButterflyNoCmp returns the inverse Harvey butterfly without reduction.
func invButterflyNoCmp(u, v, w, wS, q, twoQ uint64) (uint64, uint64) {
	u, v = u+v, u-v+twoQ
	quo, _ := bits.Mul64(v, wS)
	return u, v*w - quo*q
}

// invButterfly returns the inverse Harvey butterfly.
func invButterfly(u, v, w, wS, q, twoQ uint64) (uint64, uint64) {
	u, v = u+v, u-v+twoQ
	if u >= twoQ {
		u -= twoQ
	}
	quo, _ := bits.Mul64(v, wS)
	return u, v*w - quo*q
}

// inttInPlacePow2 computes the inverse NTT transform in-place for power-of-two length coefficients.
func inttInPlacePow2(coeffs, twInv, twInvS []uint64, q uint64) {
	if len(coeffs) < 16 {
		inttInPlacePow2Ref(coeffs, twInv, twInvS, q)
		return
	}
	inttInPlacePow2Deg16(coeffs, twInv, twInvS, q)
}

// inttInPlacePow2Ref computes the inverse NTT transform in-place for power-of-two length coefficients.
func inttInPlacePow2Ref(coeffs, twInv, twInvS []uint64, q uint64) {
	N := len(coeffs)
	twoQ := q << 1

	t := 1
	for m := N / 2; m >= 1; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w, wS := twInv[m+i], twInvS[m+i]
			for j := j1; j < j2; j++ {
				coeffs[j], coeffs[j+t] = invButterfly(coeffs[j], coeffs[j+t], w, wS, q, twoQ)
			}
		}
		t <<= 1
	}
}
