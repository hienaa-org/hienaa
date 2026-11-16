package regev

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/math/num"
)

const nttUnrollBound = 8

// fwdNTTGaussIntInPlacePow2 computes the Gaussian Integer NTT transform in-place for power-of-two length coefficients.
// twS contains the component-wise Shoup forms of tw.
// Each component of the inputs and outputs is in [0, 4q).
func fwdNTTGaussIntInPlacePow2(coeffs, tw, twS []gnum.GaussianInt, q *num.Modulus) {
	if len(coeffs) <= nttUnrollBound {
		fwdNTTGaussIntInPlacePow2Ref(coeffs, tw, twS, q)
		return
	}
	fwdNTTGaussIntInPlacePow2Unroll(coeffs, tw, twS, q)
}

// fwdButterflyGaussIntPow2 returns the Gaussian Integer Harvey butterfly.
func fwdButterflyGaussIntPow2(u, v, w, wS gnum.GaussianInt, q *num.Modulus, twoQ uint64) (gnum.GaussianInt, gnum.GaussianInt) {
	v = gnum.SMulLazy(v, w, wS, q)
	if u.Real >= twoQ {
		u.Real -= twoQ
	}
	if u.Imag >= twoQ {
		u.Imag -= twoQ
	}
	return gnum.GaussianInt{Real: u.Real + v.Real, Imag: u.Imag + v.Imag},
		gnum.GaussianInt{Real: u.Real - v.Real + twoQ, Imag: u.Imag - v.Imag + twoQ}
}

// fwdNTTGaussIntInPlacePow2Ref computes the Gaussian Integer NTT transform in-place for power-of-two length coefficients.
func fwdNTTGaussIntInPlacePow2Ref(coeffs, tw, twS []gnum.GaussianInt, q *num.Modulus) {
	N := len(coeffs)
	twoQ := q.Value() << 1

	t := N
	for m := 1; m <= N/2; m <<= 1 {
		t >>= 1
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w, wS := tw[m+i], twS[m+i]
			for j := j1; j < j2; j++ {
				coeffs[j], coeffs[j+t] = fwdButterflyGaussIntPow2(coeffs[j], coeffs[j+t], w, wS, q, twoQ)
			}
		}
	}
}

// invNTTGaussIntInPlacePow2 computes the inverse Gaussian Integer NTT transform in-place for power-of-two length coefficients.
// twInvS contains the component-wise Shoup forms of twInv.
// Each component of the inputs and outputs is in [0, 2q).
// The result is not normalized by the inverse of len(coeffs).
func invNTTGaussIntInPlacePow2(coeffs, twInv, twInvS []gnum.GaussianInt, q *num.Modulus) {
	if len(coeffs) <= nttUnrollBound {
		invNTTGaussIntInPlacePow2Ref(coeffs, twInv, twInvS, q)
		return
	}
	invNTTGaussIntInPlacePow2Unroll(coeffs, twInv, twInvS, q)
}

// invButterflyGaussIntPow2 returns the inverse Gaussian Integer Harvey butterfly.
func invButterflyGaussIntPow2(u, v, w, wS gnum.GaussianInt, q *num.Modulus, twoQ uint64) (gnum.GaussianInt, gnum.GaussianInt) {
	u, v = gnum.GaussianInt{Real: u.Real + v.Real, Imag: u.Imag + v.Imag},
		gnum.GaussianInt{Real: u.Real - v.Real + twoQ, Imag: u.Imag - v.Imag + twoQ}
	if u.Real >= twoQ {
		u.Real -= twoQ
	}
	if u.Imag >= twoQ {
		u.Imag -= twoQ
	}
	return u, gnum.SMulLazy(v, w, wS, q)
}

// invNTTGaussIntInPlacePow2Ref computes the inverse Gaussian Integer NTT transform in-place for power-of-two length coefficients.
func invNTTGaussIntInPlacePow2Ref(coeffs, twInv, twInvS []gnum.GaussianInt, q *num.Modulus) {
	N := len(coeffs)
	twoQ := q.Value() << 1

	t := 1
	for m := N / 2; m >= 1; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w, wS := twInv[m+i], twInvS[m+i]
			for j := j1; j < j2; j++ {
				coeffs[j], coeffs[j+t] = invButterflyGaussIntPow2(coeffs[j], coeffs[j+t], w, wS, q, twoQ)
			}
		}
		t <<= 1
	}
}

// fwdNTTGaussIntInPlacePow2Unroll computes the Gaussian Integer NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) > [nttUnrollBound].
func fwdNTTGaussIntInPlacePow2Unroll(coeffs, tw, twS []gnum.GaussianInt, q *num.Modulus) {
	twoQ := q.Value() << 1

	N := len(coeffs)
	L := unsafe.Sizeof(gnum.GaussianInt{})

	r := unsafe.Pointer(unsafe.SliceData(tw))
	rS := unsafe.Pointer(unsafe.SliceData(twS))
	v := unsafe.Pointer(unsafe.SliceData(coeffs))

	t := N / 2
	w := *(*gnum.GaussianInt)(unsafe.Add(r, 1*L))
	wS := *(*gnum.GaussianInt)(unsafe.Add(rS, 1*L))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j)*L))
		c1 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j+t)*L))

		c0[0], c1[0] = fwdButterflyGaussIntPow2(c0[0], c1[0], w, wS, q, twoQ)
		c0[1], c1[1] = fwdButterflyGaussIntPow2(c0[1], c1[1], w, wS, q, twoQ)
		c0[2], c1[2] = fwdButterflyGaussIntPow2(c0[2], c1[2], w, wS, q, twoQ)
		c0[3], c1[3] = fwdButterflyGaussIntPow2(c0[3], c1[3], w, wS, q, twoQ)

		c0[4], c1[4] = fwdButterflyGaussIntPow2(c0[4], c1[4], w, wS, q, twoQ)
		c0[5], c1[5] = fwdButterflyGaussIntPow2(c0[5], c1[5], w, wS, q, twoQ)
		c0[6], c1[6] = fwdButterflyGaussIntPow2(c0[6], c1[6], w, wS, q, twoQ)
		c0[7], c1[7] = fwdButterflyGaussIntPow2(c0[7], c1[7], w, wS, q, twoQ)
	}

	for m := 2; m <= N/16; m <<= 1 {
		t >>= 1
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t

			w := *(*gnum.GaussianInt)(unsafe.Add(r, uintptr(m+i)*L))
			wS := *(*gnum.GaussianInt)(unsafe.Add(rS, uintptr(m+i)*L))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j)*L))
				c1 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j+t)*L))

				c0[0], c1[0] = fwdButterflyGaussIntPow2(c0[0], c1[0], w, wS, q, twoQ)
				c0[1], c1[1] = fwdButterflyGaussIntPow2(c0[1], c1[1], w, wS, q, twoQ)
				c0[2], c1[2] = fwdButterflyGaussIntPow2(c0[2], c1[2], w, wS, q, twoQ)
				c0[3], c1[3] = fwdButterflyGaussIntPow2(c0[3], c1[3], w, wS, q, twoQ)

				c0[4], c1[4] = fwdButterflyGaussIntPow2(c0[4], c1[4], w, wS, q, twoQ)
				c0[5], c1[5] = fwdButterflyGaussIntPow2(c0[5], c1[5], w, wS, q, twoQ)
				c0[6], c1[6] = fwdButterflyGaussIntPow2(c0[6], c1[6], w, wS, q, twoQ)
				c0[7], c1[7] = fwdButterflyGaussIntPow2(c0[7], c1[7], w, wS, q, twoQ)
			}
		}
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		c := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(8*i)*L))
		w := *(*gnum.GaussianInt)(unsafe.Add(r, uintptr(i+N/8)*L))
		wS := *(*gnum.GaussianInt)(unsafe.Add(rS, uintptr(i+N/8)*L))

		c[0], c[4] = fwdButterflyGaussIntPow2(c[0], c[4], w, wS, q, twoQ)
		c[1], c[5] = fwdButterflyGaussIntPow2(c[1], c[5], w, wS, q, twoQ)
		c[2], c[6] = fwdButterflyGaussIntPow2(c[2], c[6], w, wS, q, twoQ)
		c[3], c[7] = fwdButterflyGaussIntPow2(c[3], c[7], w, wS, q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		c := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(4*i)*L))
		w := (*[2]gnum.GaussianInt)(unsafe.Add(r, uintptr(i+N/4)*L))
		wS := (*[2]gnum.GaussianInt)(unsafe.Add(rS, uintptr(i+N/4)*L))

		c[0], c[2] = fwdButterflyGaussIntPow2(c[0], c[2], w[0], wS[0], q, twoQ)
		c[1], c[3] = fwdButterflyGaussIntPow2(c[1], c[3], w[0], wS[0], q, twoQ)

		c[4], c[6] = fwdButterflyGaussIntPow2(c[4], c[6], w[1], wS[1], q, twoQ)
		c[5], c[7] = fwdButterflyGaussIntPow2(c[5], c[7], w[1], wS[1], q, twoQ)
	}

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		c := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(2*i)*L))
		w := (*[4]gnum.GaussianInt)(unsafe.Add(r, uintptr(i+N/2)*L))
		wS := (*[4]gnum.GaussianInt)(unsafe.Add(rS, uintptr(i+N/2)*L))

		c[0], c[1] = fwdButterflyGaussIntPow2(c[0], c[1], w[0], wS[0], q, twoQ)
		c[2], c[3] = fwdButterflyGaussIntPow2(c[2], c[3], w[1], wS[1], q, twoQ)
		c[4], c[5] = fwdButterflyGaussIntPow2(c[4], c[5], w[2], wS[2], q, twoQ)
		c[6], c[7] = fwdButterflyGaussIntPow2(c[6], c[7], w[3], wS[3], q, twoQ)
	}
}

// invNTTGaussIntInPlacePow2Unroll computes the inverse Gaussian Integer NTT transform in-place for power-of-two length coefficients.
// Assumes len(coeffs) > [nttUnrollBound].
func invNTTGaussIntInPlacePow2Unroll(coeffs, twInv, twInvS []gnum.GaussianInt, q *num.Modulus) {
	twoQ := q.Value() << 1

	N := len(coeffs)
	L := unsafe.Sizeof(gnum.GaussianInt{})

	r := unsafe.Pointer(unsafe.SliceData(twInv))
	rS := unsafe.Pointer(unsafe.SliceData(twInvS))
	v := unsafe.Pointer(unsafe.SliceData(coeffs))

	// t = 1, m = N / 2
	for i := 0; i < N/2; i += 4 {
		c := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(2*i)*L))
		w := (*[4]gnum.GaussianInt)(unsafe.Add(r, uintptr(i+N/2)*L))
		wS := (*[4]gnum.GaussianInt)(unsafe.Add(rS, uintptr(i+N/2)*L))

		c[0], c[1] = invButterflyGaussIntPow2(c[0], c[1], w[0], wS[0], q, twoQ)
		c[2], c[3] = invButterflyGaussIntPow2(c[2], c[3], w[1], wS[1], q, twoQ)
		c[4], c[5] = invButterflyGaussIntPow2(c[4], c[5], w[2], wS[2], q, twoQ)
		c[6], c[7] = invButterflyGaussIntPow2(c[6], c[7], w[3], wS[3], q, twoQ)
	}

	// t = 2, m = N / 4
	for i := 0; i < N/4; i += 2 {
		c := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(4*i)*L))
		w := (*[2]gnum.GaussianInt)(unsafe.Add(r, uintptr(i+N/4)*L))
		wS := (*[2]gnum.GaussianInt)(unsafe.Add(rS, uintptr(i+N/4)*L))

		c[0], c[2] = invButterflyGaussIntPow2(c[0], c[2], w[0], wS[0], q, twoQ)
		c[1], c[3] = invButterflyGaussIntPow2(c[1], c[3], w[0], wS[0], q, twoQ)

		c[4], c[6] = invButterflyGaussIntPow2(c[4], c[6], w[1], wS[1], q, twoQ)
		c[5], c[7] = invButterflyGaussIntPow2(c[5], c[7], w[1], wS[1], q, twoQ)
	}

	// t = 4, m = N / 8
	for i := 0; i < N/8; i++ {
		c := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(8*i)*L))
		w := *(*gnum.GaussianInt)(unsafe.Add(r, uintptr(i+N/8)*L))
		wS := *(*gnum.GaussianInt)(unsafe.Add(rS, uintptr(i+N/8)*L))

		c[0], c[4] = invButterflyGaussIntPow2(c[0], c[4], w, wS, q, twoQ)
		c[1], c[5] = invButterflyGaussIntPow2(c[1], c[5], w, wS, q, twoQ)
		c[2], c[6] = invButterflyGaussIntPow2(c[2], c[6], w, wS, q, twoQ)
		c[3], c[7] = invButterflyGaussIntPow2(c[3], c[7], w, wS, q, twoQ)
	}

	t := 8
	for m := N / 16; m >= 2; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t

			w := *(*gnum.GaussianInt)(unsafe.Add(r, uintptr(m+i)*L))
			wS := *(*gnum.GaussianInt)(unsafe.Add(rS, uintptr(m+i)*L))

			for j := j1; j < j2; j += 8 {
				c0 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j)*L))
				c1 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j+t)*L))

				c0[0], c1[0] = invButterflyGaussIntPow2(c0[0], c1[0], w, wS, q, twoQ)
				c0[1], c1[1] = invButterflyGaussIntPow2(c0[1], c1[1], w, wS, q, twoQ)
				c0[2], c1[2] = invButterflyGaussIntPow2(c0[2], c1[2], w, wS, q, twoQ)
				c0[3], c1[3] = invButterflyGaussIntPow2(c0[3], c1[3], w, wS, q, twoQ)

				c0[4], c1[4] = invButterflyGaussIntPow2(c0[4], c1[4], w, wS, q, twoQ)
				c0[5], c1[5] = invButterflyGaussIntPow2(c0[5], c1[5], w, wS, q, twoQ)
				c0[6], c1[6] = invButterflyGaussIntPow2(c0[6], c1[6], w, wS, q, twoQ)
				c0[7], c1[7] = invButterflyGaussIntPow2(c0[7], c1[7], w, wS, q, twoQ)
			}
		}
		t <<= 1
	}

	w := *(*gnum.GaussianInt)(unsafe.Add(r, 1*L))
	wS := *(*gnum.GaussianInt)(unsafe.Add(rS, 1*L))
	for j := 0; j < N/2; j += 8 {
		c0 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j)*L))
		c1 := (*[8]gnum.GaussianInt)(unsafe.Add(v, uintptr(j+t)*L))

		c0[0], c1[0] = invButterflyGaussIntPow2(c0[0], c1[0], w, wS, q, twoQ)
		c0[1], c1[1] = invButterflyGaussIntPow2(c0[1], c1[1], w, wS, q, twoQ)
		c0[2], c1[2] = invButterflyGaussIntPow2(c0[2], c1[2], w, wS, q, twoQ)
		c0[3], c1[3] = invButterflyGaussIntPow2(c0[3], c1[3], w, wS, q, twoQ)

		c0[4], c1[4] = invButterflyGaussIntPow2(c0[4], c1[4], w, wS, q, twoQ)
		c0[5], c1[5] = invButterflyGaussIntPow2(c0[5], c1[5], w, wS, q, twoQ)
		c0[6], c1[6] = invButterflyGaussIntPow2(c0[6], c1[6], w, wS, q, twoQ)
		c0[7], c1[7] = invButterflyGaussIntPow2(c0[7], c1[7], w, wS, q, twoQ)
	}
}
