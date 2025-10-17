package pack

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/gr"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// quoRem computes the quotient and remainder of two polynomials modulo a modulus.
func quoRem(p0, p1 []uint64, mod *num.Modulus) ([]uint64, []uint64) {
	switch {
	case len(p0) < len(p1):
		panic("quotient: dividend is shorter than divisor")
	case num.GCD(mod.Value(), p1[len(p1)-1]) != 1:
		panic("quotient: divisor is not coprime with modulus")
	}

	quo := make([]uint64, len(p0)-len(p1)+1)
	rem := make([]uint64, len(p0))
	copy(rem, p0)

	lcInv := num.Inv(p1[len(p1)-1], mod)
	for i := 0; i <= len(p0)-len(p1); i++ {
		if rem[len(rem)-i-1] != 0 {
			quo[len(quo)-i-1] = num.Mul(rem[len(rem)-i-1], lcInv, mod)
			vec.ScalarMulSubTo(rem[len(rem)-i-len(p1):len(rem)-i], p1, quo[len(quo)-i-1], mod)
		}
	}

	return quo, rem[:len(p1)-1]
}

// grNthRoot returns the n-th root of unity of Galois ring elements.
func grNthRoot(r *gr.GaloisRing, n int) *gr.Element {
	rSrc := csprng.NewUniformSamplerWithSeed(nil)

	primes, _ := num.Factor(n)
	checkOrd := make([]*big.Int, len(primes))
	for i := range checkOrd {
		checkOrd[i] = new(big.Int).Div(r.Ord(), big.NewInt(int64(primes[i])))
	}

	genEl := r.NewElement()
	root := r.NewElement()
	one := r.NewElementFromUint64(1)
	rootOrd := new(big.Int).Div(r.Ord(), big.NewInt(int64(n)))
	genElPow := r.NewElement()
	genCoeffs := genEl.Coeffs()
	for {
		for i := 0; i < r.Rank(); i++ {
			genCoeffs[i] = rSrc.SampleN(r.Modulus())
		}

		ok := true
		for i := range checkOrd {
			r.ExpBigTo(genElPow, genEl, checkOrd[i])
			if genElPow.IsEqual(one) {
				ok = false
				break
			}
		}

		if ok {
			r.ExpBigTo(root, genEl, rootOrd)
			if !root.IsEqual(one) {
				break
			}
		}
	}

	return root
}

// bitReverseInPlace computes the bit reverse of v in-place.
func bitReverseInPlace(v []*gr.Element) {
	var bit, j int
	for i := 1; i < len(v); i++ {
		bit = len(v) >> 1
		for j >= bit {
			j -= bit
			bit >>= 1
		}
		j += bit
		if i < j {
			v[i], v[j] = v[j], v[i]
		}
	}
}

// nttGaloisRingInPlacePow2 performs NTT over the Galois ring in-place.
func nttGaloisRingInPlacePow2(coeffs []*gr.Element, tw []*gr.Element, r *gr.GaloisRing) {
	t := len(coeffs) >> 1
	u := r.NewElement()
	v := r.NewElement()
	for m := 1; m < len(coeffs); m <<= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w := tw[m+i]
			for j := j1; j < j2; j++ {
				u.CopyFrom(coeffs[j])
				r.MulTo(v, coeffs[j+t], w)
				r.AddTo(coeffs[j], u, v)
				r.SubTo(coeffs[j+t], u, v)
			}
		}
		t >>= 1
	}
}

// invNTTGaloisRingInPlacePow2 performs inverse NTT over the Galois ring in-place.
func invNTTGaloisRingInPlacePow2(coeffs []*gr.Element, twInv []*gr.Element, r *gr.GaloisRing) {
	t := 1
	u := r.NewElement()
	v := r.NewElement()
	for m := len(coeffs) >> 1; m >= 1; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w := twInv[m+i]
			for j := j1; j < j2; j++ {
				r.AddTo(u, coeffs[j], coeffs[j+t])
				r.SubTo(v, coeffs[j], coeffs[j+t])
				r.MulTo(v, v, w)
				coeffs[j].CopyFrom(u)
				coeffs[j+t].CopyFrom(v)
			}
		}
		t <<= 1
	}
}
