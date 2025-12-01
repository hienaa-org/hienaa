package pack

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/fhe/internal/gr"
	"github.com/hienaa-org/hienaa/math/csprng"
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

// gIntNthRoot returns the n-th root of unity of the multiplicative group of the Gaussian integers modulo q.
func gIntNthRoot(n int, q *num.Modulus) gnum.GaussianInt {
	rSrc := csprng.NewUniformSamplerWithSeed(nil)

	qPrimes, qExps := num.Factor(q.Value())
	if len(qPrimes) != 1 {
		panic("gIntNthRoot: q must be a prime power")
	}

	ord := new(big.Int).SetUint64(qPrimes[0])
	tmp := new(big.Int).SetUint64(qPrimes[0])
	ord.Exp(ord, big.NewInt(int64(qExps[0]-1)*2), nil)
	tmp.Exp(tmp, big.NewInt(2), nil)
	tmp.Sub(tmp, big.NewInt(1))
	ord.Mul(ord, tmp)

	primes, _ := num.Factor(n)
	checkOrd := make([]*big.Int, len(primes))
	for i := range checkOrd {
		checkOrd[i] = new(big.Int).Div(ord, big.NewInt(int64(primes[i])))
	}

	var genEl, root, genElPow gnum.GaussianInt

	rootOrd := ord.Div(ord, big.NewInt(int64(n)))
	for {
		genEl.Real = rSrc.SampleN(q.Value())
		genEl.Imag = rSrc.SampleN(q.Value())

		ok := true
		for i := range checkOrd {
			genElPow = gnum.ExpBig(genEl, checkOrd[i], q)
			if genElPow.Real == 1 && genElPow.Imag == 0 {
				ok = false
				break
			}
		}

		if ok {
			root = gnum.ExpBig(genEl, rootOrd, q)
			if !(root.Real == 1 && root.Imag == 0) {
				break
			}
		}
	}

	return root
}

// bitReverseInPlace computes the bit reverse of v in-place.
func bitReverseInPlace(v []gnum.GaussianInt) {
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
func nttGaloisRingInPlacePow2(coeffs []gnum.GaussianInt, tw []gnum.GaussianInt, q *num.Modulus) {
	t := len(coeffs) >> 1

	var u, v gnum.GaussianInt
	for m := 1; m < len(coeffs); m <<= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w := tw[m+i]
			for j := j1; j < j2; j++ {
				u = coeffs[j]
				v = gnum.Mul(coeffs[j+t], w, q)
				coeffs[j] = gnum.Add(u, v, q)
				coeffs[j+t] = gnum.Sub(u, v, q)
			}
		}
		t >>= 1
	}
}

// invNTTGaloisRingInPlacePow2 performs inverse NTT over the Galois ring in-place.
func invNTTGaloisRingInPlacePow2(coeffs []gnum.GaussianInt, twInv []gnum.GaussianInt, q *num.Modulus) {
	t := 1

	var u, v gnum.GaussianInt
	for m := len(coeffs) >> 1; m >= 1; m >>= 1 {
		for i := 0; i < m; i++ {
			j1 := i * t << 1
			j2 := j1 + t
			w := twInv[m+i]
			for j := j1; j < j2; j++ {
				u = gnum.Add(coeffs[j], coeffs[j+t], q)
				v = gnum.Sub(coeffs[j], coeffs[j+t], q)
				coeffs[j] = u
				coeffs[j+t] = gnum.Mul(v, w, q)
			}
		}
		t <<= 1
	}
}
