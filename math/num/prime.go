package num

import (
	"math/bits"
	"math/rand"

	"github.com/hienaa-org/hienaa/math/mod"
)

var (
	smallPrimes = []uint64{
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
	}
)

// IsPrime checks if x is prime.
// 0 and 1 are not considered prime.
func IsPrime(x uint64) bool {
	return IsPrimeModulus(mod.NewModulus(x))
}

// IsPrimeModulus checks of x is prime.
// 0 and 1 are not considered prime.
func IsPrimeModulus(x *mod.Modulus) bool {
	xv := x.Value()

	if xv == 0 || xv == 1 {
		return false
	}

	for _, p := range smallPrimes {
		if xv%p == 0 {
			return false
		}
	}

	s := bits.TrailingZeros64(xv - 1)
	d := (xv - 1) >> s

	tests := []uint64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}
	for _, a := range tests {
		n := mod.Exp(a, d, x)
		var y uint64
		for i := 0; i < s; i++ {
			y = mod.Mul(n, n, x)
			if y == 1 && n != 1 && n != xv-1 {
				return false
			}
			n = y
		}
		if y != 1 {
			return false
		}
	}

	return true
}

// NextPrime returns the next prime number of x with skip.
// If skip == 0, it panics.
// If skip < 0, it finds a previous prime.
func NextPrime(x uint64, skip uint64) uint64 {
	if skip == 0 {
		panic("NextPrime: skip must be nonzero.")
	}

	for t := x + skip; ; t += skip {
		if IsPrime(t) {
			return t
		}
	}
}

// IsProdPowerOf checks if x can be expressed as a product of the powers of given factors.
func IsProdPowerOf(x uint64, factors []uint64) bool {
	for _, f := range factors {
		if f == 0 {
			continue
		}
		for x%f == 0 {
			x /= f
		}
	}
	return x == 1
}

// NextProdPower returns the next number of x that can be expressed as
// a product of the powers of given factors.
func NextProdPower(x uint64, factors []uint64) uint64 {
	xNext := x + 1
	for !IsProdPowerOf(xNext, factors) {
		xNext++
	}
	return xNext
}

// Factor factors x.
func Factor(x uint64) map[uint64]uint64 {
	factors := make(map[uint64]uint64)
	factorRecurse(x, factors)
	return factors
}

// factorRecurse finds a non-trivial factor of x and adds it to factors.
// It uses Brent-Pollard's rho algorithm without trivial checks.
func factorRecurse(x uint64, factors map[uint64]uint64) {
	for _, p := range smallPrimes {
		for {
			if x%p != 0 {
				break
			}
			factors[p] += 1
			x /= p
		}
	}

	switch {
	case x == 0:
		return
	case x == 1:
		if len(factors) == 0 {
			factors[1] = 1
		}
		return
	case IsPrime(x):
		factors[x] += 1
		return
	}

	n := mod.NewModulus(x)
	y, c, m := randUint64n(x), 1+randUint64n(x-3), randUint64n(x)
	g, r, q := uint64(1), uint64(1), uint64(1)

	var t, ys uint64

	for g == 1 {
		t = y
		for i := uint64(0); i < r; i++ {
			y = mod.Add(mod.Mul(y, y, n), c, n)
		}
		var k uint64
		for k < r && g == 1 {
			ys = y
			for i := uint64(0); i < min(m, r-k); i++ {
				y = mod.Add(mod.Mul(y, y, n), c, n)
				q = mod.Mul(q, subAbs(y, t), n)
			}
			g = GCD(q, x)
			k += m
		}
		r <<= 1
	}

	if g == x {
		for {
			ys = mod.Add(mod.Mul(ys, ys, n), c, n)
			g = GCD(subAbs(ys, t), x)

			if g > 1 {
				break
			}
		}
	}

	factorRecurse(g, factors)
	factorRecurse(x/g, factors)
}

// subAbs returns |x - y|.
func subAbs(x, y uint64) uint64 {
	if x > y {
		return x - y
	}
	return y - x
}

// randUint64n returns random uint64 in [0, n).
func randUint64n(n uint64) uint64 {
	return uint64(rand.Int63n(int64(n)))
}
