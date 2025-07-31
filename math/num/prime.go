package num

import (
	"math/bits"
	"math/rand"
)

var (
	smallPrimes = []uint64{
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
	}
)

// IsPrime checks if x is prime.
// 0 and 1 are not considered prime.
func IsPrime(x uint64) bool {
	return IsPrimeModulus(NewModulus(x))
}

// IsPrimeModulus checks of x is prime.
// 0 and 1 are not considered prime.
func IsPrimeModulus(x *Modulus) bool {
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
		n := Exp(a, d, x)
		var y uint64
		for i := 0; i < s; i++ {
			y = Mul(n, n, x)
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

// PrevPrime returns the previous prime number of x with skip.
// If skip == 0, or there is no prime number meets the condition, it panics.
func PrevPrime(x uint64, skip uint64) uint64 {
	if skip == 0 {
		panic("PrevPrime: skip must be nonzero.")
	}

	for t := x - skip; ; t -= skip {
		if t > MaxModulus {
			panic("PrevPrime: t must be less than or equal to MaxModulus.")
		} else if IsPrime(t) {
			return t
		}
	}
}

// NextPrime returns the next prime number of x with skip.
// If skip == 0, or there is no prime number meets the condition, it panics.
func NextPrime(x uint64, skip uint64) uint64 {
	if skip == 0 {
		panic("NextPrime: skip must be nonzero.")
	}

	for t := x + skip; ; t += skip {
		if t > MaxModulus {
			panic("NextPrime: t must be less than or equal to MaxModulus.")
		} else if IsPrime(t) {
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

	n := NewModulus(x)
	y, c, m := randUint64n(x), 1+randUint64n(x-3), randUint64n(x)
	g, r, q := uint64(1), uint64(1), uint64(1)

	var t, ys uint64

	for g == 1 {
		t = y
		for i := uint64(0); i < r; i++ {
			y = Add(Mul(y, y, n), c, n)
		}
		var k uint64
		for k < r && g == 1 {
			ys = y
			for i := uint64(0); i < min(m, r-k); i++ {
				y = Add(Mul(y, y, n), c, n)
				q = Mul(q, subAbs(y, t), n)
			}
			g = GCD(q, x)
			k += m
		}
		r <<= 1
	}

	if g == x {
		for {
			ys = Add(Mul(ys, ys, n), c, n)
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

// Order returns the multiplicative order of x modulo q.
func Order(x uint64, q *Modulus) uint64 {
	ord := uint64(1)
	acc := Reduce(x, q)
	for acc != 1 {
		acc = Mul(acc, x, q)
		ord += 1
	}
	return ord
}

// Totient returns the Euler-Phi function of x.
func Totient(x uint64) uint64 {
	if x == 0 || x == 1 {
		return x
	}

	return totientWithFactors(x, Factor(x))
}

// totientWithFactors returns the Euler-Phi function of x, given its prime factors.
func totientWithFactors(x uint64, factors map[uint64]uint64) uint64 {
	phi := x
	for f := range factors {
		phi -= phi / f
	}
	return phi
}

// PrimitiveRoot returns a primitive root of q.
func PrimitiveRoot(q *Modulus) uint64 {
	factors := Factor(q.Value())
	factorPows := make([]uint64, 0, len(factors))
	for p, e := range factors {
		pExp := uint64(1)
		for i := uint64(0); i < e; i++ {
			pExp *= p
		}
		factorPows = append(factorPows, pExp)
	}

	if len(factorPows) == 1 {
		phiQ := totientWithFactors(q.Value(), factors)
		phiQFactors := Factor(phiQ)
		testPows := make([]uint64, 0, len(phiQFactors))
		for f := range phiQFactors {
			testPows = append(testPows, phiQ/f)
		}

		g := uint64(2)
		for {
			ok := true
			for _, t := range testPows {
				if Exp(g, t, q) == 1 {
					ok = false
					break
				}
			}
			if ok && Exp(g, phiQ, q) == 1 {
				return g
			}
			g++
		}
	}

	factorPowsMod := make([]*Modulus, len(factorPows))
	gFactors := make([]uint64, len(factorPows))
	for i, pExp := range factorPows {
		factorPowsMod[i] = NewModulus(pExp)
		gFactors[i] = PrimitiveRoot(factorPowsMod[i])
	}

	g := uint64(0)
	for i, pExp := range factorPowsMod {
		t := q.Value() / pExp.Value()
		tInv := Inv(t, pExp)

		g = Add(g, Mul(Mul(gFactors[i], tInv, q), t, q), q)
	}
	return g
}

// NthRoot returns the N-th root of unity modulo q.
func NthRoot(n int, g uint64, q *Modulus) uint64 {
	factors := Factor(q.Value())
	prime := make([]uint64, 0, len(factors))
	factorPowsMod := make([]*Modulus, 0, len(factors))
	for p, e := range factors {
		prime = append(prime, p)
		pExp := uint64(1)
		for i := uint64(0); i < e; i++ {
			pExp *= p
		}
		factorPowsMod = append(factorPowsMod, NewModulus(pExp))
	}

	r := uint64(0)
	for i, pExp := range factorPowsMod {
		h := Exp(Reduce(g, pExp), (pExp.Value()-pExp.Value()/prime[i])/uint64(n), pExp)
		t := q.Value() / pExp.Value()
		tInv := Inv(t, pExp)

		r = Add(r, Mul(Mul(h, tInv, q), t, q), q)
	}
	return r
}
