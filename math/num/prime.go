package num

import (
	"math/bits"
	"math/rand/v2"
	"slices"
)

var (
	smallPrimes = []uint64{
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
	}
)

// IsPrime checks of x is prime.
// Any x <= 1 are not considered prime.
func IsPrime[T Integer](x T) bool {
	return (x > 1) && isPrimeUint64(uint64(Abs((x))))
}

// isPrime checks of x is prime.
// 0 and 1 are not considered prime.
func isPrimeUint64(x uint64) bool {
	for _, p := range smallPrimes {
		if x == p {
			return true
		} else if x%p == 0 {
			return false
		}
	}

	s := bits.TrailingZeros64(x - 1)
	d := (x - 1) >> s

	tests := []uint64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}
	for _, a := range tests {
		n := Exp(a, d, x)
		var y uint64
		for i := 0; i < s; i++ {
			y = Mul(n, n, x)
			if y == 1 && n != 1 && n != x-1 {
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

// PrevPrime returns the largest prime of the form x - k*skip, where k >= 0.
//
// Panics if skip <= 0 or no such prime exists.
func PrevPrime[T Integer](x, skip T) T {
	if skip <= 0 {
		panic("skip must be positive")
	}

	g := GCD(x, skip)
	if g > 1 {
		if x >= g && (x-g)%skip == 0 && IsPrime(g) {
			return g
		}
		panic("no previous prime exists")
	}

	for t := x; ; t -= skip {
		if IsPrime(t) {
			return t
		}

		if t < skip {
			break
		}
	}

	panic("no previous prime exists")
}

// NextPrime returns the smallest prime of the form x + k*skip, where k >= 0.
//
// Panics if skip <= 0 or no such prime exists.
func NextPrime[T Integer](x, skip T) T {
	if skip <= 0 {
		panic("skip must be positive")
	}

	g := GCD(x, skip)
	if g > 1 {
		if x <= g && (g-x%skip)%skip == 0 && IsPrime(g) {
			return g
		}
		panic("no next prime exists")
	}

	if x < 2 {
		x = T(Reduce(x, uint64(skip)))
	}

	for t := x; ; t += skip {
		if IsPrime(t) {
			return t
		}

		if t+skip <= t {
			break
		}
	}

	panic("no next prime exists")
}

// IsProdPowerOf checks if x can be expressed as a product of the powers of given factors.
func IsProdPowerOf[T Integer](x T, factors []T) bool {
	memo := make(map[T]bool)

	var search func(T) bool
	search = func(t T) bool {
		if t == 1 {
			return true
		}

		if memo[t] {
			return false
		}

		memo[t] = true
		for _, f := range factors {
			if t == f || (f != 0 && t%f == 0 && search(t/f)) {
				return true
			}
		}
		return false
	}

	return search(x)
}

// NextProdPower returns the smallest number larger or equal than x
// that can be expressed as a product of the powers of given factors.
func NextProdPower[T Integer](x T, factors []T) T {
	best, isFound := T(1), x <= 1

	memo := make(map[T]bool)
	var search func(T)
	search = func(t T) {
		if memo[t] || (isFound && best == x) {
			return
		}

		memo[t] = true
		if t >= x && (!isFound || t < best) {
			best, isFound = t, true
		}

		if (x <= 0 && Abs(t) > max(Abs(x), 1)) || (x > 0 && isFound && Abs(t) >= Abs(best)) {
			return
		}

		for _, f := range factors {
			n := t * f
			if f != 0 && (n/f != t || (t < 0 && f < 0 && n < 0)) {
				continue
			}
			search(n)
		}
	}

	search(1)

	if !isFound {
		panic("no next product found")
	}
	return best
}

// Factor factors x. The resulting primes are sorted in ascending order.
// When x == 1, it returns [1], [1].
//
// Panics when x <= 0.
func Factor[T Integer](x T) (primes []T, exps []T) {
	if x <= 0 {
		panic("x must be positive")
	}

	factors := make(map[uint64]uint64)
	factorRecurse(uint64(x), factors)

	primes = make([]T, 0, len(factors))
	for p := range factors {
		primes = append(primes, T(p))
	}
	slices.Sort(primes)

	exps = make([]T, len(primes))
	for i, p := range primes {
		exps[i] = T(factors[uint64(p)])
	}

	return primes, exps
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

	y, c, m := rand.Uint64N(x), 1+rand.Uint64N(x-3), 1+rand.Uint64N(x-1)
	g, r, q := uint64(1), uint64(1), uint64(1)

	var t, ys uint64

	for g == 1 {
		t = y
		for i := uint64(0); i < r; i++ {
			y = Add(Mul(y, y, x), c, x)
		}
		var k uint64
		for k < r && g == 1 {
			ys = y
			for i := uint64(0); i < min(m, r-k); i++ {
				y = Add(Mul(y, y, x), c, x)
				q = Mul(q, subAbs(y, t), x)
			}
			g = GCD(q, x)
			k += m
		}
		r <<= 1
	}

	if g == x {
		for {
			ys = Add(Mul(ys, ys, x), c, x)
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

// Order returns the multiplicative order of x modulo q.
//
// Panics when q <= 1 or multiplicative order does not exist.
func Order[T Integer](x, q T) T {
	if q <= 1 {
		panic("q must be greater than 1")
	} else if GCD(x, q) != 1 {
		panic("multiplicative order does not exist")
	}

	q64 := uint64(q)
	r := Reduce(x, q64)
	ord := Totient(q64)
	primes, _ := Factor(ord)
	for _, p := range primes {
		for p > 1 && ord%p == 0 && Exp(r, ord/p, q64) == 1 {
			ord /= p
		}
	}
	return T(ord)
}

// Totient returns the Euler-Phi function of x.
//
// Panics when x <= 0.
func Totient[T Integer](x T) T {
	if x <= 0 {
		panic("x must be positive")
	}

	if x == 1 {
		return 1
	}

	primes, _ := Factor(x)
	phi := x
	for _, p := range primes {
		phi -= phi / p
	}
	return phi
}

// Generators returns the generators of the subgroup of multiplicative group modulo q.
func Generators(q uint64) []uint64 {
	if q <= 1 {
		panic("q must be greater than 1")
	}

	primes, exps := Factor(q)

	primePows := make([]uint64, len(primes))
	for i := range primePows {
		primePows[i] = Exp[*Modulus](primes[i], exps[i], nil)
	}

	subGens := make([]uint64, len(primes))
	crt := make([]uint64, len(primes))
	for i := range subGens {
		subGens[i] = primitiveRoot(primes[i], primePows[i])
		crt[i] = Mul(q/primePows[i], Inv(q/primePows[i], primePows[i]), q)
	}

	gens := make([]uint64, len(subGens))
	for i := range gens {
		for j := range crt {
			if j == i {
				gens[i] = Add(gens[i], Mul(subGens[j], crt[j], q), q)
			} else {
				gens[i] = Add(gens[i], crt[j], q)
			}
		}
	}

	if primePows[0]%8 == 0 {
		gens = append([]uint64{0}, gens...)
		gens[0] = Mul(5, crt[0], q)
		for i := 1; i < len(crt); i++ {
			gens[0] = Add(gens[0], crt[i], q)
		}
	} else if gens[0] == 1 {
		gens = gens[1:]
	}

	return gens
}

// primitiveRoot returns a generator modulo p^e.
// Returns p^e - 1 if p is 2.
func primitiveRoot(p, pExp uint64) uint64 {
	if p == 2 {
		return pExp - 1
	}

	phi := pExp - pExp/p
	primes, _ := Factor(phi)
	testPows := make([]uint64, 0, len(primes))
	for _, p := range primes {
		testPows = append(testPows, phi/p)
	}

	g := uint64(2)
	for {
		ok := true
		for _, t := range testPows {
			if Exp(g, t, pExp) == 1 {
				ok = false
				break
			}
		}
		if ok && Exp(g, phi, pExp) == 1 {
			return g
		}
		g++
	}
}

// NthRoot returns a primitive n-th root of unity modulo q,
// given the generators returned by [Generators] for q.
//
// Panics when n == 0, q <= 1, or no such root exists.
func NthRoot(n uint64, g []uint64, q uint64) uint64 {
	if n == 0 {
		panic("n must be positive")
	} else if q <= 1 {
		panic("q must be greater than 1")
	}

	primes, exps := Factor(q)

	primePows := make([]uint64, len(primes))
	ords := make([]uint64, len(primes))
	for i := range primePows {
		primePows[i] = Exp[*Modulus](primes[i], exps[i], nil)
		ords[i] = primePows[i] - primePows[i]/primes[i]
	}

	if primes[0] == 2 {
		switch n {
		case 1:
			return 1
		case 2:
			if q > 2 {
				return q - 1
			}
		}
	}

	if primePows[0]%8 == 0 {
		ords[0] = 2
		ords = append([]uint64{primePows[0] / 4}, ords...)
	} else if primePows[0] == 2 {
		ords = ords[1:]
	}

	r, ord := uint64(1), uint64(1)
	for i := range ords {
		d := GCD(n, ords[i])
		r = Mul(r, Exp(g[i], ords[i]/d, q), q)
		ord = LCM(ord, d)
	}

	if ord != n {
		panic("no n-th root of unity exists")
	}
	return r
}
