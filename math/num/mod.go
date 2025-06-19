package num

import (
	"github.com/hienaa-org/hienaa/math/mod"
)

// Order returns the multiplicative order of x modulo q.
func Order(x uint64, q *mod.Modulus) uint64 {
	ord := uint64(1)
	acc := mod.Reduce(x, q)
	for acc != 1 {
		acc = mod.Mul(acc, x, q)
		ord += 1
	}
	return ord
}

// EulerPhi returns the Euler-Phi function of x.
func EulerPhi(x uint64) uint64 {
	if x == 0 || x == 1 {
		return x
	}

	return eulerPhiWithFactors(x, Factor(x))
}

// eulerPhiWithFactors returns the Euler-Phi function of x, given its prime factors.
func eulerPhiWithFactors(x uint64, factors map[uint64]uint64) uint64 {
	phi := x
	for f := range factors {
		phi -= phi / f
	}
	return phi
}

// PrimitiveRoot returns a primitive root of q.
func PrimitiveRoot(q *mod.Modulus) uint64 {
	factors := Factor(q.Value())
	factorPows := make([]uint64, 0, len(factors))
	for p, e := range factors {
		factorPows = append(factorPows, Exp(p, e))
	}

	if len(factorPows) == 1 {
		phiQ := eulerPhiWithFactors(q.Value(), factors)
		phiQFactors := Factor(phiQ)
		testPows := make([]uint64, 0, len(phiQFactors))
		for f := range phiQFactors {
			testPows = append(testPows, phiQ/f)
		}

		g := uint64(2)
		for {
			ok := true
			for _, t := range testPows {
				if mod.Exp(g, t, q) == 1 {
					ok = false
					break
				}
			}
			if ok && mod.Exp(g, phiQ, q) == 1 {
				return g
			}
			g++
		}
	}

	factorPowsMod := make([]*mod.Modulus, len(factorPows))
	gFactors := make([]uint64, len(factorPows))
	for i, pExp := range factorPows {
		factorPowsMod[i] = mod.NewModulus(pExp)
		gFactors[i] = PrimitiveRoot(factorPowsMod[i])
	}

	g := uint64(0)
	for i, pExp := range factorPowsMod {
		t := q.Value() / pExp.Value()
		tInv := mod.Inv(t, pExp)

		g = mod.Add(g, mod.Mul(mod.Mul(gFactors[i], tInv, q), t, q), q)
	}
	return g
}

// NthRoot returns the N-th root of unity modulo q.
func NthRoot(n int, g uint64, q *mod.Modulus) uint64 {
	factors := Factor(q.Value())
	prime := make([]uint64, 0, len(factors))
	exp := make([]uint64, 0, len(factors))
	factorPowsMod := make([]*mod.Modulus, 0, len(factors))
	for p, e := range factors {
		prime = append(prime, p)
		exp = append(exp, e)
		factorPowsMod = append(factorPowsMod, mod.NewModulus(Exp(p, e)))
	}

	r := uint64(0)
	for i, pExp := range factorPowsMod {
		h := mod.Exp(mod.Reduce(g, pExp), (Exp(prime[i], exp[i])-Exp(prime[i], exp[i]-1))/uint64(n), pExp)
		t := q.Value() / pExp.Value()
		tInv := mod.Inv(t, pExp)

		r = mod.Add(r, mod.Mul(mod.Mul(h, tInv, q), t, q), q)
	}
	return r
}
