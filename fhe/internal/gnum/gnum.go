package gnum

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

// GaussianInt represents a Gaussian integer.
type GaussianInt struct {
	Real uint64
	Imag uint64
}

// Add returns g0 + g1.
func Add(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Add(g0.Real, g1.Real, q),
		Imag: num.Add(g0.Imag, g1.Imag, q),
	}
}

// Sub returns g0 - g1.
func Sub(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(g0.Real, g1.Real, q),
		Imag: num.Sub(g0.Imag, g1.Imag, q),
	}
}

// Neg returns -g.
func Neg(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Neg(g.Real, q),
		Imag: num.Neg(g.Imag, q),
	}
}

// Mul returns g0 * g1.
func Mul(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(num.Mul(g0.Real, g1.Real, q), num.Mul(g0.Imag, g1.Imag, q), q),
		Imag: num.Add(num.Mul(g0.Real, g1.Imag, q), num.Mul(g0.Imag, g1.Real, q), q),
	}
}

// Exp returns g^e.
func Exp(g GaussianInt, e uint64, q *num.Modulus) GaussianInt {
	out := GaussianInt{1, 0}
	r := GaussianInt{Real: g.Real, Imag: g.Imag}

	for e > 0 {
		if e&1 == 1 {
			out = Mul(out, r, q)
		}
		e >>= 1
		r = Mul(r, r, q)
	}

	return out
}

// ExpBig returns g^e where e is a [*big.Int].
func ExpBig(g GaussianInt, e *big.Int, q *num.Modulus) GaussianInt {
	out := GaussianInt{Real: 1, Imag: 0}
	r := GaussianInt{Real: g.Real, Imag: g.Imag}

	exp := new(big.Int).Set(e)
	for exp.Sign() > 0 {
		if exp.Bit(0) == 1 {
			out = Mul(out, r, q)
		}
		exp.Rsh(exp, 1)
		r = Mul(r, r, q)
	}

	return out
}

// Inv returns the multiplicative inverse of g mod q.
func Inv(g GaussianInt, q *num.Modulus) GaussianInt {
	primes, exps := num.Factor(q.Value())
	if len(primes) != 1 {
		panic("q must be a prime power")
	}

	ord := new(big.Int).SetUint64(primes[0])
	r := new(big.Int).SetUint64(primes[0])
	ord.Exp(ord, big.NewInt(int64(exps[0]-1)*2), nil)
	r.Exp(r, big.NewInt(2), nil)
	r.Sub(r, big.NewInt(1))
	ord.Mul(ord, r)
	invExp := new(big.Int).Sub(ord, big.NewInt(1))

	return ExpBig(g, invExp, q)
}
