package gnum

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

type GaussianInt struct {
	Real uint64
	Imag uint64
}

func Add(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Add(g0.Real, g1.Real, q),
		Imag: num.Add(g0.Imag, g1.Imag, q),
	}
}

func Sub(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(g0.Real, g1.Real, q),
		Imag: num.Sub(g0.Imag, g1.Imag, q),
	}
}

func Neg(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Neg(g.Real, q),
		Imag: num.Neg(g.Imag, q),
	}
}

func Mul(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(num.Mul(g0.Real, g1.Real, q), num.Mul(g0.Imag, g1.Imag, q), q),
		Imag: num.Add(num.Mul(g0.Real, g1.Imag, q), num.Mul(g0.Imag, g1.Real, q), q),
	}
}

func Exp(g GaussianInt, e uint64, q *num.Modulus) GaussianInt {
	out := GaussianInt{1, 0}
	tmp := GaussianInt{Real: g.Real, Imag: g.Imag}

	for e > 0 {
		if e&1 == 1 {
			out = Mul(out, tmp, q)
		}
		e >>= 1
		tmp = Mul(tmp, tmp, q)
	}

	return out
}

func ExpBig(g GaussianInt, e *big.Int, q *num.Modulus) GaussianInt {
	out := GaussianInt{Real: 1, Imag: 0}
	tmp := GaussianInt{Real: g.Real, Imag: g.Imag}

	exp := new(big.Int).Set(e)
	for exp.Sign() > 0 {
		if exp.Bit(0) == 1 {
			out = Mul(out, tmp, q)
		}
		exp.Rsh(exp, 1)
		tmp = Mul(tmp, tmp, q)
	}

	return out
}

func Inv(g GaussianInt, q *num.Modulus) GaussianInt {
	primes, exps := num.Factor(q.Value())
	if len(primes) != 1 {
		panic("Inv: q must be a prime power")
	}

	ord := new(big.Int).SetUint64(primes[0])
	tmp := new(big.Int).SetUint64(primes[0])
	ord.Exp(ord, big.NewInt(int64(exps[0]-1)*2), nil)
	tmp.Exp(tmp, big.NewInt(2), nil)
	tmp.Sub(tmp, big.NewInt(1))
	ord.Mul(ord, tmp)
	invExp := new(big.Int).Sub(ord, big.NewInt(1))

	return ExpBig(g, invExp, q)
}
