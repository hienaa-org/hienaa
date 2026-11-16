// Package gnum implements various utility functions for modular and number arithmetic
// over the Gaussian Integers.
package gnum

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

// GaussianInt represents a Gaussian Integer.
type GaussianInt struct {
	Real uint64
	Imag uint64
}

// Add returns g0 + g1 mod q.
func Add(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Add(g0.Real, g1.Real, q),
		Imag: num.Add(g0.Imag, g1.Imag, q),
	}
}

// Sub returns g0 - g1 mod q.
func Sub(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(g0.Real, g1.Real, q),
		Imag: num.Sub(g0.Imag, g1.Imag, q),
	}
}

// Neg returns -g mod q.
func Neg(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Neg(g.Real, q),
		Imag: num.Neg(g.Imag, q),
	}
}

// Mul returns g0 * g1 mod q.
func Mul(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(num.Mul(g0.Real, g1.Real, q), num.Mul(g0.Imag, g1.Imag, q), q),
		Imag: num.Add(num.Mul(g0.Real, g1.Imag, q), num.Mul(g0.Imag, g1.Real, q), q),
	}
}

// MulLazy returns g0 * g1 mod q,
// but each component of the result is in [0, 2q).
func MulLazy(g0, g1 GaussianInt, q *num.Modulus) GaussianInt {
	twoQ := q.Value() << 1
	real := num.MulLazy(g0.Real, g1.Real, q) + twoQ - num.MulLazy(g0.Imag, g1.Imag, q)
	imag := num.MulLazy(g0.Real, g1.Imag, q) + num.MulLazy(g0.Imag, g1.Real, q)
	if real >= twoQ {
		real -= twoQ
	}
	if imag >= twoQ {
		imag -= twoQ
	}
	return GaussianInt{Real: real, Imag: imag}
}

// MForm transforms each component of g into Montgomery form.
func MForm(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.MForm(g.Real, q),
		Imag: num.MForm(g.Imag, q),
	}
}

// InvMForm transforms each component of gM to Standard form.
func InvMForm(gM GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.InvMForm(gM.Real, q),
		Imag: num.InvMForm(gM.Imag, q),
	}
}

// MMul returns g0 * g1 mod q in Montgomery form.
// g0M and g1M must be a valid Montgomery form.
func MMul(g0M, g1M GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(num.MMul(g0M.Real, g1M.Real, q), num.MMul(g0M.Imag, g1M.Imag, q), q),
		Imag: num.Add(num.MMul(g0M.Real, g1M.Imag, q), num.MMul(g0M.Imag, g1M.Real, q), q),
	}
}

// MMulLazy returns g0 * g1 mod q in Montgomery form,
// but each component of the result is in [0, 2q).
// g0M and g1M must be a valid Montgomery form.
func MMulLazy(g0M, g1M GaussianInt, q *num.Modulus) GaussianInt {
	twoQ := q.Value() << 1
	real := num.MMulLazy(g0M.Real, g1M.Real, q) + twoQ - num.MMulLazy(g0M.Imag, g1M.Imag, q)
	imag := num.MMulLazy(g0M.Real, g1M.Imag, q) + num.MMulLazy(g0M.Imag, g1M.Real, q)
	if real >= twoQ {
		real -= twoQ
	}
	if imag >= twoQ {
		imag -= twoQ
	}
	return GaussianInt{Real: real, Imag: imag}
}

// SForm transforms each component of g into Shoup form.
func SForm(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.SForm(g.Real, q),
		Imag: num.SForm(g.Imag, q),
	}
}

// SMul returns g0 * g1 mod q using Shoup multiplication.
// g1S must be a valid Shoup form.
func SMul(g0, g1, g1S GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Sub(num.SMul(g0.Real, g1.Real, g1S.Real, q), num.SMul(g0.Imag, g1.Imag, g1S.Imag, q), q),
		Imag: num.Add(num.SMul(g0.Real, g1.Imag, g1S.Imag, q), num.SMul(g0.Imag, g1.Real, g1S.Real, q), q),
	}
}

// SMulLazy returns g0 * g1 mod q using Shoup multiplication,
// but each component of the result is in [0, 2q).
// g1S must be a valid Shoup form.
func SMulLazy(g0, g1, g1S GaussianInt, q *num.Modulus) GaussianInt {
	twoQ := q.Value() << 1
	real := num.SMulLazy(g0.Real, g1.Real, g1S.Real, q) + twoQ - num.SMulLazy(g0.Imag, g1.Imag, g1S.Imag, q)
	imag := num.SMulLazy(g0.Real, g1.Imag, g1S.Imag, q) + num.SMulLazy(g0.Imag, g1.Real, g1S.Real, q)
	if real >= twoQ {
		real -= twoQ
	}
	if imag >= twoQ {
		imag -= twoQ
	}
	return GaussianInt{Real: real, Imag: imag}
}

// Reduce returns g mod q using Barrett reduction on each component.
func Reduce(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Reduce(g.Real, q),
		Imag: num.Reduce(g.Imag, q),
	}
}

// Reduce128 returns g mod q.
// gHi and gLo contain the high and low words of each component of g.
func Reduce128(gHi, gLo GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Reduce128(gHi.Real, gLo.Real, q),
		Imag: num.Reduce128(gHi.Imag, gLo.Imag, q),
	}
}

// Reduce128Lazy returns g mod q,
// but each component of the result is in [0, 2q).
// gHi and gLo contain the high and low words of each component of g.
func Reduce128Lazy(gHi, gLo GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Reduce128Lazy(gHi.Real, gLo.Real, q),
		Imag: num.Reduce128Lazy(gHi.Imag, gLo.Imag, q),
	}
}

// Reduce2Q reduces g assuming each component is in [0, 2q).
func Reduce2Q(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Reduce2Q(g.Real, q),
		Imag: num.Reduce2Q(g.Imag, q),
	}
}

// Reduce4Q reduces g assuming each component is in [0, 4q).
func Reduce4Q(g GaussianInt, q *num.Modulus) GaussianInt {
	return GaussianInt{
		Real: num.Reduce4Q(g.Real, q),
		Imag: num.Reduce4Q(g.Imag, q),
	}
}

// Exp returns g^e mod q.
//
// Panics if e is negative or q is nil.
func Exp[T uint64 | *big.Int](g GaussianInt, e T, q *num.Modulus) GaussianInt {
	e64, eIsUint64 := any(e).(uint64)
	eBig, _ := any(e).(*big.Int)

	out := GaussianInt{Real: 1, Imag: 0}
	r := GaussianInt{Real: g.Real, Imag: g.Imag}
	switch {
	case eIsUint64:
		for e64 > 0 {
			if e64&1 == 1 {
				out = Mul(out, r, q)
			}
			e64 >>= 1
			if e64 == 0 {
				break
			}
			r = Mul(r, r, q)
		}

	case !eIsUint64:
		if eBig.Sign() < 0 {
			panic("e must be nonnegative")
		}

		exp := new(big.Int).Set(eBig)
		for exp.Sign() > 0 {
			if exp.Bit(0) == 1 {
				out = Mul(out, r, q)
			}
			exp.Rsh(exp, 1)
			if exp.Sign() == 0 {
				break
			}
			r = Mul(r, r, q)
		}
	}

	return out
}

// Inv returns the multiplicative inverse of g mod q.
//
// Panics if no inverse exists.
func Inv(g GaussianInt, q *num.Modulus) GaussianInt {
	norm := num.Add(num.Mul(g.Real, g.Real, q), num.Mul(g.Imag, g.Imag, q), q)
	normInv := num.Inv(norm, q)
	return GaussianInt{
		Real: num.Mul(g.Real, normInv, q),
		Imag: num.Neg(num.Mul(g.Imag, normInv, q), q),
	}
}
