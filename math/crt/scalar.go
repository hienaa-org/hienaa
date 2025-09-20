package crt

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

// Scalar represents a scalar in CRT basis.
type Scalar []uint64

// NewScalar creates a new [Scalar].
func NewScalar[T int64 | uint64 | *big.Int](x T, mod []*num.Modulus) Scalar {
	r := make(Scalar, len(mod))

	var z T
	switch any(z).(type) {
	case *big.Int:
		u := any(x).(*big.Int)
		q, t := new(big.Int), new(big.Int)
		for i := range r {
			q.SetUint64(mod[i].Value())
			r[i] = t.Mod(u, q).Uint64()
		}

	case int64:
		u := any(x).(int64)
		for i := range r {
			r[i] = num.Reduce(u, mod[i])
		}

	case uint64:
		u := any(x).(uint64)
		for i := range r {
			r[i] = num.Reduce(u, mod[i])
		}
	}

	return r
}

// isTernaryScalarToOperable checks if inputs are consistent.
func isScalarToOperable(modLen int, x ...Scalar) bool {
	for i := range x {
		if len(x[i]) != modLen {
			return false
		}
	}
	return true
}

// AddScalar returns x0 + x1.
func AddScalar(x0, x1 Scalar, mod []*num.Modulus) Scalar {
	xOut := make(Scalar, len(mod))
	AddScalarTo(xOut, x0, x1, mod)
	return xOut
}

// AddScalarTo computes xOut = x0 + x1.
func AddScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
	if !isScalarToOperable(len(mod), xOut, x0, x1) {
		panic("AddScalarTo: inputs not consistent")
	}

	for i := range mod {
		xOut[i] = num.Add(x0[i], x1[i], mod[i])
	}
}

// SubScalar returns x0 - x1.
func SubScalar(x0, x1 Scalar, mod []*num.Modulus) Scalar {
	xOut := make(Scalar, len(mod))
	SubScalarTo(xOut, x0, x1, mod)
	return xOut
}

// SubScalarTo computes xOut = x0 - x1.
func SubScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
	if !isScalarToOperable(len(mod), xOut, x0, x1) {
		panic("SubScalarTo: inputs not consistent")
	}

	for i := range mod {
		xOut[i] = num.Sub(x0[i], x1[i], mod[i])
	}
}

// NegScalar returns -x.
func NegScalar(x Scalar, mod []*num.Modulus) Scalar {
	xOut := make(Scalar, len(mod))
	NegScalarTo(xOut, x, mod)
	return xOut
}

// NegScalarTo computes xOut = -x.
func NegScalarTo(xOut, x Scalar, mod []*num.Modulus) {
	if !isScalarToOperable(len(mod), xOut, x) {
		panic("NegTo: inputs not consistent")
	}

	for i := range mod {
		xOut[i] = num.Neg(x[i], mod[i])
	}
}

// MulScalar returns x0 * x1.
func MulScalar(x0, x1 Scalar, mod []*num.Modulus) Scalar {
	xOut := make(Scalar, len(mod))
	MulScalarTo(xOut, x0, x1, mod)
	return xOut
}

// MulScalarTo computes xOut = x0 * x1.
func MulScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
	if !isScalarToOperable(len(mod), xOut, x0, x1) {
		panic("MulScalarTo: inputs not consistent")
	}

	for i := range mod {
		xOut[i] = num.Mul(x0[i], x1[i], mod[i])
	}
}

// MulAddScalarTo computes xOut += x0 * x1.
func MulAddScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
	if !isScalarToOperable(len(mod), xOut, x0, x1) {
		panic("MulAddScalarTo: inputs not consistent")
	}

	for i := range mod {
		xOut[i] = num.Add(xOut[i], num.Mul(x0[i], x1[i], mod[i]), mod[i])
	}
}

// MulSubScalarTo computes xOut -= x0 * x1.
func MulSubScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
	if !isScalarToOperable(len(mod), xOut, x0, x1) {
		panic("MulSubScalarTo: inputs not consistent")
	}

	for i := range mod {
		xOut[i] = num.Sub(xOut[i], num.Mul(x0[i], x1[i], mod[i]), mod[i])
	}
}
