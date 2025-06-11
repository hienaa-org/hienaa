package num

import (
	"encoding/binary"
	"math/big"
	"math/bits"
)

const (
	// MaxModulus is the maximum possible modulus value for the reduction.
	// Currently, this is set to 60 bits, due to various lazy reduction used in NTT/InvNTT.
	MaxModulus = 1 << 62
)

// Modulus holds precomputed constants for efficient modulus reduction.
type Modulus struct {
	// modulus is the raw modulus value.
	modulus uint64

	// inv is a constant used for Montgomery multiplication.
	// Equals to the modular inverse of modulus modulo 2^64.
	// Zero if modulus is even.
	inv uint64

	// divHi is a constant used for Barrett reduction.
	// Equals to floor(2^128 / modulus).
	divHi uint64
	// divLo is a constant used for Barrett reduction.
	// Equals to floor(2^128 / modulus).
	divLo uint64
}

// Value returns the modulus value.
func (q Modulus) Value() uint64 {
	return q.modulus
}

// Inv is a constant used for Montgomery multiplication.
// Equals to the modular inverse of modulus modulo 2^64.
// Zero if modulus is even.
func (q Modulus) Inv() uint64 {
	return q.inv
}

// Div is a constant used for Barrett reduction.
// Equals to floor(2^128 / modulus).
func (q Modulus) Div() (hi, lo uint64) {
	return q.divHi, q.divLo
}

// NewModulus creates a new Modulus.
func NewModulus(modulus uint64) *Modulus {
	switch {
	case modulus == 0:
		panic("NewModulus: modulus cannot be zero")
	case modulus > MaxModulus:
		panic("NewModulus: modulus exceeds MaxModulus")
	}

	exp128, exp64 := big.NewInt(1), big.NewInt(1)
	exp128.Lsh(exp128, 128)
	exp64.Lsh(exp64, 64)

	modulusbig := new(big.Int).SetUint64(modulus)

	modulusInvBig := big.NewInt(0)
	if modulus%2 != 0 {
		modulusInvBig.ModInverse(modulusbig, exp64)
	}

	var barConstBytes [16]byte
	barConstBig := new(big.Int).Div(exp128, modulusbig)
	barConstBig.FillBytes(barConstBytes[:])

	return &Modulus{
		modulus: modulus,

		inv: modulusInvBig.Uint64(),

		divHi: binary.BigEndian.Uint64(barConstBytes[:8]),
		divLo: binary.BigEndian.Uint64(barConstBytes[8:]),
	}
}

// Add computes x + y mod q.
func Add(x, y uint64, q *Modulus) uint64 {
	z := x + y
	if z >= q.modulus {
		z -= q.modulus
	}
	return z
}

// Sub computes x - y mod q.
func Sub(x, y uint64, q *Modulus) uint64 {
	z := x - y
	if z >= q.modulus {
		z += q.modulus
	}
	return z
}

// BMul computes x * y mod q using Barrett reduction.
func BMul(x, y uint64, q *Modulus) uint64 {
	zHi, zLo := bits.Mul64(x, y)
	return BMod(zHi, zLo, q)
}

// BMulLazy computes x * y mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulLazy(x, y uint64, q *Modulus) uint64 {
	zHi, zLo := bits.Mul64(x, y)
	return BModLazy(zHi, zLo, q)
}

// BMod computes x mod q using Barrett reduction.
func BMod(xHi, xLo uint64, q *Modulus) uint64 {
	quo := xHi * q.divHi

	quoLo, _ := bits.Mul64(xLo, q.divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, q.divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, q.divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	rem := xLo - quo*q.modulus
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// BModLazy computes x mod q using Barret reduction,
// but the result is in [0, 2q).
func BModLazy(xHi, xLo uint64, q *Modulus) uint64 {
	quo := xHi * q.divHi

	quoLo, _ := bits.Mul64(xLo, q.divLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, q.divHi)
	quo += quoMid0

	quoMid1, quoMid1Lo := bits.Mul64(xHi, q.divLo)
	quo += quoMid1

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	return xLo - quo*q.modulus
}

// MForm transforms x into Montgomery form.
func MForm(x uint64, q *Modulus) uint64 {
	xM, _ := bits.Mul64(x, q.divLo)
	xM += x * q.divHi

	rem := -xM * q.modulus
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// InvMForm transforms xM to Normal form.
func InvMForm(xM uint64, q *Modulus) uint64 {
	x, _ := bits.Mul64(xM*q.inv, q.modulus)
	rem := q.modulus - x
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// MMul computes x * y mod q in Montgomery form.
func MMul(xM, yM uint64, q *Modulus) uint64 {
	zMHi, zMLo := bits.Mul64(xM, yM)

	wHi, _ := bits.Mul64(zMLo*q.inv, q.modulus)

	rem := zMHi - wHi + q.modulus
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// MMulLazy computes x * y mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(xM, yM uint64, q *Modulus) uint64 {
	zMHi, zMLo := bits.Mul64(xM, yM)

	wHi, _ := bits.Mul64(zMLo*q.inv, q.modulus)

	return zMHi - wHi + q.modulus
}
