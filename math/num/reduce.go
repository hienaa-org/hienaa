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

	// modulusInv is the modular inverse of modulus modulo 2^64.
	// Zero if modulus is even.
	modulusInv uint64

	// barConst is a constant used for Barrett reduction.
	// Equals to floor(2^128 / modulus).
	barConstHi uint64
	// barConst is a constant used for Barrett reduction.
	// Equals to floor(2^128 / modulus).
	barConstLo uint64
}

// Value returns the modulus value.
func (q Modulus) Value() uint64 {
	return q.modulus
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

		modulusInv: modulusInvBig.Uint64(),

		barConstHi: binary.BigEndian.Uint64(barConstBytes[:8]),
		barConstLo: binary.BigEndian.Uint64(barConstBytes[8:]),
	}
}

// Add computes x + y mod q.
func Add(x, y uint64, q *Modulus) uint64 {
	rem := x + y
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// Sub computes x - y mod q.
func Sub(x, y uint64, q *Modulus) uint64 {
	rem := x - y
	if rem >= q.modulus {
		rem += q.modulus
	}
	return rem
}

// BMod computes x mod q using Barrett reduction.
func BMod(xHi, xLo uint64, q *Modulus) uint64 {
	quo := xHi * q.barConstHi

	quoLo, _ := bits.Mul64(xLo, q.barConstLo)

	quoMid0, quoMid0Lo := bits.Mul64(xLo, q.barConstHi)
	quo += quoMid0

	quoMidSum, quoMidCarry := bits.Add64(quoMid0Lo, quoLo, 0)
	quo += quoMidCarry

	quoMid1, quoMid1Lo := bits.Mul64(xHi, q.barConstLo)
	quo += quoMid1

	_, quoMidCarry = bits.Add64(quoMid1Lo, quoMidSum, 0)
	quo += quoMidCarry

	rem := xLo - quo*q.modulus
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// BMul computes x * y mod q using Barrett reduction.
func BMul(x, y uint64, q *Modulus) uint64 {
	zHi, zLo := bits.Mul64(x, y)
	return BMod(zHi, zLo, q)
}

// MForm transforms x into Montgomery form.
func MForm(x uint64, q *Modulus) uint64 {
	xM, _ := bits.Mul64(x, q.barConstLo)
	xM += x * q.barConstHi

	rem := -xM * q.modulus
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// InvMForm transforms xM to Normal form.
func InvMForm(xM uint64, q *Modulus) uint64 {
	x, _ := bits.Mul64(xM*q.modulusInv, q.modulus)
	rem := q.modulus - x
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}

// MMul computes x * y mod q in Montgomery form.
func MMul(xM, yM uint64, q *Modulus) uint64 {
	zMHi, zMLo := bits.Mul64(xM, yM)

	wHi, _ := bits.Mul64(zMLo*q.modulusInv, q.modulus)

	rem := zMHi - wHi + q.modulus
	if rem >= q.modulus {
		rem -= q.modulus
	}
	return rem
}
