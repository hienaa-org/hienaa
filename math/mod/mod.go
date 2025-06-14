// Package mod implements modular arithmetic.
//
// Operations usually take two forms: for example,
//   - Add(v0, v1) adds v0, v1, allocates a new vector to store the result and returns it.
//   - AddTo(v0, v1, vOut) adds v0, v1 and writes the result to pre-allocated vOut without returning.
//
// Note that in most cases, v0, v1, and vOut can overlap.
// However, for operations that cannot, InPlace methods are implemented separately.
//
// For performance reasons, most functions in this package don't implement bound checks.
// If length mismatch happens, it may panic or produce wrong results.
package mod

import (
	"encoding/binary"
	"math/big"
)

const (
	// MaxModulusBits equals to log2(MaxModulus).
	// See [MaxModulus] for details.
	MaxModulusBits = 62
	// MaxModulus is the maximum possible modulus value for the reduction.
	// All numbers in HEINAA are assumed to be less than this value.
	// Currently, this is set to 60 bits, due to various lazy reduction used in NTT/InvNTT.
	MaxModulus = 1 << MaxModulusBits
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
func (q *Modulus) Value() uint64 {
	return q.modulus
}

// Inv is a constant used for Montgomery multiplication.
// Equals to the modular inverse of modulus modulo 2^64.
// Zero if modulus is even.
func (q *Modulus) Inv() uint64 {
	return q.inv
}

// Div is a constant used for Barrett reduction.
// Equals to floor(2^128 / modulus).
func (q *Modulus) Div() (hi, lo uint64) {
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

// Add computes x0 + x1 mod q.
func Add(x0, x1 uint64, q *Modulus) uint64 {
	return add(x0, x1, q.modulus)
}

// Sub computes x0 - x1 mod q.
func Sub(x0, x1 uint64, q *Modulus) uint64 {
	return sub(x0, x1, q.modulus)
}

// Mul computes x * y mod q using Barrett reduction.
func Mul(x0, y0 uint64, q *Modulus) uint64 {
	return bMul(x0, y0, q.modulus, q.divHi, q.divLo)
}

// MulLazy computes x * y mod q using Barrett reduction,
// but the result is in [0, 2q).
func MulLazy(x0, y0 uint64, q *Modulus) uint64 {
	return bMulLazy(x0, y0, q.modulus, q.divHi, q.divLo)
}

// Reduce128 computes x mod q using Barrett reduction.
func Reduce128(xHi, xLo uint64, q *Modulus) uint64 {
	return bMod128(xHi, xLo, q.modulus, q.divHi, q.divLo)
}

// Reduce computes x mod q using Barrett reduction.
func Reduce(x uint64, q *Modulus) uint64 {
	return bMod64(x, q.modulus, q.divHi)
}

// Reduce128Lazy computes x mod q using Barret reduction,
// but the result is in [0, 2q).
func Reduce128Lazy(xHi, xLo uint64, q *Modulus) uint64 {
	return bMod128Lazy(xHi, xLo, q.modulus, q.divHi, q.divLo)
}

// MForm transforms x into Montgomery form.
func MForm(x uint64, q *Modulus) uint64 {
	return mForm(x, q.modulus, q.divHi, q.divLo)
}

// InvMForm transforms xM to Normal form.
func InvMForm(xM uint64, q *Modulus) uint64 {
	return invMForm(xM, q.modulus, q.inv)
}

// MMul computes x0 * x1 mod q in Montgomery form.
func MMul(x0M, x1M uint64, q *Modulus) uint64 {
	return mMul(x0M, x1M, q.modulus, q.inv)
}

// MMulLazy computes x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(x0M, y0M uint64, q *Modulus) uint64 {
	return mMulLazy(x0M, y0M, q.modulus, q.inv)
}

// Exp computes x ** e mod q.
func Exp(x, e uint64, q *Modulus) uint64 {
	switch e {
	case 0:
		return 1
	case 1:
		return x
	}

	r := uint64(1)
	for e > 0 {
		if e&1 == 1 {
			r = Mul(r, x, q)
		}
		e >>= 1
		x = Mul(x, x, q)
	}
	return r
}
