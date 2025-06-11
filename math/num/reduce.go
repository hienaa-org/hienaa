package num

import (
	"encoding/binary"
	"math/big"

	"github.com/hienaa-org/hienaa/internal/mod"
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
	return mod.Add(x0, x1, q.modulus)
}

// Sub computes x0 - x1 mod q.
func Sub(x0, x1 uint64, q *Modulus) uint64 {
	return mod.Sub(x0, x1, q.modulus)
}

// BMul computes x * y mod q using Barrett reduction.
func BMul(x0, y0 uint64, q *Modulus) uint64 {
	return mod.BMul(x0, y0, q.modulus, q.divHi, q.divLo)
}

// BMulLazy computes x * y mod q using Barrett reduction,
// but the result is in [0, 2q).
func BMulLazy(x0, y0 uint64, q *Modulus) uint64 {
	return mod.BMulLazy(x0, y0, q.modulus, q.divHi, q.divLo)
}

// BMod computes x mod q using Barrett reduction.
func BMod(xHi, xLo uint64, q *Modulus) uint64 {
	return mod.BMod(xHi, xLo, q.modulus, q.divHi, q.divLo)
}

// BModLazy computes x mod q using Barret reduction,
// but the result is in [0, 2q).
func BModLazy(xHi, xLo uint64, q *Modulus) uint64 {
	return mod.BModLazy(xHi, xLo, q.modulus, q.divHi, q.divLo)
}

// MForm transforms x into Montgomery form.
func MForm(x uint64, q *Modulus) uint64 {
	return mod.MForm(x, q.modulus, q.divHi, q.divLo)
}

// InvMForm transforms xM to Normal form.
func InvMForm(xM uint64, q *Modulus) uint64 {
	return mod.InvMForm(xM, q.modulus, q.inv)
}

// MMul computes x0 * x1 mod q in Montgomery form.
func MMul(x0M, x1M uint64, q *Modulus) uint64 {
	return mod.MMul(x0M, x1M, q.modulus, q.inv)
}

// MMulLazy computes x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
func MMulLazy(x0M, y0M uint64, q *Modulus) uint64 {
	return mod.MMulLazy(x0M, y0M, q.modulus, q.inv)
}
