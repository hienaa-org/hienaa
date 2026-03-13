package crt

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/math/num"
)

func f128Acc(x, yHi, yLo, accHi, accLo uint64) (hi, lo uint64) {
	rLo := x * yHi
	rHi, _ := bits.Mul64(x, yLo)

	lo, carry := bits.Add64(accLo, rLo+rHi, 0)
	hi, _ = bits.Add64(accHi, 0, carry)
	return hi, lo
}

const (
	// fixedPrec is the precision of fixed-point 128-bit real number.
	fixedPrec = 124
	// floatPrec is the precision of floating-point 128-bit real number.
	floatPrec = 62
	// roundMask is the mask for rounding.
	roundMask = 1<<(fixedPrec-floatPrec) - 1
)

// mulAndFloor returns floor((x * (yHi * 2^64 + yLo)) / 2^float_prec).
// Output is in [0, 2^128).
func mulAndFloor(x, yHi, yLo uint64) (uint64, uint64) {
	rHi, rLo := bits.Mul64(x, yHi)
	rHi <<= 64 - floatPrec
	rHi += rLo >> floatPrec
	rLo <<= 64 - floatPrec

	tHi, tLo := bits.Mul64(x, yLo)
	tLo >>= floatPrec
	tLo += tHi << (64 - floatPrec)
	tHi >>= floatPrec

	rLo, carry := bits.Add64(rLo, tLo, 0)
	rHi, _ = bits.Add64(rHi, tHi, carry)

	return rHi, rLo
}

// roundTo64 returns round((xHi * 2^64 + xLo) / 2^fixed_prec).
func roundTo64(xHi, xLo uint64) uint64 {
	r := (xLo >> (fixedPrec - floatPrec)) + (xHi << (64 - fixedPrec + floatPrec))
	carry := (xLo & roundMask) >> (fixedPrec - floatPrec - 1)
	r += carry

	return r
}

// roundTo128 returns round((xHi * 2^64 + xLo) / 2^fixed_prec).
// Output is in [0, 2^128).
func roundTo128(xHi, xLo uint64) (uint64, uint64) {
	rLo := (xLo >> (fixedPrec - floatPrec)) + (xHi << (64 - fixedPrec + floatPrec))
	rHi := xHi >> (fixedPrec - floatPrec)
	carry := (xLo & roundMask) >> (fixedPrec - floatPrec - 1)

	rLo, carry = bits.Add64(rLo, carry, 0)
	rHi, _ = bits.Add64(rHi, 0, carry)

	return rHi, rLo
}

// roundTo128Signed returns round((xHi * 2^64 + xLo) / 2^fixed_prec) for int128.
func roundTo128Signed(xHi, xLo uint64) (uint64, uint64) {
	if xHi>>63 != 0 {
		rHi, rLo := roundTo128(-xHi, -xLo)
		return -rHi, -rLo
	} else {
		return roundTo128(xHi, xLo)
	}
}

// add128 returns x0 + x1.
func add128(x0Hi, x0Lo, x1Hi, x1Lo uint64) (uint64, uint64) {
	rLo, carry := bits.Add64(x0Lo, x1Lo, 0)
	rHi, _ := bits.Add64(x0Hi, x1Hi, carry)
	return rHi, rLo
}

// sub128 returns x0 - x1.
func sub128(x0Hi, x0Lo, x1Hi, x1Lo uint64) (uint64, uint64) {
	rLo, borrow := bits.Sub64(x0Lo, x1Lo, 0)
	rHi, _ := bits.Sub64(x0Hi, x1Hi, borrow)
	return rHi, rLo
}

// reduceModInToModOutSigned returns sign(x) mod qOut for x in [0, qIn).
func reduceModInToModOutSigned(x uint64, qOut *num.Modulus, qIn, halfQIn uint64) uint64 {
	if x <= halfQIn {
		return num.Reduce(x, qOut)
	}
	return num.Neg(num.Reduce(qIn-x, qOut), qOut)
}

// add64To128Signed returns x0 + int128(x1) mod q.
func add64To128Signed(x0, x1Hi, x1Lo uint64, q *num.Modulus) uint64 {
	if x1Hi>>63 != 0 {
		return num.Sub(x0, num.Reduce128(-x1Hi, -x1Lo, q), q)
	}
	return num.Add(x0, num.Reduce128(x1Hi, x1Lo, q), q)
}
