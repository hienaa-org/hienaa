package float128

import (
	"math"
	"math/big"
)

// Float128 represents a 128-bit floating point number using double-double arithmetic.
type Float128 struct {
	Hi, Lo float64
}

// FromInt64 converts a int64 to [Float128].
func FromInt64(x int64) Float128 {
	hi := float64(x)
	lo := float64(x - int64(hi))
	return Float128{hi, lo}
}

// FromRat converts a [*big.Rat] to [Float128].
func FromRat(x *big.Rat) Float128 {
	hi, _ := x.Float64()
	rHi, _ := big.NewFloat(hi).Rat(nil)
	rLo := big.NewRat(0, 1).Sub(x, rHi)
	lo, _ := rLo.Float64()
	return Float128{hi, lo}
}

// ToUint64 rounds x to uint64.
// Assumes x is positive.
func ToUint64(x Float128) uint64 {
	return uint64(int64(math.Round(x.Hi)) + int64(math.Round(x.Lo)))
}

// ToUint128 rounds x to uint128.
// Assumes x is positive.
func ToUint128(x Float128) (hi, lo uint64) {
	r := math.Pow(2, 64)
	hi = uint64(math.Round(x.Hi/r)) + uint64(math.Round(x.Lo/r))

	u64 := FromInt64(int64(hi))
	u := Sub(x, u64)
	lo = ToUint64(u)

	return
}

// Add returns x0 + x1.
func Add(x0, x1 Float128) Float128 {
	x0Hi, x0Lo := x0.Hi, x0.Lo
	x1Hi, x1Lo := x1.Hi, x1.Lo

	sHi := x0Hi + x1Hi
	sLo := x0Lo + x1Lo

	xHiF := sHi - x1Hi
	xLoF := sLo - x1Lo

	eHi := (x0Hi - xHiF) + (x1Hi - (sHi - xHiF))
	eLo := (x0Lo - xLoF) + (x1Lo - (sLo - xLoF))

	ssHi := sHi + sLo
	esLo := sLo - (ssHi - sHi)
	eHi = (eHi + eLo) + esLo

	hi := ssHi + eHi
	lo := eHi - (hi - ssHi)

	return Float128{hi, lo}
}

// Sub returns x0 - x1.
func Sub(x0, x1 Float128) Float128 {
	return Add(x0, Float128{-x1.Hi, -x1.Lo})
}

// Mul returns x0 * x1.
func Mul(x0, x1 Float128) Float128 {
	x0Hi, x0Lo := x0.Hi, x0.Lo
	x1Hi, x1Lo := x1.Hi, x1.Lo

	p00 := x0Hi * x1Hi
	e00 := math.FMA(x0Hi, x1Hi, -p00) + (x0Hi*x1Lo + x0Lo*x1Hi)

	hi := p00 + e00
	lo := e00 - (hi - p00)

	return Float128{hi, lo}
}
