package crt

import (
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

// qFloatFromBigRat converts a big.Rat to a quadruple float.
func qFloatFromBigRat(f *big.Rat) (hi, lo float64) {
	hi, _ = f.Float64()
	rHi, _ := big.NewFloat(hi).Rat(nil)
	rLo := big.NewRat(0, 1).Sub(f, rHi)
	lo, _ = rLo.Float64()

	return
}

// qFloatFromUint64 converts a uint64 to a quadruple float.
func qFloatFromUint64(x uint64) (hi, lo float64) {
	xint := int64(x)

	hi = float64(xint)
	xint -= int64(hi)
	lo = float64(xint)

	return
}

// qFloatAdd adds two quadruple floats.
func qFloatAdd(x0Hi, x0Lo, x1Hi, x1Lo float64) (hi, lo float64) {
	sHi := x0Hi + x1Hi
	sLo := x0Lo + x1Lo

	xHiF := sHi - x1Hi
	xLoF := sLo - x1Lo

	eHi := (x0Hi - xHiF) + (x1Hi - (sHi - xHiF))
	eLo := (x0Lo - xLoF) + (x1Lo - (sLo - xLoF))

	ssHi := sHi + sLo
	esLo := sLo - (ssHi - sHi)
	eHi = (eHi + eLo) + esLo

	hi = ssHi + eHi
	lo = eHi - (hi - ssHi)

	return
}

// qFloatSub subtracts two quadruple floats.
func qFloatSub(x0Hi, x0Lo, x1Hi, x1Lo float64) (hi, lo float64) {
	sHi := x0Hi - x1Hi
	sLo := x0Lo - x1Lo

	xHiF := sHi + x1Hi
	xLoF := sLo + x1Lo

	eHi := (x0Hi - xHiF) + (-x1Hi - (sHi - xHiF))
	eLo := (x0Lo - xLoF) + (-x1Lo - (sLo - xLoF))

	ssHi := sHi + sLo
	esLo := sLo - (ssHi - sHi)
	eHi = (eHi + eLo) + esLo

	hi = ssHi + eHi
	lo = eHi - (hi - ssHi)

	return
}

// qFloatMul multiplies two quadruple floats.
func qFloatMul(x0Hi, x0Lo, x1Hi, x1Lo float64) (hi, lo float64) {
	p00 := x0Hi * x1Hi
	e00 := math.FMA(x0Hi, x1Hi, -p00) + (x0Hi*x1Lo + x0Lo*x1Hi)
	hi = p00 + e00
	lo = e00 - (hi - p00)

	return
}

// qFloatRoundAsUint64 rounds a quadruple float to a uint64.
// This function assumes that the input quadruple float is positive.
func qFloatRoundAsUint64(xHi, xLo float64) uint64 {
	return uint64(int64(math.Round(xHi)) + int64(math.Round(xLo)))
}

// qFloatRoundAsUint128 rounds a quadruple float to a uint128.
// This function assumes that the input quadruple float is positive.
func qFloatRoundAsUint128(xHi, xLo float64) (hi, lo uint64) {
	hi = uint64(math.Round(xHi/math.Pow(2, 64))) + uint64(math.Round(xLo/math.Pow(2, 64)))

	x64Hi, x64Lo := qFloatFromUint64(hi)
	xHi, xLo = qFloatSub(xHi, xLo, x64Hi, x64Lo)
	lo = uint64(math.Round(xHi)) + uint64(math.Round(xLo))

	return
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
