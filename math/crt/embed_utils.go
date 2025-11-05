package crt

import (
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/num"
)

func qFloatFromBigRat(f *big.Rat) (hi, lo float64) {
	hi, _ = f.Float64()
	rHi, _ := big.NewFloat(hi).Rat(nil)
	rLo := big.NewRat(0, 1).Sub(f, rHi)
	lo, _ = rLo.Float64()

	return
}

func qFloatFromInt[T num.Integer](x T) (hi, lo float64) {
	hi = float64(x)
	x -= T(hi)
	lo = float64(x)

	return
}

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

func qFloatMul(x0Hi, x0Lo, x1Hi, x1Lo float64) (hi, lo float64) {
	p00 := x0Hi * x1Hi
	e00 := math.FMA(x0Hi, x1Hi, -p00) + (x0Hi*x1Lo + x0Lo*x1Hi)
	hi = p00 + e00
	lo = e00 - (hi - p00)

	return
}

func qFloatRoundAsInt[T num.Integer](xHi, xLo float64) T {
	return T(math.Round(xHi)) + T(math.Round(xLo))
}

func qFloatRoundAsUint128(xHi, xLo float64) (hi, lo uint64) {
	hi = uint64(math.Round(xHi/math.Pow(2, 64))) + uint64(math.Round(xLo/math.Pow(2, 64)))

	x64Hi, x64Lo := qFloatFromInt(hi)
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
