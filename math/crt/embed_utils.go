package crt

import (
	"github.com/hienaa-org/hienaa/math/num"
)

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
