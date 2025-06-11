package vec

import "github.com/hienaa-org/hienaa/math/num"

// Add returns v0 + v1 mod q.
func Add(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	AddTo(v0, v1, q, vOut)
	return vOut
}

// AddLazy returns v0 + v1.
func AddLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	AddLazyTo(v0, v1, vOut)
	return vOut
}

// Sub returns v0 - v1 mod q.
func Sub(v0, v1 []uint64, q *num.Modulus) []uint64 {
	vOut := make([]uint64, len(v0))
	SubTo(v0, v1, q, vOut)
	return vOut
}

// SubLazy returns v0 - v1.
func SubLazy(v0, v1 []uint64) []uint64 {
	vOut := make([]uint64, len(v0))
	SubLazyTo(v0, v1, vOut)
	return vOut
}
