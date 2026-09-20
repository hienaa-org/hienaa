package num

import (
	"cmp"
	"fmt"
	"math/bits"

	"github.com/hienaa-org/hienaa/math/internal/modops"
)

const (
	// MaxModulusBits equals to log2(MaxModulus).
	// See [MaxModulus] for details.
	MaxModulusBits = 62
	// MaxModulus is the maximum possible modulus value for the reduction.
	// All numbers in HIENAA are assumed to be less than this value.
	MaxModulus = 1 << MaxModulusBits
)

// Modulus holds precomputed constants for efficient modulus reduction.
type Modulus struct {
	// value is the raw value value.
	value uint64

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

// NewModulus creates a new [Modulus].
func NewModulus[T Integer](mod T) *Modulus {
	if mod <= 1 {
		panic("modulus must be greater than 1")
	} else if uint64(mod) > MaxModulus {
		panic("modulus must be less or equal than MaxModulus")
	}

	q := uint64(mod)

	var divHi, divLo, rem uint64
	divHi, rem = bits.Div64(1, 0, q)
	divLo, _ = bits.Div64(rem, 0, q)

	var inv uint64
	if q%2 == 1 {
		inv = 1
		acc := q
		for range 63 {
			inv *= acc
			acc *= acc
		}
	}

	return &Modulus{
		value: q,

		inv: inv,

		divHi: divHi,
		divLo: divLo,
	}
}

// Value returns the modulus value.
func (q *Modulus) Value() uint64 {
	return q.value
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

// String implements the [fmt.Stringer] interface.
func (q *Modulus) String() string {
	return fmt.Sprintf("%v", q.value)
}

// Add returns x0 + x1 mod q.
// If q is [*Modulus], x0 and x1 must be in [0, q).
// If q is nil, then it returns x0 + x1.
func Add[Q uint64 | *Modulus](x0, x1 uint64, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		return modops.AnyAdd(x0, x1, q)
	case *Modulus:
		if q == nil {
			return x0 + x1
		}
		return modops.Add(x0, x1, q.value)
	}
	return 0
}

// Sub returns x0 - x1 mod q.
// If q is [*Modulus], x0 and x1 must be in [0, q).
// If q is nil, then it returns x0 - x1.
func Sub[Q uint64 | *Modulus](x0, x1 uint64, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		return modops.AnySub(x0, x1, q)
	case *Modulus:
		if q == nil {
			return x0 - x1
		}
		return modops.Sub(x0, x1, q.value)
	}
	return 0
}

// Neg returns -x mod q.
// If q is [*Modulus], x0 and x1 must be in [0, q).
// If q is nil, then it returns -x.
func Neg[Q uint64 | *Modulus](x uint64, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		return modops.AnyNeg(x, q)
	case *Modulus:
		if q == nil {
			return -x
		}
		return modops.Neg(x, q.value)
	}
	return 0
}

// Mul returns x0 * x1 mod q.
// If q is [*Modulus], it uses Barrett reduction.
// If q is nil, then it returns x0 * x1.
func Mul[Q uint64 | *Modulus](x0, x1 uint64, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		return modops.AnyMul(x0, x1, q)
	case *Modulus:
		if q == nil {
			return x0 * x1
		}
		return modops.BMul(x0, x1, q.value, q.divHi, q.divLo)
	}
	return 0
}

// MulLazy returns x0 * x1 mod q using Barrett reduction,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func MulLazy(x0, x1 uint64, q *Modulus) uint64 {
	return modops.BMulLazy(x0, x1, q.value, q.divHi, q.divLo)
}

// Reduce returns x mod q using Barrett reduction.
// If q is [*Modulus], it uses Barrett reduction.
//
// Panics if q is nil.
func Reduce[T Integer, Q uint64 | *Modulus](x T, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		if x < 0 {
			return modops.AnyNeg(Abs(x)%q, q)
		}
		return uint64(x) % q
	case *Modulus:
		return modops.BMod(x, q.value, q.divHi)
	}
	return 0
}

// Reduce128 returns x mod q.
// If q is [*Modulus], it uses Barrett reduction.
//
// Panics if q is nil.
func Reduce128[Q uint64 | *Modulus](xHi, xLo uint64, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		return bits.Rem64(xHi, xLo, q)
	case *Modulus:
		return modops.BMod128(xHi, xLo, q.value, q.divHi, q.divLo)
	}
	return 0
}

// Reduce128Lazy returns x mod q using Barrett reduction,
// but the result is in [0, 2q).
//
// Panics if q is nil.
func Reduce128Lazy(xHi, xLo uint64, q *Modulus) uint64 {
	return modops.BMod128Lazy(xHi, xLo, q.value, q.divHi, q.divLo)
}

// Reduce2Q reduces x assuming it is in [0, 2q).
func Reduce2Q(x uint64, q *Modulus) uint64 {
	return modops.Reduce2Q(x, q.value)
}

// Reduce4Q reduces x assuming it is in [0, 4q).
func Reduce4Q(x uint64, q *Modulus) uint64 {
	return modops.Reduce4Q(x, q.value, q.value<<1)
}

// MForm transforms x into Montgomery form.
//
// Panics if q is even or nil.
func MForm(x uint64, q *Modulus) uint64 {
	if q.inv == 0 {
		panic("modulus must be odd")
	}
	return modops.MForm(x, q.value, q.divHi, q.divLo)
}

// InvMForm transforms xM to Normal form.
//
// Panics if q is even or nil.
func InvMForm(xM uint64, q *Modulus) uint64 {
	if q.inv == 0 {
		panic("modulus must be odd")
	}
	return modops.InvMForm(xM, q.value, q.inv)
}

// MMul returns x0 * x1 mod q in Montgomery form.
// x0M and x1M must be a valid Montgomery form.
//
// Panics if q is nil.
func MMul(x0M, x1M uint64, q *Modulus) uint64 {
	return modops.MMul(x0M, x1M, q.value, q.inv)
}

// MMulLazy returns x0 * x1 mod q in Montgomery form,
// but the result is in [0, 2q).
// x0M and x1M must be a valid Montgomery form.
//
// Panics if q is nil.
func MMulLazy(x0M, y0M uint64, q *Modulus) uint64 {
	return modops.MMulLazy(x0M, y0M, q.value, q.inv)
}

// SForm transforms x into Shoup form.
// x must be in [0, q).
//
// Panics if q is nil.
func SForm(x uint64, q *Modulus) uint64 {
	return modops.SForm(x, q.value)
}

// SMul returns x0 * x1 mod q using Shoup multiplication.
// x1S must be a valid Shoup form.
//
// Panics if q is nil.
func SMul(x0, x1, x1S uint64, q *Modulus) uint64 {
	return modops.SMul(x0, x1, x1S, q.value)
}

// SMulLazy returns x0 * x1 mod q using Shoup multiplication,
// but the result is in [0, 2q).
// x1S must be a valid Shoup form.
//
// Panics if q is nil.
func SMulLazy(x0, x1, x1S uint64, q *Modulus) uint64 {
	return modops.SMulLazy(x0, x1, x1S, q.value)
}

// Exp returns x^e mod q.
// If q is nil, then it returns x^e.
func Exp[Q uint64 | *Modulus](x, e uint64, q Q) uint64 {
	switch q := any(q).(type) {
	case uint64:
		if q == 1 {
			return 0
		}
	}

	r := uint64(1)
	for e > 0 {
		if e%2 == 1 {
			r = Mul(r, x, q)
		}
		e >>= 1
		x = Mul(x, x, q)
	}

	return r
}

// Inv returns the inverse of x modulo q.
//
// Panics if no inverse exists or q is nil.
func Inv[Q uint64 | *Modulus](x uint64, q Q) uint64 {
	rr := x
	var r uint64
	switch q := any(q).(type) {
	case uint64:
		r = q
	case *Modulus:
		r = q.value
	}

	ssSign, sSign := true, true
	ss, s := uint64(1), uint64(0)

	for r != 0 {
		quo := rr / r
		rr, r = r, rr-quo*r
		if sSign != ssSign {
			ss, s = s, ss+quo*s
			ssSign, sSign = sSign, ssSign
		} else {
			if ss > quo*s {
				ss, s = s, ss-quo*s
				ssSign, sSign = sSign, ssSign
			} else {
				ss, s = s, quo*s-ss
				ssSign, sSign = sSign, !sSign
			}
		}
	}

	if rr != 1 {
		panic("input not invertible")
	}

	if !ssSign {
		return Neg(ss, q)
	}
	return ss
}

// CmpModulus implements [cmp.Ordered] functionality for [Modulus].
func CmpModulus(a, b *Modulus) int {
	return cmp.Compare(a.Value(), b.Value())
}
