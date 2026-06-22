package num

import (
	"cmp"
	"fmt"
	"math"
	"math/bits"

	"github.com/hienaa-org/hienaa/math/internal/modops"
)

const (
	// MaxModulusBits equals to log2(MaxModulus).
	// See [MaxModulus] for details.
	MaxModulusBits = 50
	// MaxModulus is the maximum possible modulus value for the reduction.
	// All numbers in HIENAA are assumed to be less than this value.
	MaxModulus = 1 << MaxModulusBits
	// MinAmbModulusBits is the minimum possible bits for ambient modulus limbs.
	// Currently this is set to 48, because we use [vec.Reduce4Q] for embedding to ambient modulus.
	MinAmbModulusBits = 48
)

// Modulus holds precomputed constants for efficient modulus reduction.
type Modulus struct {
	// modulus is the raw modulus value.
	modulus uint64

	// float equals float64(modulus).
	float float64
	// floatInv equals 1 / float64(modulus).
	floatInv float64

	// div is a constant used for Barrett reduction.
	div uint64
	// log equals floor(log2(modulus)).
	log uint64
}

// NewModulus creates a new [Modulus].
func NewModulus[T Integer](mod T) *Modulus {
	if mod <= 1 {
		panic("modulus must be greater than 1")
	} else if uint64(mod) >= MaxModulus {
		panic("modulus must be less than MaxModulus")
	}

	q := uint64(mod)

	log := uint64(bits.Len64(q) - 1)
	if IsPowerOfTwo(q) {
		log -= 1
	}
	div, _ := bits.Div64(1<<log, 0, q)

	return &Modulus{
		modulus: q,

		float:    float64(q),
		floatInv: math.Nextafter(1/float64(q), math.Inf(1)),

		div: div,
		log: log,
	}
}

// Value returns the modulus value.
func (q *Modulus) Value() uint64 {
	return q.modulus
}

// Div is a constant used for Barrett reduction.
func (q *Modulus) Div() uint64 {
	return q.div
}

// Log equals floor(log2(modulus)).
func (q *Modulus) Log() uint64 {
	return q.log
}

// Float equals float64(q).
func (q *Modulus) Float() float64 {
	return q.float
}

// FloatInv equals 1 / float64(q).
func (q *Modulus) FloatInv() float64 {
	return q.floatInv
}

// String implements the [fmt.Stringer] interface.
func (q *Modulus) String() string {
	return fmt.Sprintf("%v", q.modulus)
}

// MulForm stores precomputed values for modular multiplication.
type MulForm struct {
	// Float is the float64 representation of the value.
	Float float64
	// SForm is the Shoup Form of the value.
	SForm uint64
}

// Add returns x0 + x1 mod q.
// x0 and x1 must be in [0, q).
// If q is nil, then it returns x0 + x1.
func Add(x0, x1 uint64, q *Modulus) uint64 {
	if q != nil {
		return modops.Add(x0, x1, q.modulus)
	}
	return x0 + x1
}

// Sub returns x0 - x1 mod q.
// x0 and x1 must be in [0, q).
// If q is nil, then it returns x0 - x1.
func Sub(x0, x1 uint64, q *Modulus) uint64 {
	if q != nil {
		return modops.Sub(x0, x1, q.modulus)
	}
	return x0 - x1
}

// Neg returns -x mod q.
// x must be in [0, q).
// If q is nil, then it returns -x.
func Neg(x uint64, q *Modulus) uint64 {
	if q != nil {
		return modops.Neg(x, q.modulus)
	}
	return -x
}

// Mul returns x0 * x1 mod q using Barrett reduction.
// x0 and x1 must be in [0, q).
// If q is nil, then it returns x0 * x1.
func Mul(x0, x1 uint64, q *Modulus) uint64 {
	if q != nil {
		return modops.Mul(x0, x1, q.modulus, q.div, q.log)
	}
	return x0 * x1
}

// Reduce returns x mod q using Barrett reduction.
//
// Panics if q is nil.
func Reduce[T Integer](x T, q *Modulus) uint64 {
	return modops.BMod(x, q.modulus, q.div, q.log)
}

// MForm transforms x into [MulForm].
//
// Panics if q is nil.
func ToMulForm(x uint64, q *Modulus) MulForm {
	return MulForm{
		Float: float64(x),
		SForm: modops.SForm(x, q.modulus),
	}
}

// FMul returns x0 * x1 mod q using [MulForm] of x1.
//
// Panics if q is nil.
func FMul(x0, x1 uint64, x1M MulForm, q *Modulus) uint64 {
	return modops.SMul(x0, x1, x1M.SForm, q.modulus)
}

// Exp returns x^e mod q.
// If q is nil, then it returns x^e.
func Exp(x, e uint64, q *Modulus) uint64 {
	switch e {
	case 0:
		return 1
	case 1:
		return x
	}

	r := uint64(1)
	if q == nil {
		for e > 0 {
			if e%2 == 1 {
				r = r * x
			}
			e >>= 1
			x = x * x
		}
	} else {
		for e > 0 {
			if e%2 == 1 {
				r = Mul(r, x, q)
			}
			e >>= 1
			x = Mul(x, x, q)
		}
	}

	return r
}

// Inv returns the inverse of x modulo q.
//
// Panics if q is nil.
func Inv(x uint64, q *Modulus) uint64 {
	rr, r := x, q.Value()

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
		return q.Value() - ss
	}
	return ss
}

// CmpModulus implements [cmp.Ordered] functionality for [Modulus].
func CmpModulus(a, b *Modulus) int {
	return cmp.Compare(a.Value(), b.Value())
}
