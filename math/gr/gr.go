// Package gr implements galois ring arithmetic.
package gr

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// GaloisRing represents a galois ring as [*crt.PolyEvaluator].
type GaloisRing struct {
	polyEval *crt.PolyEvaluator
	// ord is the order of the multiplicative group of the Galois ring.
	// Equals prime^((exp-1)*rank) * (prime^rank - 1).
	ord *big.Int

	buf galoisRingBuffer
}

// galoisRingBuffer is a buffer for [GaloisRing].
type galoisRingBuffer struct {
	xOut *Element
	xTmp *Element
	exp  *big.Int
}

// NewGaloisRing creates a new [GaloisRing].
//
// Panics if the correspnding conway polynomial for given modulus and rank is not found.
// In this case, find a irreducible polynomial of given rank + 1, and use
// [NewGaloisRingCustom].
func NewGaloisRing(modulus uint64, rank int) *GaloisRing {
	primes, exps := num.Factor(modulus)
	if len(primes) != 1 {
		panic("NewGaloisRing: modulus must be a prime power")
	} else if !num.IsPrime(primes[0]) {
		panic("NewGaloisRing: modulus must be a prime")
	}
	prime, _ := primes[0], exps[0]
	return NewGaloisRingCustom(modulus, findConway(prime, rank))
}

// NewGaloisRingCustom creates a new [GaloisRing] with a custom irreducible polynomial.
func NewGaloisRingCustom(modulus uint64, modPoly []int64) *GaloisRing {
	primes, exps := num.Factor(modulus)
	if len(primes) != 1 {
		panic("NewGaloisRingCustom: modulus must be a prime power")
	} else if !num.IsPrime(primes[0]) {
		panic("NewGaloisRingCustom: modulus must be a prime")
	}
	prime, exp := primes[0], exps[0]

	rank := len(modPoly) - 1
	ord := new(big.Int).SetUint64(prime)
	tmp := new(big.Int).SetUint64(prime)
	ord.Exp(ord, big.NewInt(int64(exp-1)*int64(rank)), nil)
	tmp.Exp(tmp, big.NewInt(int64(rank)), nil)
	tmp.Sub(tmp, big.NewInt(1))
	ord.Mul(ord, tmp)

	return &GaloisRing{
		polyEval: crt.NewPolyEvaluatorWithModPoly([]*num.Modulus{num.NewModulus(modulus)}, modPoly),
		ord:      ord,

		buf: newGaloisRingBuffer(modulus, rank),
	}
}

// newGaloisRingBuffer creates a new [finiteFieldBuffer].
func newGaloisRingBuffer(modulus uint64, rank int) galoisRingBuffer {
	exp := new(big.Int).SetUint64(modulus)
	exp.Exp(exp, big.NewInt(int64(rank)), nil)

	return galoisRingBuffer{
		xOut: &Element{poly: crt.NewNTTPoly(rank, 1)},
		xTmp: &Element{poly: crt.NewNTTPoly(rank, 1)},
		exp:  exp,
	}
}

// NewElement creates a new [Element] in the finite field.
func (gr *GaloisRing) NewElement() *Element {
	return &Element{
		poly: gr.polyEval.NewNTTPoly(),
	}
}

// NewElementFromUint64 creates a new [Element] from a uint64 value.
func (gr *GaloisRing) NewElementFromUint64(x uint64) *Element {
	e := gr.NewElement()
	e.poly.Coeffs[0][0] = num.Reduce(x, gr.polyEval.Modulus()[0])
	return e
}

// Modulus returns the modulus of the finite field.
func (gr *GaloisRing) Modulus() uint64 {
	return gr.polyEval.Modulus()[0].Value()
}

// Rank returns the rank of the finite field.
func (gr *GaloisRing) Rank() int {
	return gr.polyEval.Params().Rank()
}

// Ord returns the order of the multiplicative group of the Galois ring.
func (gr *GaloisRing) Ord() *big.Int {
	return gr.ord
}

// Add returns x0 + x1.
func (gr *GaloisRing) Add(x0, x1 *Element) *Element {
	xOut := gr.NewElement()
	gr.AddTo(xOut, x0, x1)
	return xOut
}

// AddTo computes xOut = x0 + x1.
func (gr *GaloisRing) AddTo(xOut, x0, x1 *Element) {
	gr.polyEval.AddTo(xOut.poly, x0.poly, x1.poly)
}

// Sub returns x0 - x1.
func (gr *GaloisRing) Sub(x0, x1 *Element) *Element {
	xOut := gr.NewElement()
	gr.SubTo(xOut, x0, x1)
	return xOut
}

// SubTo computes xOut = x0 - x1.
func (gr *GaloisRing) SubTo(xOut, x0, x1 *Element) {
	gr.polyEval.SubTo(xOut.poly, x0.poly, x1.poly)
}

// Neg returns -x.
func (gr *GaloisRing) Neg(x *Element) *Element {
	xOut := gr.NewElement()
	gr.polyEval.NegTo(xOut.poly, x.poly)
	return xOut
}

// NegTo computes xOut = -x.
func (gr *GaloisRing) NegTo(xOut, x *Element) {
	gr.polyEval.NegTo(xOut.poly, x.poly)
}

// ScalarMul returns x * c.
func (gr *GaloisRing) ScalarMul(x *Element, c uint64) *Element {
	xOut := gr.NewElement()
	gr.ScalarMulTo(xOut, x, c)
	return xOut
}

// ScalarMulTo computes xOut = x * c.
func (gr *GaloisRing) ScalarMulTo(xOut, x *Element, c uint64) {
	gr.polyEval.ScalarMulTo(xOut.poly, x.poly, crt.Scalar{c})
}

// ScalarMulAddTo computes xOut += x * c.
func (gr *GaloisRing) ScalarMulAddTo(xOut, x *Element, c uint64) {
	gr.polyEval.ScalarMulAddTo(xOut.poly, x.poly, crt.Scalar{c})
}

// ScalarMulSubTo computes xOut -= x * c.
func (gr *GaloisRing) ScalarMulSubTo(xOut, x *Element, c uint64) {
	gr.polyEval.ScalarMulSubTo(xOut.poly, x.poly, crt.Scalar{c})
}

// Mul returns x0 * x1.
func (gr *GaloisRing) Mul(x0, x1 *Element) *Element {
	xOut := gr.NewElement()
	gr.MulTo(xOut, x0, x1)
	return xOut
}

// MulTo computes xOut = x0 * x1.
func (gr *GaloisRing) MulTo(xOut, x0, x1 *Element) {
	gr.polyEval.MulTo(xOut.poly, x0.poly, x1.poly)
}

// MulAddTo computes xOut += x0 * x1.
func (gr *GaloisRing) MulAddTo(xOut, x0, x1 *Element) {
	gr.polyEval.MulAddTo(xOut.poly, x0.poly, x1.poly)
}

// MulSubTo computes xOut -= x0 * x1.
func (gr *GaloisRing) MulSubTo(xOut, x0, x1 *Element) {
	gr.polyEval.MulSubTo(xOut.poly, x0.poly, x1.poly)
}

// Exp computes xOut = x^e.
func (gr *GaloisRing) Exp(xOut, x *Element, e uint64) *Element {
	gr.ExpTo(xOut, x, e)
	return xOut
}

// ExpTo computes xOut = x^e.
func (gr *GaloisRing) ExpTo(xOut, x *Element, e uint64) {
	gr.buf.xOut.Clear()
	gr.buf.xOut.Coeffs()[0] = 1
	gr.buf.xTmp.CopyFrom(x)

	for e > 0 {
		if e&1 == 1 {
			gr.MulTo(gr.buf.xOut, gr.buf.xOut, gr.buf.xTmp)
		}
		e >>= 1
		gr.MulTo(gr.buf.xTmp, gr.buf.xTmp, gr.buf.xTmp)
	}

	xOut.CopyFrom(gr.buf.xOut)
}

// ExpBig computes xOut = x^e.
func (gr *GaloisRing) ExpBig(xOut, x *Element, e *big.Int) *Element {
	gr.ExpBigTo(xOut, x, e)
	return xOut
}

// ExpBigTo computes xOut = x^e.
func (gr *GaloisRing) ExpBigTo(xOut, x *Element, e *big.Int) {
	gr.buf.xOut.Clear()
	gr.buf.xOut.Coeffs()[0] = 1
	gr.buf.xTmp.CopyFrom(x)

	gr.buf.exp.Set(e)
	for gr.buf.exp.Sign() > 0 {
		if gr.buf.exp.Bit(0) == 1 {
			gr.MulTo(gr.buf.xOut, gr.buf.xOut, gr.buf.xTmp)
		}
		gr.buf.exp.Rsh(gr.buf.exp, 1)
		gr.MulTo(gr.buf.xTmp, gr.buf.xTmp, gr.buf.xTmp)
	}

	xOut.CopyFrom(gr.buf.xOut)
}

// Inv returns xOut = x^-1.
func (gr *GaloisRing) Inv(x *Element) *Element {
	xOut := gr.NewElement()
	gr.InvTo(xOut, x)
	return xOut
}

// InvTo computes xOut = x^-1.
func (gr *GaloisRing) InvTo(xOut, x *Element) {
	invExp := new(big.Int).Sub(gr.ord, big.NewInt(1))
	gr.ExpBigTo(xOut, x, invExp)
}

// SafeCopy returns a thread-safe copy.
func (gr *GaloisRing) SafeCopy() *GaloisRing {
	return &GaloisRing{
		polyEval: gr.polyEval.SafeCopy(),
		ord:      gr.ord,

		buf: newGaloisRingBuffer(gr.Modulus(), gr.Rank()),
	}
}
