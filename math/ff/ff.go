// Package ff implements finite field arithmetic.
package ff

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// FiniteField represents a finite field as [*crt.PolyEvaluator].
type FiniteField struct {
	polyEval *crt.PolyEvaluator
	// invExp is the exponent for the inverse operation.
	// Equals modulus^rank - 2.
	invExp *big.Int

	buf finiteFieldBuffer
}

// finiteFieldBuffer is a buffer for [FiniteField].
type finiteFieldBuffer struct {
	xOut *Element
	xTmp *Element
	exp  *big.Int
}

// NewFiniteField creates a new [FiniteField].
//
// Panics if the correspnding conway polynomial for given modulus and rank is not found.
// In this case, find a irreducible polynomial of given rank + 1, and use
// [NewFiniteFieldCustom].
func NewFiniteField(modulus uint64, rank int) *FiniteField {
	return NewFiniteFieldCustom(modulus, findConway(modulus, rank))
}

// NewFiniteFieldCustom creates a new [FiniteField] with a custom irreducible polynomial.
func NewFiniteFieldCustom(modulus uint64, modPoly []int64) *FiniteField {
	rank := len(modPoly) - 1
	invExp := new(big.Int).SetUint64(modulus)
	invExp.Exp(invExp, big.NewInt(int64(rank)), nil)
	invExp.Sub(invExp, big.NewInt(2))

	return &FiniteField{
		polyEval: crt.NewPolyEvaluatorWithModPoly([]*num.Modulus{num.NewModulus(modulus)}, modPoly),
		invExp:   invExp,

		buf: newFiniteFieldBuffer(modulus, rank),
	}
}

// newFiniteFieldBuffer creates a new [finiteFieldBuffer].
func newFiniteFieldBuffer(modulus uint64, rank int) finiteFieldBuffer {
	exp := new(big.Int).SetUint64(modulus)
	exp.Exp(exp, big.NewInt(int64(rank)), nil)

	return finiteFieldBuffer{
		xOut: &Element{poly: crt.NewNTTPoly(rank, 1)},
		xTmp: &Element{poly: crt.NewNTTPoly(rank, 1)},
		exp:  exp,
	}
}

// NewElement creates a new [Element] in the finite field.
func (f *FiniteField) NewElement() *Element {
	return &Element{
		poly: f.polyEval.NewNTTPoly(),
	}
}

// NewElementFromUint64 creates a new [Element] from a uint64 value.
func (f *FiniteField) NewElementFromUint64(x uint64) *Element {
	e := f.NewElement()
	e.poly.Coeffs[0][0] = num.Reduce(x, f.polyEval.Modulus()[0])
	return e
}

// Modulus returns the modulus of the finite field.
func (f *FiniteField) Modulus() uint64 {
	return f.polyEval.Modulus()[0].Value()
}

// Rank returns the rank of the finite field.
func (f *FiniteField) Rank() int {
	return f.polyEval.Params().Rank()
}

// Add returns x0 + x1.
func (ff *FiniteField) Add(x0, x1 *Element) *Element {
	xOut := ff.NewElement()
	ff.AddTo(xOut, x0, x1)
	return xOut
}

// AddTo computes xOut = x0 + x1.
func (ff *FiniteField) AddTo(xOut, x0, x1 *Element) {
	ff.polyEval.AddTo(xOut.poly, x0.poly, x1.poly)
}

// Sub returns x0 - x1.
func (ff *FiniteField) Sub(x0, x1 *Element) *Element {
	xOut := ff.NewElement()
	ff.SubTo(xOut, x0, x1)
	return xOut
}

// SubTo computes xOut = x0 - x1.
func (ff *FiniteField) SubTo(xOut, x0, x1 *Element) {
	ff.polyEval.SubTo(xOut.poly, x0.poly, x1.poly)
}

// Neg returns -x.
func (ff *FiniteField) Neg(x *Element) *Element {
	xOut := ff.NewElement()
	ff.polyEval.NegTo(xOut.poly, x.poly)
	return xOut
}

// NegTo computes xOut = -x.
func (ff *FiniteField) NegTo(xOut, x *Element) {
	ff.polyEval.NegTo(xOut.poly, x.poly)
}

// ScalarMul returns x * c.
func (ff *FiniteField) ScalarMul(x *Element, c uint64) *Element {
	xOut := ff.NewElement()
	ff.ScalarMulTo(xOut, x, c)
	return xOut
}

// ScalarMulTo computes xOut = x * c.
func (ff *FiniteField) ScalarMulTo(xOut, x *Element, c uint64) {
	ff.polyEval.ScalarMulTo(xOut.poly, x.poly, crt.Scalar{c})
}

// ScalarMulAddTo computes xOut += x * c.
func (ff *FiniteField) ScalarMulAddTo(xOut, x *Element, c uint64) {
	ff.polyEval.ScalarMulAddTo(xOut.poly, x.poly, crt.Scalar{c})
}

// ScalarMulSubTo computes xOut -= x * c.
func (ff *FiniteField) ScalarMulSubTo(xOut, x *Element, c uint64) {
	ff.polyEval.ScalarMulSubTo(xOut.poly, x.poly, crt.Scalar{c})
}

// Mul returns x0 * x1.
func (ff *FiniteField) Mul(x0, x1 *Element) *Element {
	xOut := ff.NewElement()
	ff.MulTo(xOut, x0, x1)
	return xOut
}

// MulTo computes xOut = x0 * x1.
func (ff *FiniteField) MulTo(xOut, x0, x1 *Element) {
	ff.polyEval.MulTo(xOut.poly, x0.poly, x1.poly)
}

// MulAddTo computes xOut += x0 * x1.
func (ff *FiniteField) MulAddTo(xOut, x0, x1 *Element) {
	ff.polyEval.MulAddTo(xOut.poly, x0.poly, x1.poly)
}

// MulSubTo computes xOut -= x0 * x1.
func (ff *FiniteField) MulSubTo(xOut, x0, x1 *Element) {
	ff.polyEval.MulSubTo(xOut.poly, x0.poly, x1.poly)
}

// Exp computes xOut = x^e.
func (ff *FiniteField) Exp(xOut, x *Element, e uint64) *Element {
	ff.ExpTo(xOut, x, e)
	return xOut
}

// ExpTo computes xOut = x^e.
func (ff *FiniteField) ExpTo(xOut, x *Element, e uint64) {
	ff.buf.xOut.Clear()
	ff.buf.xOut.Coeffs()[0] = 1
	ff.buf.xTmp.CopyFrom(x)

	for e > 0 {
		if e&1 == 1 {
			ff.MulTo(ff.buf.xOut, ff.buf.xOut, ff.buf.xTmp)
		}
		e >>= 1
		ff.MulTo(ff.buf.xTmp, ff.buf.xTmp, ff.buf.xTmp)
	}

	xOut.CopyFrom(ff.buf.xOut)
}

// ExpBig computes xOut = x^e.
func (ff *FiniteField) ExpBig(xOut, x *Element, e *big.Int) *Element {
	ff.ExpBigTo(xOut, x, e)
	return xOut
}

// ExpBigTo computes xOut = x^e.
func (ff *FiniteField) ExpBigTo(xOut, x *Element, e *big.Int) {
	ff.buf.xOut.Clear()
	ff.buf.xOut.Coeffs()[0] = 1
	ff.buf.xTmp.CopyFrom(x)

	ff.buf.exp.Set(e)
	for ff.buf.exp.Sign() > 0 {
		if ff.buf.exp.Bit(0) == 1 {
			ff.MulTo(ff.buf.xOut, ff.buf.xOut, ff.buf.xTmp)
		}
		ff.buf.exp.Rsh(ff.buf.exp, 1)
		ff.MulTo(ff.buf.xTmp, ff.buf.xTmp, ff.buf.xTmp)
	}

	xOut.CopyFrom(ff.buf.xOut)
}

// Inv returns xOut = x^-1.
func (ff *FiniteField) Inv(x *Element) *Element {
	xOut := ff.NewElement()
	ff.InvTo(xOut, x)
	return xOut
}

// InvTo computes xOut = x^-1.
func (ff *FiniteField) InvTo(xOut, x *Element) {
	ff.ExpBigTo(xOut, x, ff.invExp)
}

// SafeCopy returns a thread-safe copy.
func (ff *FiniteField) SafeCopy() *FiniteField {
	return &FiniteField{
		polyEval: ff.polyEval.SafeCopy(),
		invExp:   ff.invExp,

		buf: newFiniteFieldBuffer(ff.Modulus(), ff.Rank()),
	}
}
