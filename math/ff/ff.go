// Package ff implements finite field arithmetic.
package ff

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// FiniteField represents a finite field as [*crt.PolyEvaluator].
type FiniteField struct {
	polyEval *crt.PolyEvaluator
}

// NewFiniteField creates a new [FiniteField].
//
// Panics if the correspnding conway polynomial for given modulus and rank is not found.
// In this case, find a irreducible polynomial of given rank + 1, and use
// [NewFiniteFieldCustom].
func NewFiniteField(modulus uint64, rank int) *FiniteField {
	conway := findConway(modulus, rank)
	return &FiniteField{
		crt.NewPolyEvaluatorWithModPoly([]*num.Modulus{num.NewModulus(modulus)}, conway),
	}
}

// NewFiniteFieldCustom creates a new [FiniteField] with a custom irreducible polynomial.
func NewFiniteFieldCustom(modulus uint64, modPoly []int64) *FiniteField {
	return &FiniteField{
		crt.NewPolyEvaluatorWithModPoly([]*num.Modulus{num.NewModulus(modulus)}, modPoly),
	}
}

// NewElement creates a new [Element] in the finite field.
func (f *FiniteField) NewElement() *Element {
	return &Element{
		poly: f.polyEval.NewNTTPoly(),
	}
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
	ff.polyEval.AddTo(xOut.poly, x0.poly, x1.poly)
	return xOut
}

// AddTo computes xOut = x0 + x1.
func (ff *FiniteField) AddTo(xOut, x0, x1 *Element) {
	ff.polyEval.AddTo(xOut.poly, x0.poly, x1.poly)
}

// Sub returns x0 - x1.
func (ff *FiniteField) Sub(x0, x1 *Element) *Element {
	xOut := ff.NewElement()
	ff.polyEval.SubTo(xOut.poly, x0.poly, x1.poly)
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

// Mul returns x0 * x1.
func (ff *FiniteField) Mul(x0, x1 *Element) *Element {
	xOut := ff.NewElement()
	ff.polyEval.MulTo(xOut.poly, x0.poly, x1.poly)
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

// SafeCopy returns a thread-safe copy.
func (ff *FiniteField) SafeCopy() *FiniteField {
	return &FiniteField{
		ff.polyEval.SafeCopy(),
	}
}
