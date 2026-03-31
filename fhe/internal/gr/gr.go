// Package gr implements galois ring arithmetic.
package gr

import (
	"math/big"
	"sync"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// GaloisRing represents a galois ring as [*crt.Operator].
type GaloisRing struct {
	op *crt.Operator
	// ord is the order of the multiplicative group of the Galois ring.
	// Equals prime^((exp-1)*rank) * (prime^rank - 1).
	ord *big.Int
	// invExp is the exponent for computing inverses.
	// Equals ord - 1.
	invExp *big.Int

	pool *sync.Pool
}

// NewGaloisRing creates a new [GaloisRing].
//
// Panics if the correspnding conway polynomial for given modulus and rank is not found.
// In this case, find a irreducible polynomial of given rank + 1, and use
// [NewGaloisRingCustom].
func NewGaloisRing(modulus uint64, rank int) *GaloisRing {
	primes, _ := num.Factor(modulus)
	return NewGaloisRingCustom(modulus, findConway(primes[0], rank))
}

// NewGaloisRingCustom creates a new [GaloisRing] with a custom irreducible polynomial.
func NewGaloisRingCustom(modulus uint64, modPoly []int64) *GaloisRing {
	primes, exps := num.Factor(modulus)
	if len(primes) != 1 {
		panic("modulus must be a prime power")
	} else if !num.IsPrime(primes[0]) {
		panic("modulus must be a prime")
	}

	prime, exp := primes[0], exps[0]

	rank := len(modPoly) - 1
	ord := new(big.Int).SetUint64(prime)
	ord.Exp(ord, big.NewInt(int64(exp-1)*int64(rank)), nil)

	modSubOne := new(big.Int).SetUint64(prime)
	modSubOne.Exp(modSubOne, big.NewInt(int64(rank)), nil)
	modSubOne.Sub(modSubOne, big.NewInt(1))
	ord.Mul(ord, modSubOne)

	invExp := new(big.Int).Sub(ord, big.NewInt(1))

	return &GaloisRing{
		op:     crt.NewOperator(dft.NewOtherParameters(modPoly), []*num.Modulus{num.NewModulus(modulus)}),
		ord:    ord,
		invExp: invExp,

		pool: &sync.Pool{
			New: func() any {
				return &Element{poly: crt.NewNTTPoly(rank, 1)}
			},
		},
	}
}

// NewElement creates a new [Element] in the finite field.
func (gr *GaloisRing) NewElement() *Element {
	return &Element{
		poly: gr.op.NewNTTPoly(),
	}
}

// NewElementFromUint64 creates a new [Element] from a uint64 value.
func (gr *GaloisRing) NewElementFromUint64(x uint64) *Element {
	e := gr.NewElement()
	e.poly.Coeffs[0][0] = num.Reduce(x, gr.op.Modulus()[0])
	return e
}

// Modulus returns the modulus of the finite field.
func (gr *GaloisRing) Modulus() uint64 {
	return gr.op.Modulus()[0].Value()
}

// Rank returns the rank of the finite field.
func (gr *GaloisRing) Rank() int {
	return gr.op.Params().Rank()
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
	gr.op.AddTo(xOut.poly, x0.poly, x1.poly)
}

// Sub returns x0 - x1.
func (gr *GaloisRing) Sub(x0, x1 *Element) *Element {
	xOut := gr.NewElement()
	gr.SubTo(xOut, x0, x1)
	return xOut
}

// SubTo computes xOut = x0 - x1.
func (gr *GaloisRing) SubTo(xOut, x0, x1 *Element) {
	gr.op.SubTo(xOut.poly, x0.poly, x1.poly)
}

// Neg returns -x.
func (gr *GaloisRing) Neg(x *Element) *Element {
	xOut := gr.NewElement()
	gr.op.NegTo(xOut.poly, x.poly)
	return xOut
}

// NegTo computes xOut = -x.
func (gr *GaloisRing) NegTo(xOut, x *Element) {
	gr.op.NegTo(xOut.poly, x.poly)
}

// Mul returns x0 * x1.
func (gr *GaloisRing) Mul(x0, x1 *Element) *Element {
	xOut := gr.NewElement()
	gr.MulTo(xOut, x0, x1)
	return xOut
}

// MulTo computes xOut = x0 * x1.
func (gr *GaloisRing) MulTo(xOut, x0, x1 *Element) {
	gr.op.MulTo(xOut.poly, x0.poly, x1.poly)
}

// MulAddTo computes xOut += x0 * x1.
func (gr *GaloisRing) MulAddTo(xOut, x0, x1 *Element) {
	gr.op.MulAddTo(xOut.poly, x0.poly, x1.poly)
}

// MulSubTo computes xOut -= x0 * x1.
func (gr *GaloisRing) MulSubTo(xOut, x0, x1 *Element) {
	gr.op.MulSubTo(xOut.poly, x0.poly, x1.poly)
}

// Exp returns xOut = x^e.
func (gr *GaloisRing) Exp(x *Element, e uint64) *Element {
	xOut := gr.NewElement()
	gr.ExpTo(xOut, x, e)
	return xOut
}

// ExpTo computes xOut = x^e.
func (gr *GaloisRing) ExpTo(xOut, x *Element, e uint64) {
	xOutBuf := gr.pool.Get().(*Element)
	defer gr.pool.Put(xOutBuf)

	xBuf := gr.pool.Get().(*Element)
	defer gr.pool.Put(xBuf)

	xOutBuf.Clear()
	xOutBuf.poly.Coeffs[0][0] = 1
	xBuf.CopyFrom(x)

	for e > 0 {
		if e&1 == 1 {
			gr.MulTo(xOutBuf, xOutBuf, xBuf)
		}
		e >>= 1
		gr.MulTo(xBuf, xBuf, xBuf)
	}

	xOut.CopyFrom(xOutBuf)
}

// ExpBig returns xOut = x^e.
func (gr *GaloisRing) ExpBig(xOut, x *Element, e *big.Int) *Element {
	gr.ExpBigTo(xOut, x, e)
	return xOut
}

// ExpBigTo computes xOut = x^e.
func (gr *GaloisRing) ExpBigTo(xOut, x *Element, e *big.Int) {
	exp := new(big.Int).Set(e)

	xBuf := gr.pool.Get().(*Element)
	defer gr.pool.Put(xBuf)

	xOutBuf := gr.pool.Get().(*Element)
	defer gr.pool.Put(xOutBuf)

	xOutBuf.Clear()
	xOutBuf.poly.Coeffs[0][0] = 1
	xBuf.CopyFrom(x)

	exp.Set(e)
	for exp.Sign() > 0 {
		if exp.Bit(0) == 1 {
			gr.MulTo(xOutBuf, xOutBuf, xBuf)
		}
		exp.Rsh(exp, 1)
		gr.MulTo(xBuf, xBuf, xBuf)
	}

	xOut.CopyFrom(xOutBuf)
}

// Inv returns xOut = x^-1.
func (gr *GaloisRing) Inv(x *Element) *Element {
	xOut := gr.NewElement()
	gr.InvTo(xOut, x)
	return xOut
}

// InvTo computes xOut = x^-1.
func (gr *GaloisRing) InvTo(xOut, x *Element) {
	gr.ExpBigTo(xOut, x, gr.invExp)
}
