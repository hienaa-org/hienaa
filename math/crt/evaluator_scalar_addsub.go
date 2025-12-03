package crt

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// polyScalarAddSubEvaluator is the evaluator for addition/subtraction.
// All polyScalarAddSubEvaluators are assumed to be thread-safe.
type polyScalarAddSubEvaluator interface {
	// ScalarAdd returns p + c.
	ScalarAdd(p *Poly, c Scalar) *Poly
	// ScalarAddTo computes pOut = p + c.
	ScalarAddTo(pOut, p *Poly, c Scalar)
	// ScalarSub returns p - c.
	ScalarSub(p *Poly, c Scalar) *Poly
	// ScalarSubTo computes pOut = p - c.
	ScalarSubTo(pOut, p *Poly, c Scalar)
}

// defaultPolyScalarAddSubEvaluator is a [polyScalarAddSubEvaluator] for every ring
// except prime-order autfixed ring.
type defaultPolyScalarAddSubEvaluator struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newDefaultPolyScalarAddSubEvaluator creates a new [defaultPolyScalarAddSubEvaluator].
func newDefaultPolyScalarAddSubEvaluator(params dft.RingParameters, mod []*num.Modulus) defaultPolyScalarAddSubEvaluator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return defaultPolyScalarAddSubEvaluator{
		rank:          params.Rank(),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

// ScalarAdd returns p + c.
func (e *defaultPolyScalarAddSubEvaluator) ScalarAdd(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarAddTo(pOut, p, c)
	return pOut
}

// ScalarAddTo computes pOut = p + c.
func (e *defaultPolyScalarAddSubEvaluator) ScalarAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarAddTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
			vec.ScalarAddTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
			pOut.Coeffs[i][0] = num.Add(pOut.Coeffs[i][0], c[i], e.mod[i])
		}
	}

	pOut.IsNTT = p.IsNTT
}

// ScalarSub returns p - c.
func (e *defaultPolyScalarAddSubEvaluator) ScalarSub(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarSubTo(pOut, p, c)
	return pOut
}

// ScalarSubTo computes pOut = p - c.
func (e *defaultPolyScalarAddSubEvaluator) ScalarSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarSubTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
			vec.ScalarSubTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
			pOut.Coeffs[i][0] = num.Sub(pOut.Coeffs[i][0], c[i], e.mod[i])
		}
	}

	pOut.IsNTT = p.IsNTT
}

func (e *defaultPolyScalarAddSubEvaluator) subEvaluator(idx ...int) defaultPolyScalarAddSubEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return defaultPolyScalarAddSubEvaluator{
		rank:          e.rank,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,
	}
}

// primeAutFixedPolyScalarAddSubEvaluator is a [polyScalarAddSubEvaluator] for prime-order autfixed ring.
type primeAutFixedPolyScalarAddSubEvaluator struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newPrimeAutFixedPolyScalarAddSubEvaluator creates a new [primeAutFixedPolyScalarAddSubEvaluator].
func newPrimeAutFixedPolyScalarAddSubEvaluator(params dft.RingParameters, mod []*num.Modulus) primeAutFixedPolyScalarAddSubEvaluator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return primeAutFixedPolyScalarAddSubEvaluator{
		rank:          params.Rank(),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

// ScalarAdd returns p + c.
func (e *primeAutFixedPolyScalarAddSubEvaluator) ScalarAdd(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarAddTo(pOut, p, c)
	return pOut
}

// ScalarAddTo computes pOut = p + c.
func (e *primeAutFixedPolyScalarAddSubEvaluator) ScalarAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarAddTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
			vec.ScalarAddTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
			pOut.Coeffs[i][0] = num.Add(pOut.Coeffs[i][0], c[i], e.mod[i])
		}
	}

	pOut.IsNTT = p.IsNTT
}

// ScalarSub returns p - c.
func (e *primeAutFixedPolyScalarAddSubEvaluator) ScalarSub(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarSubTo(pOut, p, c)
	return pOut
}

// ScalarSubTo computes pOut = p - c.
func (e *primeAutFixedPolyScalarAddSubEvaluator) ScalarSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarSubTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.IsNTT && e.isNTTFriendly[i] {
			vec.ScalarSubTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
			pOut.Coeffs[i][0] = num.Sub(pOut.Coeffs[i][0], c[i], e.mod[i])
		}
	}

	pOut.IsNTT = p.IsNTT
}

func (e *primeAutFixedPolyScalarAddSubEvaluator) subEvaluator(idx ...int) primeAutFixedPolyScalarAddSubEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return primeAutFixedPolyScalarAddSubEvaluator{
		rank:          e.rank,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,
	}
}
