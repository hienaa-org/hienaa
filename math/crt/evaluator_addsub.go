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
	// subEvaluator returns a evaluator for modulus of given indices.
	subEvaluator(idx ...int) polyScalarAddSubEvaluator
}

// polyScalarAddSubEvaluatorDefault is a [polyScalarAddSubEvaluator] for every ring
// except prime-order autfixed ring.
type polyScalarAddSubEvaluatorDefault struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newPolyScalarAddSubEvaluatorDefault creates a new [polyScalarEvaluatorDefault].
func newPolyScalarAddSubEvaluatorDefault(params dft.RingParameters, mod []*num.Modulus) *polyScalarAddSubEvaluatorDefault {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyScalarAddSubEvaluatorDefault{
		rank: params.Rank(),
		mod:  mod,
	}
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarAdd(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarAddTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarAddTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			vec.ScalarAddTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			pOut.Coeffs[i][0] = num.Add(p.Coeffs[i][0], c[i], e.mod[i])
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarSub(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarSubTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || isScalarToOperable(len(e.mod), c) {
		panic("ScalarSubTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			vec.ScalarSubTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			pOut.Coeffs[i][0] = num.Sub(p.Coeffs[i][0], c[i], e.mod[i])
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyScalarAddSubEvaluatorDefault) subEvaluator(idx ...int) polyScalarAddSubEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return &polyScalarAddSubEvaluatorDefault{
		rank:          e.rank,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,
	}
}

// polyScalarAddSubEvaluatorAutFixedPrime is a [polyScalarAddSubEvaluator] for prime-order autfixed ring.
type polyScalarAddSubEvaluatorAutFixedPrime struct {
	rank          int
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newPolyScalarAddSubEvaluatorAutFixedPrime creates a new [polyScalarAddSubEvaluatorAutFixedPrime].
func newPolyScalarAddSubEvaluatorAutFixedPrime(params dft.RingParameters, mod []*num.Modulus) *polyScalarAddSubEvaluatorAutFixedPrime {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyScalarAddSubEvaluatorAutFixedPrime{
		rank:          params.Rank(),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarAdd(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarAddTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarAddTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			vec.ScalarAddTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			vec.ScalarSubTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarSub(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.rank, len(e.mod))
	e.ScalarSubTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.rank, len(e.mod), pOut, p) || isScalarToOperable(len(e.mod), c) {
		panic("ScalarSubTo: inputs not consistent")
	}

	for i := range e.mod {
		if p.isNTT && e.isNTTFriendly[i] {
			vec.ScalarSubTo(pOut.Coeffs[i], p.Coeffs[i], num.MForm(c[i], e.mod[i]), e.mod[i])
		} else {
			vec.ScalarAddTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
		}
	}

	pOut.isNTT = p.isNTT
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) subEvaluator(idx ...int) polyScalarAddSubEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
	}

	return &polyScalarAddSubEvaluatorAutFixedPrime{
		rank:          e.rank,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,
	}
}
