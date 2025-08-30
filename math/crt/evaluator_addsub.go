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

// polyScalarAddSubEvaluatorDefault is a [polyScalarAddSubEvaluator] for every ring
// except prime-order autfixed ring.
type polyScalarAddSubEvaluatorDefault struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newPolyScalarAddSubEvaluatorDefault creates a new [polyScalarEvaluatorDefault].
func newPolyScalarAddSubEvaluatorDefault(params dft.RingParameters, mod []*num.Modulus) *polyScalarAddSubEvaluatorDefault {
	isNTTFriendly := make([]bool, len(mod))
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyScalarAddSubEvaluatorDefault{
		params: params,
		mod:    mod,
	}
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarAdd(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.ScalarAddTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
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
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.ScalarSubTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorDefault) ScalarSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || isScalarToOperable(len(e.mod), c) {
		panic("ScalarAddTo: inputs not consistent")
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

// polyScalarAddSubEvaluatorAutFixedPrime is a [polyScalarAddSubEvaluator] for prime-order autfixed ring.
type polyScalarAddSubEvaluatorAutFixedPrime struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool
}

// newPolyScalarAddSubEvaluatorAutFixedPrime creates a new [polyScalarAddSubEvaluatorAutFixedPrime].
func newPolyScalarAddSubEvaluatorAutFixedPrime(params dft.RingParameters, mod []*num.Modulus) *polyScalarAddSubEvaluatorAutFixedPrime {
	isNTTFriendly := make([]bool, len(mod))
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &polyScalarAddSubEvaluatorAutFixedPrime{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,
	}
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarAdd(p *Poly, c Scalar) *Poly {
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.ScalarAddTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
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
	pOut := NewPoly(e.params.Rank(), len(e.mod))
	e.ScalarSubTo(pOut, p, c)
	return pOut
}

func (e *polyScalarAddSubEvaluatorAutFixedPrime) ScalarSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || isScalarToOperable(len(e.mod), c) {
		panic("ScalarAddTo: inputs not consistent")
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
