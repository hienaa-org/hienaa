package crt

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// PolyEvaluator evaluates ring operations over [Poly].
//
// Operations usually take two forms: for example,
//   - Add(p0, p1) adds p0, p1, allocates a new vector to store the result and returns it.
//   - AddTo(pOut, p0, p1) adds p0, p1 and writes the result to pre-allocated pOut without returning.
//
// Moreover, operations panics when inputs are not consistent with the
// PolyEvaluator's parameters, or operations itself are not valid.
type PolyEvaluator struct {
	*polyEvaluatorBase
	polyScalarAddSubEvaluator
	polyMulEvaluator
	polyAutEvaluator
}

// NewPolyEvaluator creates a new [PolyEvaluator].
func NewPolyEvaluator(params dft.RingParameters, mod []*num.Modulus) *PolyEvaluator {
	switch {
	case !isCoprime(mod):
		panic("NewPolyEvaluator: modulus not coprime")
	}

	switch params.RingType() {
	case dft.Cyclotomic:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
				polyMulEvaluator:          newPolyMulEvaluatorNoReduce(params, mod),
				polyAutEvaluator:          newPolyAutEvaluatorCyclotomicPow2(params, mod),
			}
		default:
			reducer := NewCyclotomicReducer(params, mod)
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
				polyMulEvaluator:          newPolyMulEvaluatorCyclotomicNonPow2(params, mod, reducer),
				polyAutEvaluator:          newPolyAutEvaluatorCyclotomicNonPow2(params, mod, reducer),
			}
		}

	case dft.Cyclic:
		return &PolyEvaluator{
			polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
			polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
			polyMulEvaluator:          newPolyMulEvaluatorNoReduce(params, mod),
			polyAutEvaluator:          &polyAutEvaluatorPanic{},
		}

	case dft.AutFixed:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
				polyMulEvaluator:          newPolyMulEvaluatorNoReduce(params, mod),
				polyAutEvaluator:          newPolyAutEvaluatorAutFixedPow2(params, mod),
			}
		case num.IsPrime(params.CycloOrder()):
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorAutFixedPrime(params, mod),
				polyMulEvaluator:          newPolyMulEvaluatorNoReduce(params, mod),
				polyAutEvaluator:          newPolyAutEvaluatorAutFixedPrime(params, mod),
			}
		}
	}

	panic("NewPolyEvaluator: unsupported parameters")
}

// NewPolyEvaluatorWithModPoly creates a new [PolyEvaluator] with a given modulus polynomial.
func NewPolyEvaluatorWithModPoly(mod []*num.Modulus, modPoly []int64) *PolyEvaluator {
	params := dft.NewOtherParameters(modPoly)

	reducer := NewReducer(num.NextProdPower(2*len(modPoly)-1, []int{2}), mod, modPoly)

	return &PolyEvaluator{
		polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
		polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
		polyMulEvaluator:          newPolyMulEvaluatorReduce(mod, modPoly, reducer),
		polyAutEvaluator:          &polyAutEvaluatorPanic{},
	}
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *PolyEvaluator) SubEvaluator(idx ...int) *PolyEvaluator {
	return &PolyEvaluator{
		polyEvaluatorBase:         e.polyEvaluatorBase.subEvaluator(idx...),
		polyScalarAddSubEvaluator: e.polyScalarAddSubEvaluator.subEvaluator(idx...),
		polyMulEvaluator:          e.polyMulEvaluator.subEvaluator(idx...),
		polyAutEvaluator:          e.polyAutEvaluator.subEvaluator(idx...),
	}
}

// SafeCopy returns a thread-safe copy.
func (e *PolyEvaluator) SafeCopy() *PolyEvaluator {
	return &PolyEvaluator{
		polyEvaluatorBase:         e.polyEvaluatorBase.safeCopy(),
		polyScalarAddSubEvaluator: e.polyScalarAddSubEvaluator,
		polyMulEvaluator:          e.polyMulEvaluator.safeCopy(),
		polyAutEvaluator:          e.polyAutEvaluator.safeCopy(),
	}
}
