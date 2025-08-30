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
	gcd := uint64(1)
	for i := range mod {
		gcd = num.GCD(gcd, mod[i].Value())
	}
	if gcd != 1 {
		panic("NewPolyEvaluator: moduli not coprime")
	}

	switch params.RingType() {
	case dft.Cyclotomic:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
				polyMulEvaluator:          NewPolyMulEvaluatorNoReduce(params, mod),
				polyAutEvaluator:          newPolyAutEvaluatorCyclotomicPow2(params, mod),
			}
		default:
			reducers := make([]reducer, len(mod))
			for i := range reducers {
				reducers[i] = newCyclotomicReducer(params, mod[i])
			}
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
				polyMulEvaluator:          newPolyMulEvaluatorCyclotomicNonPow2(params, mod, reducers),
				polyAutEvaluator:          newPolyAutEvaluatorCyclotomicNonPow2(params, mod, reducers),
			}
		}

	case dft.Cyclic:
		return &PolyEvaluator{
			polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
			polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
			polyMulEvaluator:          NewPolyMulEvaluatorNoReduce(params, mod),
			polyAutEvaluator:          &polyAutEvaluatorPanic{},
		}

	case dft.AutFixed:
		switch {
		case num.IsPrime(params.CycloOrder()):
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorAutFixedPrime(params, mod),
				polyMulEvaluator:          NewPolyMulEvaluatorNoReduce(params, mod),
				polyAutEvaluator:          newPolyAutEvaluatorAutFixedPrime(params, mod),
			}
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &PolyEvaluator{
				polyEvaluatorBase:         newPolyEvaluatorBase(params, mod),
				polyScalarAddSubEvaluator: newPolyScalarAddSubEvaluatorDefault(params, mod),
				polyMulEvaluator:          NewPolyMulEvaluatorNoReduce(params, mod),
				polyAutEvaluator:          newPolyAutEvaluatorAutFixedPow2(params, mod),
			}
		}
	}

	panic("NewPolyEvaluator: unsupported parameters")
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
