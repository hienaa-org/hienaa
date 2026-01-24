package crt

import (
	"math/big"

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
type PolyEvaluator interface {
	// NewPoly creates a new [Poly].
	NewPoly() *Poly
	// NewNTTPoly creates a new [Poly] in NTT form.
	NewNTTPoly() *Poly
	// NewPolyCustom creates a new [Poly] with the given parameters.
	NewPolyCustom(isNTT bool) *Poly
	// Params returns the ring parameters.
	Params() dft.RingParameters
	// Modulus returns the modulus.
	Modulus() []*num.Modulus
	// ModulusPoly returns the quotient polynomial of the ring.
	ModulusPoly() []int64

	// FwdNTT returns FwdNTT(p).
	FwdNTT(p *Poly) *Poly
	// FwdNTTTo computes pOut = NTT(p).
	FwdNTTTo(pNTT, p *Poly)
	// InvNTT returns InvNTT(p).
	InvNTT(pNTT *Poly) *Poly
	// InvNTTTo computes pOut = InvNTT(p).
	InvNTTTo(p, pNTT *Poly)

	// Add returns p0 + p1.
	Add(p0, p1 *Poly) *Poly
	// AddTo computes pOut = p0 + p1.
	AddTo(pOut, p0, p1 *Poly)

	// Sub returns p0 - p1.
	Sub(p0, p1 *Poly) *Poly
	// SubTo computes pOut = p0 - p1.
	SubTo(pOut, p0, p1 *Poly)

	// Neg returns -p.
	Neg(p *Poly) *Poly
	// NegTo computes pOut = -p.
	NegTo(pOut, p *Poly)

	// ScalarMul returns p * c.
	ScalarMul(p *Poly, c Scalar) *Poly
	// ScalarMulTo computes pOut = p * c.
	ScalarMulTo(pOut, p *Poly, c Scalar)
	// ScalarMulAddTo computes pOut += p * c.
	ScalarMulAddTo(pOut, p *Poly, c Scalar)
	// ScalarMulSubTo computes pOut -= p * c.
	ScalarMulSubTo(pOut, p *Poly, c Scalar)

	polyScalarAddSubEvaluator

	polyMulEvaluator

	polyAutEvaluator

	// AsBig returns p as *[big.Int] vector.
	AsBig(p *Poly) []*big.Int

	// SubEvaluator returns a evaluator for modulus of given indices.
	// Useful for "levelled" operations, especially with [vec.Range].
	SubEvaluator(idx ...int) PolyEvaluator
	// SafeCopy returns a thread-safe copy.
	SafeCopy() PolyEvaluator
}

// NewPolyEvaluator creates a new [PolyEvaluator].
func NewPolyEvaluator(params dft.RingParameters, mod []*num.Modulus) PolyEvaluator {
	if !isCoprime(mod) {
		panic("modulus must be coprime")
	}

	switch params.RingType() {
	case dft.Cyclotomic:
		modPoly := dft.CyclotomicPolynomial(params.CycloOrder())
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &pow2CyclotomicPolyEvaluator{
				polyBaseEvaluator:                newPolyBaseEvaluator(params, mod, modPoly),
				defaultPolyScalarAddSubEvaluator: newDefaultPolyScalarAddSubEvaluator(params, mod),
				noReducePolyMulEvaluator:         newNoReducePolyMulEvaluator(params, mod),
				pow2CyclotomicPolyAutEvaluator:   newPow2CyclotomicPolyAutEvaluator(params, mod),
			}
		default:
			reducer := NewCyclotomicReducer(params, mod)
			return &anyCyclotomicPolyEvaluator{
				polyBaseEvaluator:                newPolyBaseEvaluator(params, mod, modPoly),
				defaultPolyScalarAddSubEvaluator: newDefaultPolyScalarAddSubEvaluator(params, mod),
				anyCyclotomicPolyMulEvaluator:    newAnyCyclotomicPolyMulEvaluator(params, mod, reducer),
				anyCyclotomicPolyAutEvaluator:    newAnyCyclotomicPolyAutEvaluator(params, mod, reducer),
			}
		}

	case dft.Cyclic:
		modPoly := make([]int64, params.Rank()+1)
		return &pow235CyclicPolyEvaluator{
			polyBaseEvaluator:                newPolyBaseEvaluator(params, mod, modPoly),
			defaultPolyScalarAddSubEvaluator: newDefaultPolyScalarAddSubEvaluator(params, mod),
			noReducePolyMulEvaluator:         newNoReducePolyMulEvaluator(params, mod),
			noAutPolyAutEvaluator:            noAutPolyAutEvaluator{},
		}

	case dft.AutFixed:
		modPoly := dft.CyclotomicPolynomial(params.CycloOrder())
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &pow2AutFixedPolyEvaluator{
				polyBaseEvaluator:                newPolyBaseEvaluator(params, mod, modPoly),
				defaultPolyScalarAddSubEvaluator: newDefaultPolyScalarAddSubEvaluator(params, mod),
				noReducePolyMulEvaluator:         newNoReducePolyMulEvaluator(params, mod),
				pow2AutFixedPolyAutEvaluator:     newPow2AutFixedPolyAutEvaluator(params, mod),
			}
		case num.IsPrime(params.CycloOrder()):
			return &primeAutFixedPolyEvaluator{
				polyBaseEvaluator:                      newPolyBaseEvaluator(params, mod, modPoly),
				primeAutFixedPolyScalarAddSubEvaluator: newPrimeAutFixedPolyScalarAddSubEvaluator(params, mod),
				noReducePolyMulEvaluator:               newNoReducePolyMulEvaluator(params, mod),
				primeAutFixedPolyAutEvaluator:          newPrimeAutFixedPolyAutEvaluator(params, mod),
			}
		}
	}

	panic("unsupported parameters")
}

// NewPolyEvaluatorWithModPoly creates a new [PolyEvaluator] with a given modulus polynomial.
func NewPolyEvaluatorWithModPoly(mod []*num.Modulus, modPoly []int64) PolyEvaluator {
	params := dft.NewOtherParameters(modPoly)

	reducer := NewReducer(num.NextProdPower(2*len(modPoly)-1, []int{2}), mod, modPoly)

	return &anyPolyEvaluator{
		polyBaseEvaluator:                newPolyBaseEvaluator(params, mod, modPoly),
		defaultPolyScalarAddSubEvaluator: newDefaultPolyScalarAddSubEvaluator(params, mod),
		reducePolyMulEvaluator:           newReducePolyMulEvaluator(mod, modPoly, reducer),
		noAutPolyAutEvaluator:            noAutPolyAutEvaluator{},
	}
}

// pow2CyclotomicPolyEvaluator is a PolyEvaluator for power-of-two cyclotomic ring.
type pow2CyclotomicPolyEvaluator struct {
	polyBaseEvaluator
	defaultPolyScalarAddSubEvaluator
	noReducePolyMulEvaluator
	pow2CyclotomicPolyAutEvaluator
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *pow2CyclotomicPolyEvaluator) SubEvaluator(idx ...int) PolyEvaluator {
	return &pow2CyclotomicPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.subEvaluator(idx...),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator.subEvaluator(idx...),
		noReducePolyMulEvaluator:         e.noReducePolyMulEvaluator.subEvaluator(idx...),
		pow2CyclotomicPolyAutEvaluator:   e.pow2CyclotomicPolyAutEvaluator.subEvaluator(idx...),
	}
}

// SafeCopy returns a thread-safe copy.
func (e *pow2CyclotomicPolyEvaluator) SafeCopy() PolyEvaluator {
	return &pow2CyclotomicPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.safeCopy(),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator,
		noReducePolyMulEvaluator:         e.noReducePolyMulEvaluator.safeCopy(),
		pow2CyclotomicPolyAutEvaluator:   e.pow2CyclotomicPolyAutEvaluator.safeCopy(),
	}
}

// anyCyclotomicPolyEvaluator is a PolyEvaluator for arbitrary order cyclotomic ring.
type anyCyclotomicPolyEvaluator struct {
	polyBaseEvaluator
	defaultPolyScalarAddSubEvaluator
	anyCyclotomicPolyMulEvaluator
	anyCyclotomicPolyAutEvaluator
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *anyCyclotomicPolyEvaluator) SubEvaluator(idx ...int) PolyEvaluator {
	return &anyCyclotomicPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.subEvaluator(idx...),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator.subEvaluator(idx...),
		anyCyclotomicPolyMulEvaluator:    e.anyCyclotomicPolyMulEvaluator.subEvaluator(idx...),
		anyCyclotomicPolyAutEvaluator:    e.anyCyclotomicPolyAutEvaluator.subEvaluator(idx...),
	}
}

// SafeCopy returns a thread-safe copy.
func (e *anyCyclotomicPolyEvaluator) SafeCopy() PolyEvaluator {
	return &anyCyclotomicPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.safeCopy(),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator,
		anyCyclotomicPolyMulEvaluator:    e.anyCyclotomicPolyMulEvaluator.safeCopy(),
		anyCyclotomicPolyAutEvaluator:    e.anyCyclotomicPolyAutEvaluator.safeCopy(),
	}
}

// pow235CyclicPolyEvaluator is a PolyEvaluator for cyclic ring with ranks multiple of 2, 3 and 5.
type pow235CyclicPolyEvaluator struct {
	polyBaseEvaluator
	defaultPolyScalarAddSubEvaluator
	noReducePolyMulEvaluator
	noAutPolyAutEvaluator
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *pow235CyclicPolyEvaluator) SubEvaluator(idx ...int) PolyEvaluator {
	return &pow235CyclicPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.subEvaluator(idx...),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator.subEvaluator(idx...),
		noReducePolyMulEvaluator:         e.noReducePolyMulEvaluator.subEvaluator(idx...),
		noAutPolyAutEvaluator:            e.noAutPolyAutEvaluator,
	}
}

// SafeCopy returns a thread-safe copy.
func (e *pow235CyclicPolyEvaluator) SafeCopy() PolyEvaluator {
	return &pow235CyclicPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.safeCopy(),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator,
		noReducePolyMulEvaluator:         e.noReducePolyMulEvaluator.safeCopy(),
		noAutPolyAutEvaluator:            e.noAutPolyAutEvaluator,
	}
}

// pow2AutFixedPolyEvaluator is a PolyEvaluator for power-of-two conjugate invariant ring.
type pow2AutFixedPolyEvaluator struct {
	polyBaseEvaluator
	defaultPolyScalarAddSubEvaluator
	noReducePolyMulEvaluator
	pow2AutFixedPolyAutEvaluator
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *pow2AutFixedPolyEvaluator) SubEvaluator(idx ...int) PolyEvaluator {
	return &pow2AutFixedPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.subEvaluator(idx...),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator.subEvaluator(idx...),
		noReducePolyMulEvaluator:         e.noReducePolyMulEvaluator.subEvaluator(idx...),
		pow2AutFixedPolyAutEvaluator:     e.pow2AutFixedPolyAutEvaluator.subEvaluator(idx...),
	}
}

// SafeCopy returns a thread-safe copy.
func (e *pow2AutFixedPolyEvaluator) SafeCopy() PolyEvaluator {
	return &pow2AutFixedPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.safeCopy(),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator,
		noReducePolyMulEvaluator:         e.noReducePolyMulEvaluator.safeCopy(),
		pow2AutFixedPolyAutEvaluator:     e.pow2AutFixedPolyAutEvaluator.safeCopy(),
	}
}

// primeAutFixedPolyEvaluator is a PolyEvaluator for prime order autfixed ring.
type primeAutFixedPolyEvaluator struct {
	polyBaseEvaluator
	primeAutFixedPolyScalarAddSubEvaluator
	noReducePolyMulEvaluator
	primeAutFixedPolyAutEvaluator
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *primeAutFixedPolyEvaluator) SubEvaluator(idx ...int) PolyEvaluator {
	return &primeAutFixedPolyEvaluator{
		polyBaseEvaluator:                      e.polyBaseEvaluator.subEvaluator(idx...),
		primeAutFixedPolyScalarAddSubEvaluator: e.primeAutFixedPolyScalarAddSubEvaluator.subEvaluator(idx...),
		noReducePolyMulEvaluator:               e.noReducePolyMulEvaluator.subEvaluator(idx...),
		primeAutFixedPolyAutEvaluator:          e.primeAutFixedPolyAutEvaluator.subEvaluator(idx...),
	}
}

// SafeCopy returns a thread-safe copy.
func (e *primeAutFixedPolyEvaluator) SafeCopy() PolyEvaluator {
	return &primeAutFixedPolyEvaluator{
		polyBaseEvaluator:                      e.polyBaseEvaluator.safeCopy(),
		primeAutFixedPolyScalarAddSubEvaluator: e.primeAutFixedPolyScalarAddSubEvaluator,
		noReducePolyMulEvaluator:               e.noReducePolyMulEvaluator.safeCopy(),
		primeAutFixedPolyAutEvaluator:          e.primeAutFixedPolyAutEvaluator.safeCopy(),
	}
}

// anyPolyEvaluator is a PolyEvaluator for aribtrary modulus polynomial.
type anyPolyEvaluator struct {
	polyBaseEvaluator
	defaultPolyScalarAddSubEvaluator
	reducePolyMulEvaluator
	noAutPolyAutEvaluator
}

// SubEvaluator returns a evaluator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *anyPolyEvaluator) SubEvaluator(idx ...int) PolyEvaluator {
	return &anyPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.subEvaluator(idx...),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator.subEvaluator(idx...),
		reducePolyMulEvaluator:           e.reducePolyMulEvaluator.subEvaluator(idx...),
		noAutPolyAutEvaluator:            e.noAutPolyAutEvaluator,
	}
}

// SafeCopy returns a thread-safe copy.
func (e *anyPolyEvaluator) SafeCopy() PolyEvaluator {
	return &anyPolyEvaluator{
		polyBaseEvaluator:                e.polyBaseEvaluator.safeCopy(),
		defaultPolyScalarAddSubEvaluator: e.defaultPolyScalarAddSubEvaluator,
		reducePolyMulEvaluator:           e.reducePolyMulEvaluator.safeCopy(),
		noAutPolyAutEvaluator:            e.noAutPolyAutEvaluator,
	}
}
