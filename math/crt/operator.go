package crt

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates ring operations over [Element].
//
// Operations usually take two forms: for example,
//   - Add(p0, p1) adds p0, p1, allocates a new vector to store the result and returns it.
//   - AddTo(pOut, p0, p1) adds p0, p1 and writes the result to pre-allocated pOut without returning.
//
// Moreover, operations panics when inputs are not consistent with the
// Operator's parameters, or operations itself are not valid.
type Operator interface {
	// Params returns the ring parameters.
	Params() dft.RingParameters
	// Modulus returns the modulus.
	Modulus() []*num.Modulus

	// NewPoly creates a new polynomial [Element] in Standard form.
	NewPoly() *Element
	// NewNTTPoly creates a new polynomial [Element] in NTT form.
	NewNTTPoly() *Element
	// NewPolyCustom creates a new polynomial [Element].
	NewPolyCustom(isNTT bool) *Element

	// FwdNTT returns FwdNTT(e).
	FwdNTT(e *Element) *Element
	// FwdNTTTo computes eOut = NTT(e).
	FwdNTTTo(eOut, e *Element)
	// InvNTT returns InvNTT(e).
	InvNTT(e *Element) *Element
	// InvNTTTo computes eOut = InvNTT(e).
	InvNTTTo(eOut, e *Element)

	addSubOperator

	// Neg returns -e.
	Neg(e *Element) *Element
	// NegTo computes eOut = -e.
	NegTo(eOut, e *Element)

	mulOperator

	autOperator

	// AsBig returns e as *[big.Int] vector.
	AsBig(e *Element) []*big.Int

	// SubOperator returns a operator for modulus of given indices.
	// Useful for "levelled" operations, especially with [vec.Range].
	SubOperator(idx ...int) Operator
}

// NewOperator creates a new [Operator].
func NewOperator(params dft.RingParameters, mod []*num.Modulus) Operator {
	if !isCoprime(mod) {
		panic("modulus must be coprime")
	}

	switch params.RingType() {
	case dft.TypeCyclotomic:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &pow2CyclotomicOperator{
				baseOperator:              newBaseOperator(params, mod),
				baseAddSubOperator:        newBaseAddSubOperator(params, mod),
				baseMulOperator:           newBaseMulOperator(params, mod),
				pow2CyclotomicAutOperator: newPow2CyclotomicAutOperator(params, mod),
			}
		default:
			reducer := NewCyclotomicReducer(params, mod)
			return &anyCyclotomicOperator{
				baseOperator:             newBaseOperator(params, mod),
				baseAddSubOperator:       newBaseAddSubOperator(params, mod),
				anyCyclotomicMulOperator: newAnyCyclotomicMulOperator(params, mod, reducer),
				anyCyclotomicAutOperator: newAnyCyclotomicAutOperator(params, mod, reducer),
			}
		}

	case dft.TypeCyclic:
		return &pow235CyclicOperator{
			baseOperator:       newBaseOperator(params, mod),
			baseAddSubOperator: newBaseAddSubOperator(params, mod),
			baseMulOperator:    newBaseMulOperator(params, mod),
			noAutOperator:      noAutOperator{},
		}

	case dft.TypeAutFixed:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &pow2AutFixedOperator{
				baseOperator:            newBaseOperator(params, mod),
				baseAddSubOperator:      newBaseAddSubOperator(params, mod),
				baseMulOperator:         newBaseMulOperator(params, mod),
				pow2AutFixedAutOperator: newPow2AutFixedAutOperator(params, mod),
			}
		case num.IsPrime(params.CycloOrder()):
			return &primeAutFixedOperator{
				baseOperator:                newBaseOperator(params, mod),
				primeAutFixedAddSubOperator: newPrimeAutFixedAddSubOperator(params, mod),
				baseMulOperator:             newBaseMulOperator(params, mod),
				primeAutFixedAutOperator:    newPrimeAutFixedAutOperator(params, mod),
			}
		}
	}

	panic("unsupported parameters")
}

// NewOperatorWithModPoly creates a new [Operator] with a given modulus polynomial.
func NewOperatorWithModPoly(mod []*num.Modulus, modPoly []int64) Operator {
	params := dft.NewOtherParameters(modPoly)

	reducer := NewReducer(num.NextProdPower(2*len(modPoly)-1, []int{2}), mod, modPoly)

	return &anyOperator{
		baseOperator:       newBaseOperator(params, mod),
		baseAddSubOperator: newBaseAddSubOperator(params, mod),
		reduceMulOperator:  newReduceMulOperator(mod, modPoly, reducer),
		noAutOperator:      noAutOperator{},
	}
}

// pow2CyclotomicOperator is a [Operator] for power-of-two cyclotomic ring.
type pow2CyclotomicOperator struct {
	baseOperator
	baseAddSubOperator
	baseMulOperator
	pow2CyclotomicAutOperator
}

// SubOperator returns a operator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *pow2CyclotomicOperator) SubOperator(idx ...int) Operator {
	return &pow2CyclotomicOperator{
		baseOperator:              e.baseOperator.subOperator(idx...),
		baseAddSubOperator:        e.baseAddSubOperator.subOperator(idx...),
		baseMulOperator:           e.baseMulOperator.subOperator(idx...),
		pow2CyclotomicAutOperator: e.pow2CyclotomicAutOperator.subOperator(idx...),
	}
}

// anyCyclotomicOperator is a [Operator] for arbitrary order cyclotomic ring.
type anyCyclotomicOperator struct {
	baseOperator
	baseAddSubOperator
	anyCyclotomicMulOperator
	anyCyclotomicAutOperator
}

// SubOperator returns a operator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *anyCyclotomicOperator) SubOperator(idx ...int) Operator {
	return &anyCyclotomicOperator{
		baseOperator:             e.baseOperator.subOperator(idx...),
		baseAddSubOperator:       e.baseAddSubOperator.subOperator(idx...),
		anyCyclotomicMulOperator: e.anyCyclotomicMulOperator.subOperator(idx...),
		anyCyclotomicAutOperator: e.anyCyclotomicAutOperator.subOperator(idx...),
	}
}

// pow235CyclicOperator is a [Operator] for cyclic ring with ranks multiple of 2, 3 and 5.
type pow235CyclicOperator struct {
	baseOperator
	baseAddSubOperator
	baseMulOperator
	noAutOperator
}

// SubOperator returns a operator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *pow235CyclicOperator) SubOperator(idx ...int) Operator {
	return &pow235CyclicOperator{
		baseOperator:       e.baseOperator.subOperator(idx...),
		baseAddSubOperator: e.baseAddSubOperator.subOperator(idx...),
		baseMulOperator:    e.baseMulOperator.subOperator(idx...),
		noAutOperator:      e.noAutOperator,
	}
}

// pow2AutFixedOperator is a [Operator] for power-of-two conjugate invariant ring.
type pow2AutFixedOperator struct {
	baseOperator
	baseAddSubOperator
	baseMulOperator
	pow2AutFixedAutOperator
}

// SubOperator returns a operator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *pow2AutFixedOperator) SubOperator(idx ...int) Operator {
	return &pow2AutFixedOperator{
		baseOperator:            e.baseOperator.subOperator(idx...),
		baseAddSubOperator:      e.baseAddSubOperator.subOperator(idx...),
		baseMulOperator:         e.baseMulOperator.subOperator(idx...),
		pow2AutFixedAutOperator: e.pow2AutFixedAutOperator.subOperator(idx...),
	}
}

// primeAutFixedOperator is a [Operator] for prime order autfixed ring.
type primeAutFixedOperator struct {
	baseOperator
	primeAutFixedAddSubOperator
	baseMulOperator
	primeAutFixedAutOperator
}

// SubOperator returns a operator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *primeAutFixedOperator) SubOperator(idx ...int) Operator {
	return &primeAutFixedOperator{
		baseOperator:                e.baseOperator.subOperator(idx...),
		primeAutFixedAddSubOperator: e.primeAutFixedAddSubOperator.subOperator(idx...),
		baseMulOperator:             e.baseMulOperator.subOperator(idx...),
		primeAutFixedAutOperator:    e.primeAutFixedAutOperator.subOperator(idx...),
	}
}

// anyOperator is a [Operator] for aribtrary modulus polynomial.
type anyOperator struct {
	baseOperator
	baseAddSubOperator
	reduceMulOperator
	noAutOperator
}

// SubOperator returns a operator for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (e *anyOperator) SubOperator(idx ...int) Operator {
	return &anyOperator{
		baseOperator:       e.baseOperator.subOperator(idx...),
		baseAddSubOperator: e.baseAddSubOperator.subOperator(idx...),
		reduceMulOperator:  e.reduceMulOperator.subOperator(idx...),
		noAutOperator:      e.noAutOperator,
	}
}
