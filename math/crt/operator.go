package crt

import (
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
type Operator struct {
	*baseOperator
	addSubOperator
	mulOperator
	autOperator
}

// NewOperator creates a new [Operator].
func NewOperator(params dft.RingParameters, mod []*num.Modulus) *Operator {
	if !isCoprime(mod) {
		panic("modulus must be coprime")
	}

	switch params.RingType() {
	case dft.TypeCyclotomic:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &Operator{
				baseOperator:   newBaseOperator(params, mod),
				addSubOperator: newBaseAddSubOperator(params, mod),
				mulOperator:    newBaseMulOperator(params, mod),
				autOperator:    newPow2CyclotomicAutOperator(params, mod),
			}
		default:
			reducer := NewCyclotomicReducer(params, mod)
			return &Operator{
				baseOperator:   newBaseOperator(params, mod),
				addSubOperator: newBaseAddSubOperator(params, mod),
				mulOperator:    newAnyCyclotomicMulOperator(params, mod, reducer),
				autOperator:    newAnyCyclotomicAutOperator(params, mod, reducer),
			}
		}

	case dft.TypeCyclic:
		return &Operator{
			baseOperator:   newBaseOperator(params, mod),
			addSubOperator: newBaseAddSubOperator(params, mod),
			mulOperator:    newBaseMulOperator(params, mod),
			autOperator:    noAutOperator{},
		}

	case dft.TypeAutFixed:
		switch {
		case num.IsPowerOfTwo(params.CycloOrder()):
			return &Operator{
				baseOperator:   newBaseOperator(params, mod),
				addSubOperator: newBaseAddSubOperator(params, mod),
				mulOperator:    newBaseMulOperator(params, mod),
				autOperator:    newPow2AutFixedAutOperator(params, mod),
			}
		case num.IsPrime(params.CycloOrder()):
			return &Operator{
				baseOperator:   newBaseOperator(params, mod),
				addSubOperator: newPrimeAutFixedAddSubOperator(params, mod),
				mulOperator:    newBaseMulOperator(params, mod),
				autOperator:    newPrimeAutFixedAutOperator(params, mod),
			}
		}

	case dft.TypeOther:
		modPoly := params.ModulusPoly()
		reducer := NewReducer(num.NextProdPower(2*len(modPoly)-1, []int{2}), mod, modPoly)

		return &Operator{
			baseOperator:   newBaseOperator(params, mod),
			addSubOperator: newBaseAddSubOperator(params, mod),
			mulOperator:    newReduceMulOperator(mod, modPoly, reducer),
			autOperator:    noAutOperator{},
		}
	}

	panic("unsupported parameters")
}

// WithModIdx returns an [Operator] for modulus of given indices.
// Useful for "levelled" operations, especially with [vec.Range].
func (op *Operator) WithModIdx(idx ...int) *Operator {
	return &Operator{
		baseOperator:   op.baseOperator.withModIdx(idx...),
		addSubOperator: op.addSubOperator.withModIdx(idx...),
		mulOperator:    op.mulOperator.withModIdx(idx...),
		autOperator:    op.autOperator.withModIdx(idx...),
	}
}

// Append appends a new [Operator] and returns the new [Operator].
func (op *Operator) Append(op0 *Operator) *Operator {
	return &Operator{
		baseOperator:   op.baseOperator.append(op0.baseOperator),
		addSubOperator: op.addSubOperator.append(op0.addSubOperator),
		mulOperator:    op.mulOperator.append(op0.mulOperator),
		autOperator:    op.autOperator.append(op0.autOperator),
	}
}

// AppendAuxModulus appends "auxillary" modulus to the moduli chain and returns the new [Operator].
// This assumes that modulus is NTT-unfriendly, trading the appending performance with operation performance.
func (op *Operator) AppendAuxModulus(mod *num.Modulus) *Operator {
	return &Operator{
		baseOperator:   op.baseOperator.appendAuxModulus(mod),
		addSubOperator: op.addSubOperator.appendAuxModulus(mod),
		mulOperator:    op.mulOperator.appendAuxModulus(mod),
		autOperator:    op.autOperator.appendAuxModulus(mod),
	}
}
