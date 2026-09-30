package crt

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Operator evaluates ring operations over [Scalar] and [Poly].
//
// Operations usually take two forms: for example,
//   - Add(p0, p1) adds p0, p1, allocates a new vector to store the result and returns it.
//   - AddTo(pOut, p0, p1) adds p0, p1 and writes the result to pre-allocated pOut without returning.
//
// Moreover, operations panics when inputs are not consistent with the
// Operator's parameters, or operations itself are not valid.
type Operator struct {
	params dft.RingParameters
	mod    []*num.Modulus

	isNTTFriendly []bool
	ntt           []*dft.Transformer

	addSubType addSubType
	*trivialAddSubOperator
	*primeAutFixedAddSubOperator

	mulType mulType
	*trivialMulOperator
	*anyCyclotomicMulOperator
	*reduceMulOperator

	autType autType
	*pow2CyclotomicAutOperator
	*anyCyclotomicAutOperator
	*pow2AutFixedAutOperator
	*primeAutFixedAutOperator
}

// NewOperator creates a new [Operator].
func NewOperator(params dft.RingParameters, mod []*num.Modulus) *Operator {
	isNTTFriendly := make([]bool, len(mod))
	ntt := make([]*dft.Transformer, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
		if isNTTFriendly[i] {
			ntt[i] = dft.NewTransformer(params, mod[i])
		}
	}

	op := &Operator{
		params: params,
		mod:    mod,

		isNTTFriendly: isNTTFriendly,
		ntt:           ntt,

		addSubType: addSubTypeOf(params),
		mulType:    mulTypeOf(params),
		autType:    autTypeOf(params),
	}

	var reducer *CyclotomicReducer
	if params.RingType() == dft.TypeCyclotomic && !num.IsPowerOfTwo(params.CycloIndex()) {
		reducer = NewCyclotomicReducer(params, mod)
	}

	switch op.addSubType {
	case typeTrivialAddSub:
		op.trivialAddSubOperator = newTrivialAddSubOperator(params, mod)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSubOperator = newPrimeAutFixedAddSubOperator(params, mod)
	}

	switch op.mulType {
	case typeTrivialMul:
		op.trivialMulOperator = newTrivialMulOperator(params, mod)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMulOperator = newAnyCyclotomicMulOperator(params, mod, reducer)
	case typeReduceMul:
		modPoly := params.ModulusPoly()
		maxRank := num.NextProdPower(2*len(modPoly)-1, []int{2})
		op.reduceMulOperator = newReduceMulOperator(mod, modPoly, NewReducer(maxRank, mod, modPoly))
	}

	switch op.autType {
	case typePow2CyclotomicAut:
		op.pow2CyclotomicAutOperator = newPow2CyclotomicAutOperator(params, mod)
	case typeAnyCyclotomicAut:
		op.anyCyclotomicAutOperator = newAnyCyclotomicAutOperator(params, mod, reducer)
	case typePow2AutFixedAut:
		op.pow2AutFixedAutOperator = newPow2AutFixedAutOperator(params, mod)
	case typePrimeAutFixedAut:
		op.primeAutFixedAutOperator = newPrimeAutFixedAutOperator(params, mod)
	}

	return op
}

// Params returns the ring parameters.
func (op *Operator) Params() dft.RingParameters {
	return op.params
}

// Modulus returns the modulus.
func (op *Operator) Modulus() []*num.Modulus {
	return op.mod
}

// NewPoly creates a new [Poly] in Coefficient form.
func (op *Operator) NewPoly() *Poly {
	return NewPoly(op.params.Rank(), len(op.mod))
}

// NewNTTPoly creates a new [Poly] in NTT form.
func (op *Operator) NewNTTPoly() *Poly {
	return NewNTTPoly(op.params.Rank(), len(op.mod))
}

// NewPolyCustom creates a new [Poly].
func (op *Operator) NewPolyCustom(form dft.Form) *Poly {
	return NewPolyCustom(op.params.Rank(), len(op.mod), form)
}

// FwdNTT returns FwdNTT(p).
func (op *Operator) FwdNTT(p *Poly) *Poly {
	pOut := op.NewPoly()
	op.FwdNTTTo(pOut, p)
	return pOut
}

// FwdNTTTo computes pOut = NTT(p).
func (op *Operator) FwdNTTTo(pOut, p *Poly) {
	isUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if p.Form == dft.FormNTT {
		panic("input(s) must be in coefficient form")
	}

	for i := range op.ntt {
		if op.ntt[i] != nil {
			op.ntt[i].ForwardTo(pOut.Coeffs[i], p.Coeffs[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
		}
	}

	pOut.Form = dft.FormNTT
}

// InvNTT returns InvNTT(p).
func (op *Operator) InvNTT(p *Poly) *Poly {
	pOut := op.NewPoly()
	op.InvNTTTo(pOut, p)
	return pOut
}

// InvNTTTo computes pOut = InvNTT(p).
func (op *Operator) InvNTTTo(pOut, p *Poly) {
	isUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if p.Form != dft.FormNTT {
		panic("input(s) must be in NTT form")
	}

	for i := range op.ntt {
		if op.ntt[i] != nil {
			op.ntt[i].InverseTo(pOut.Coeffs[i], p.Coeffs[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
		}
	}

	pOut.Form = dft.FormCoeff
}

// Add returns eOut = e0 + e1.
func (op *Operator) Add[TOut, T0, T1 *Scalar | *Poly](e0 T0, e1 T1) TOut {
	eOut := newBinaryOpOut[TOut](e0, e1)
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *Operator) AddTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.addSubType {
	case typeTrivialAddSub:
		op.trivialAddSubOperator.addTo(eOut, e0, e1)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSubOperator.addTo(eOut, e0, e1)
	}
}

// Sub returns eOut = e0 - e1.
func (op *Operator) Sub[TOut, T0, T1 *Scalar | *Poly](e0 T0, e1 T1) TOut {
	eOut := newBinaryOpOut[TOut](e0, e1)
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *Operator) SubTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.addSubType {
	case typeTrivialAddSub:
		op.trivialAddSubOperator.subTo(eOut, e0, e1)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSubOperator.subTo(eOut, e0, e1)
	}
}

// Neg returns eOut = -e.
func (op *Operator) Neg[TOut, T *Scalar | *Poly](e T) TOut {
	eOut := newUnaryOpOut[TOut](e)
	op.NegTo(eOut, e)
	return eOut
}

// NegTo computes eOut = -e.
func (op *Operator) NegTo[TOut, T *Scalar | *Poly](eOut TOut, e T) {
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	switch {
	case eIsScalar:
		for i := range op.mod {
			eOutScalar.Value[i] = num.Neg(eScalar.Value[i], op.mod[i])
		}

	case !eIsScalar:
		for i := range op.mod {
			vec.NegTo(eOutPoly.Coeffs[i], ePoly.Coeffs[i], op.mod[i])
		}
		eOutPoly.Form = ePoly.Form
	}
}

// Mul returns eOut = e0 * e1.
func (op *Operator) Mul[TOut, T0, T1 *Scalar | *Poly](e0 T0, e1 T1) TOut {
	eOut := newBinaryOpOut[TOut](e0, e1)
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
func (op *Operator) MulTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.mulType {
	case typeTrivialMul:
		op.trivialMulOperator.mulTo(eOut, e0, e1)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMulOperator.mulTo(eOut, e0, e1)
	case typeReduceMul:
		op.reduceMulOperator.mulTo(eOut, e0, e1)
	}
}

// MulAddTo computes eOut += e0 * e1.
func (op *Operator) MulAddTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.mulType {
	case typeTrivialMul:
		op.trivialMulOperator.mulAddTo(eOut, e0, e1)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMulOperator.mulAddTo(eOut, e0, e1)
	case typeReduceMul:
		op.reduceMulOperator.mulAddTo(eOut, e0, e1)
	}
}

// MulSubTo computes eOut -= e0 * e1.
func (op *Operator) MulSubTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.mulType {
	case typeTrivialMul:
		op.trivialMulOperator.mulSubTo(eOut, e0, e1)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMulOperator.mulSubTo(eOut, e0, e1)
	case typeReduceMul:
		op.reduceMulOperator.mulSubTo(eOut, e0, e1)
	}
}

// CanAut returns if automorphism X -> X^idx is valid.
func (op *Operator) CanAut(idx int) bool {
	switch op.autType {
	case typePow2CyclotomicAut:
		return op.pow2CyclotomicAutOperator.canAut(idx)
	case typeAnyCyclotomicAut:
		return op.anyCyclotomicAutOperator.canAut(idx)
	case typePow2AutFixedAut:
		return op.pow2AutFixedAutOperator.canAut(idx)
	case typePrimeAutFixedAut:
		return op.primeAutFixedAutOperator.canAut(idx)
	}
	return false
}

// Aut returns pOut = p(X^idx).
func (op *Operator) Aut(p *Poly, idx int) *Poly {
	pOut := op.NewPoly()
	op.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = p(X^idx).
func (op *Operator) AutTo(pOut, p *Poly, idx int) {
	switch op.autType {
	case typePow2CyclotomicAut:
		op.pow2CyclotomicAutOperator.autTo(pOut, p, idx)
	case typeAnyCyclotomicAut:
		op.anyCyclotomicAutOperator.autTo(pOut, p, idx)
	case typePow2AutFixedAut:
		op.pow2AutFixedAutOperator.autTo(pOut, p, idx)
	case typePrimeAutFixedAut:
		op.primeAutFixedAutOperator.autTo(pOut, p, idx)
	default:
		panic("automorphism unsupported")
	}
}

// WithModIdx returns an [Operator] for modulus of given indices.
// Useful for "levelled" operations with arbitrary modulus indices.
// Use [Operator.Slice] for contiguous ranges.
func (op *Operator) WithModIdx(idx ...int) *Operator {
	opOut := &Operator{
		params: op.params,
		mod:    vec.Gather(op.mod, idx...),

		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),
		ntt:           vec.Gather(op.ntt, idx...),

		addSubType: op.addSubType,
		mulType:    op.mulType,
		autType:    op.autType,
	}

	switch opOut.addSubType {
	case typeTrivialAddSub:
		opOut.trivialAddSubOperator = op.trivialAddSubOperator.withModIdx(idx...)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSubOperator = op.primeAutFixedAddSubOperator.withModIdx(idx...)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMulOperator = op.trivialMulOperator.withModIdx(idx...)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMulOperator = op.anyCyclotomicMulOperator.withModIdx(idx...)
	case typeReduceMul:
		opOut.reduceMulOperator = op.reduceMulOperator.withModIdx(idx...)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAutOperator = op.pow2CyclotomicAutOperator.withModIdx(idx...)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAutOperator = op.anyCyclotomicAutOperator.withModIdx(idx...)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAutOperator = op.pow2AutFixedAutOperator.withModIdx(idx...)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAutOperator = op.primeAutFixedAutOperator.withModIdx(idx...)
	}

	return opOut
}

// Slice slices the modulus of [Operator] and returns the new [Operator].
// Equal to op.WithModIdx(vec.Range(lo, hi)...).
func (op *Operator) Slice(lo, hi int) *Operator {
	opOut := &Operator{
		params: op.params,
		mod:    op.mod[lo:hi],

		isNTTFriendly: op.isNTTFriendly[lo:hi],
		ntt:           op.ntt[lo:hi],

		addSubType: op.addSubType,
		mulType:    op.mulType,
		autType:    op.autType,
	}

	switch opOut.addSubType {
	case typeTrivialAddSub:
		opOut.trivialAddSubOperator = op.trivialAddSubOperator.slice(lo, hi)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSubOperator = op.primeAutFixedAddSubOperator.slice(lo, hi)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMulOperator = op.trivialMulOperator.slice(lo, hi)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMulOperator = op.anyCyclotomicMulOperator.slice(lo, hi)
	case typeReduceMul:
		opOut.reduceMulOperator = op.reduceMulOperator.slice(lo, hi)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAutOperator = op.pow2CyclotomicAutOperator.slice(lo, hi)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAutOperator = op.anyCyclotomicAutOperator.slice(lo, hi)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAutOperator = op.pow2AutFixedAutOperator.slice(lo, hi)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAutOperator = op.primeAutFixedAutOperator.slice(lo, hi)
	}

	return opOut
}

// Append appends a new [Operator] and returns the new [Operator].
//
// Panics when the ring parameters do not match.
func (op *Operator) Append(op0 *Operator) *Operator {
	if !op.params.IsEqual(op0.params) {
		panic("inconsistent input(s)")
	}

	opOut := &Operator{
		params: op.params,
		mod:    slices.Concat(op.mod, op0.mod),

		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),
		ntt:           slices.Concat(op.ntt, op0.ntt),

		addSubType: op.addSubType,
		mulType:    op.mulType,
		autType:    op.autType,
	}

	switch opOut.addSubType {
	case typeTrivialAddSub:
		opOut.trivialAddSubOperator = op.trivialAddSubOperator.append(op0.trivialAddSubOperator)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSubOperator = op.primeAutFixedAddSubOperator.append(op0.primeAutFixedAddSubOperator)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMulOperator = op.trivialMulOperator.append(op0.trivialMulOperator)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMulOperator = op.anyCyclotomicMulOperator.append(op0.anyCyclotomicMulOperator)
	case typeReduceMul:
		opOut.reduceMulOperator = op.reduceMulOperator.append(op0.reduceMulOperator)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAutOperator = op.pow2CyclotomicAutOperator.append(op0.pow2CyclotomicAutOperator)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAutOperator = op.anyCyclotomicAutOperator.append(op0.anyCyclotomicAutOperator)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAutOperator = op.pow2AutFixedAutOperator.append(op0.pow2AutFixedAutOperator)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAutOperator = op.primeAutFixedAutOperator.append(op0.primeAutFixedAutOperator)
	}

	return opOut
}

// AppendTmpModulus appends "temporary" modulus to the moduli chain and returns the new [Operator].
// This assumes that modulus is NTT-unfriendly, trading the appending performance with operation performance.
func (op *Operator) AppendTmpModulus(mod *num.Modulus) *Operator {
	opOut := &Operator{
		params: op.params,
		mod:    slices.Concat(op.mod, []*num.Modulus{mod}),

		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),
		ntt:           slices.Concat(op.ntt, []*dft.Transformer{nil}),

		addSubType: op.addSubType,
		mulType:    op.mulType,
		autType:    op.autType,
	}

	switch opOut.addSubType {
	case typeTrivialAddSub:
		opOut.trivialAddSubOperator = op.trivialAddSubOperator.appendTmpModulus(mod)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSubOperator = op.primeAutFixedAddSubOperator.appendTmpModulus(mod)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMulOperator = op.trivialMulOperator.appendTmpModulus(mod)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMulOperator = op.anyCyclotomicMulOperator.appendTmpModulus(mod)
	case typeReduceMul:
		opOut.reduceMulOperator = op.reduceMulOperator.appendTmpModulus(mod)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAutOperator = op.pow2CyclotomicAutOperator.appendTmpModulus(mod)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAutOperator = op.anyCyclotomicAutOperator.appendTmpModulus(mod)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAutOperator = op.pow2AutFixedAutOperator.appendTmpModulus(mod)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAutOperator = op.primeAutFixedAutOperator.appendTmpModulus(mod)
	}

	return opOut
}
