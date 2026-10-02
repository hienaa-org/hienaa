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

	addSubType          addSubType
	trivialAddSub       *trivialAddSubOperator
	primeAutFixedAddSub *primeAutFixedAddSubOperator

	mulType          mulType
	trivialMul       *trivialMulOperator
	anyCyclotomicMul *anyCyclotomicMulOperator
	reduceMul        *reduceMulOperator

	autType           autType
	pow2CyclotomicAut *pow2CyclotomicAutOperator
	anyCyclotomicAut  *anyCyclotomicAutOperator
	pow2AutFixedAut   *pow2AutFixedAutOperator
	primeAutFixedAut  *primeAutFixedAutOperator
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
		op.trivialAddSub = newTrivialAddSubOperator(params, mod)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSub = newPrimeAutFixedAddSubOperator(params, mod)
	}

	switch op.mulType {
	case typeTrivialMul:
		op.trivialMul = newTrivialMulOperator(params, mod)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMul = newAnyCyclotomicMulOperator(params, mod, reducer)
	case typeReduceMul:
		modPoly := params.ModulusPoly()
		maxRank := num.NextProdPower(2*len(modPoly)-1, []int{2})
		op.reduceMul = newReduceMulOperator(mod, modPoly, NewReducer(maxRank, mod, modPoly))
	}

	switch op.autType {
	case typePow2CyclotomicAut:
		op.pow2CyclotomicAut = newPow2CyclotomicAutOperator(params, mod)
	case typeAnyCyclotomicAut:
		op.anyCyclotomicAut = newAnyCyclotomicAutOperator(params, mod, reducer)
	case typePow2AutFixedAut:
		op.pow2AutFixedAut = newPow2AutFixedAutOperator(params, mod)
	case typePrimeAutFixedAut:
		op.primeAutFixedAut = newPrimeAutFixedAutOperator(params, mod)
	}

	return op
}

// Params returns the [dft.RingParameters] of the operator.
func (op *Operator) Params() dft.RingParameters {
	return op.params
}

// Modulus returns the slice of *[num.Modulus] of the operator.
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
		panic("input(s) must be in Coeff form")
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
		op.trivialAddSub.addTo(eOut, e0, e1)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSub.addTo(eOut, e0, e1)
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
		op.trivialAddSub.subTo(eOut, e0, e1)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSub.subTo(eOut, e0, e1)
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
		op.trivialMul.mulTo(eOut, e0, e1)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMul.mulTo(eOut, e0, e1)
	case typeReduceMul:
		op.reduceMul.mulTo(eOut, e0, e1)
	}
}

// MulAddTo computes eOut += e0 * e1.
func (op *Operator) MulAddTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.mulType {
	case typeTrivialMul:
		op.trivialMul.mulAddTo(eOut, e0, e1)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMul.mulAddTo(eOut, e0, e1)
	case typeReduceMul:
		op.reduceMul.mulAddTo(eOut, e0, e1)
	}
}

// MulSubTo computes eOut -= e0 * e1.
func (op *Operator) MulSubTo[TOut, T0, T1 *Scalar | *Poly](eOut TOut, e0 T0, e1 T1) {
	switch op.mulType {
	case typeTrivialMul:
		op.trivialMul.mulSubTo(eOut, e0, e1)
	case typeAnyCyclotomicMul:
		op.anyCyclotomicMul.mulSubTo(eOut, e0, e1)
	case typeReduceMul:
		op.reduceMul.mulSubTo(eOut, e0, e1)
	}
}

// CanAut returns if automorphism X -> X^idx is valid.
func (op *Operator) CanAut(idx int) bool {
	switch op.autType {
	case typePow2CyclotomicAut:
		return op.pow2CyclotomicAut.canAut(idx)
	case typeAnyCyclotomicAut:
		return op.anyCyclotomicAut.canAut(idx)
	case typePow2AutFixedAut:
		return op.pow2AutFixedAut.canAut(idx)
	case typePrimeAutFixedAut:
		return op.primeAutFixedAut.canAut(idx)
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
		op.pow2CyclotomicAut.autTo(pOut, p, idx)
	case typeAnyCyclotomicAut:
		op.anyCyclotomicAut.autTo(pOut, p, idx)
	case typePow2AutFixedAut:
		op.pow2AutFixedAut.autTo(pOut, p, idx)
	case typePrimeAutFixedAut:
		op.primeAutFixedAut.autTo(pOut, p, idx)
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
		opOut.trivialAddSub = op.trivialAddSub.withModIdx(idx...)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSub = op.primeAutFixedAddSub.withModIdx(idx...)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMul = op.trivialMul.withModIdx(idx...)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMul = op.anyCyclotomicMul.withModIdx(idx...)
	case typeReduceMul:
		opOut.reduceMul = op.reduceMul.withModIdx(idx...)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAut = op.pow2CyclotomicAut.withModIdx(idx...)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAut = op.anyCyclotomicAut.withModIdx(idx...)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAut = op.pow2AutFixedAut.withModIdx(idx...)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAut = op.primeAutFixedAut.withModIdx(idx...)
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
		opOut.trivialAddSub = op.trivialAddSub.slice(lo, hi)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSub = op.primeAutFixedAddSub.slice(lo, hi)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMul = op.trivialMul.slice(lo, hi)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMul = op.anyCyclotomicMul.slice(lo, hi)
	case typeReduceMul:
		opOut.reduceMul = op.reduceMul.slice(lo, hi)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAut = op.pow2CyclotomicAut.slice(lo, hi)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAut = op.anyCyclotomicAut.slice(lo, hi)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAut = op.pow2AutFixedAut.slice(lo, hi)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAut = op.primeAutFixedAut.slice(lo, hi)
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
		opOut.trivialAddSub = op.trivialAddSub.append(op0.trivialAddSub)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSub = op.primeAutFixedAddSub.append(op0.primeAutFixedAddSub)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMul = op.trivialMul.append(op0.trivialMul)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMul = op.anyCyclotomicMul.append(op0.anyCyclotomicMul)
	case typeReduceMul:
		opOut.reduceMul = op.reduceMul.append(op0.reduceMul)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAut = op.pow2CyclotomicAut.append(op0.pow2CyclotomicAut)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAut = op.anyCyclotomicAut.append(op0.anyCyclotomicAut)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAut = op.pow2AutFixedAut.append(op0.pow2AutFixedAut)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAut = op.primeAutFixedAut.append(op0.primeAutFixedAut)
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
		opOut.trivialAddSub = op.trivialAddSub.appendTmpModulus(mod)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSub = op.primeAutFixedAddSub.appendTmpModulus(mod)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMul = op.trivialMul.appendTmpModulus(mod)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMul = op.anyCyclotomicMul.appendTmpModulus(mod)
	case typeReduceMul:
		opOut.reduceMul = op.reduceMul.appendTmpModulus(mod)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAut = op.pow2CyclotomicAut.appendTmpModulus(mod)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAut = op.anyCyclotomicAut.appendTmpModulus(mod)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAut = op.pow2AutFixedAut.appendTmpModulus(mod)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAut = op.primeAutFixedAut.appendTmpModulus(mod)
	}

	return opOut
}
