package crt

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Operator evaluates ring operations over scalar and polynomial [Element] values.
//
// Operations usually take two forms: for example,
//   - Add(e0, e1) adds e0, e1, allocates a new element to store the result and returns it.
//   - AddTo(eOut, e0, e1) adds e0, e1 and writes the result to pre-allocated eOut without returning.
//
// Operations panic when inputs are not consistent with the
// Operator's parameters, or the operations themselves are not valid.
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
	if len(mod) == 0 {
		panic("empty modulus")
	} else if !IsCoprime(mod) {
		panic("modulus must be coprime")
	}

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
		modPoly := params.ModPoly()
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

// Mod returns the slice of *[num.Modulus] of the operator.
func (op *Operator) Mod() []*num.Modulus {
	return op.mod
}

// NewScalar creates a new scalar [Element].
func (op *Operator) NewScalar() *Element {
	return NewScalar(len(op.mod))
}

// NewPoly creates a new polynomial [Element].
func (op *Operator) NewPoly(form dft.Form) *Element {
	return NewPoly(op.params.Rank(), len(op.mod), form)
}

// FwdNTT returns FwdNTT(p).
//
// Panics if p is [TypeScalar].
func (op *Operator) FwdNTT(p *Element) *Element {
	pOut := op.NewPoly(dft.FormNTT)
	op.FwdNTTTo(pOut, p)
	return pOut
}

// FwdNTTTo computes pOut = NTT(p).
//
// Panics if p is [TypeScalar].
func (op *Operator) FwdNTTTo(pOut, p *Element) {
	checkUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if p.Form == dft.FormNTT {
		panic("inputs must be in Coeff form")
	} else if p.Type == TypeScalar {
		panic("inputs must be polynomial")
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
//
// Panics if p is [TypeScalar].
func (op *Operator) InvNTT(p *Element) *Element {
	pOut := op.NewPoly(dft.FormCoeff)
	op.InvNTTTo(pOut, p)
	return pOut
}

// InvNTTTo computes pOut = InvNTT(p).
//
// Panics if p is [TypeScalar].
func (op *Operator) InvNTTTo(pOut, p *Element) {
	checkUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if p.Form != dft.FormNTT {
		panic("inputs must be in NTT form")
	} else if p.Type == TypeScalar {
		panic("inputs must be polynomial")
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
func (op *Operator) Add(e0, e1 *Element) *Element {
	eOut := newBinaryOpOut(e0, e1)
	op.AddTo(eOut, e0, e1)
	return eOut
}

// AddTo computes eOut = e0 + e1.
func (op *Operator) AddTo(eOut, e0, e1 *Element) {
	switch op.addSubType {
	case typeTrivialAddSub:
		op.trivialAddSub.addTo(eOut, e0, e1)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSub.addTo(eOut, e0, e1)
	}
}

// Sub returns eOut = e0 - e1.
func (op *Operator) Sub(e0, e1 *Element) *Element {
	eOut := newBinaryOpOut(e0, e1)
	op.SubTo(eOut, e0, e1)
	return eOut
}

// SubTo computes eOut = e0 - e1.
func (op *Operator) SubTo(eOut, e0, e1 *Element) {
	switch op.addSubType {
	case typeTrivialAddSub:
		op.trivialAddSub.subTo(eOut, e0, e1)
	case typePrimeAutFixedAddSub:
		op.primeAutFixedAddSub.subTo(eOut, e0, e1)
	}
}

// Neg returns eOut = -e.
func (op *Operator) Neg(e *Element) *Element {
	eOut := newUnaryOpOut(e)
	op.NegTo(eOut, e)
	return eOut
}

// NegTo computes eOut = -e.
func (op *Operator) NegTo(eOut, e *Element) {
	checkUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	switch {
	case e.Type == TypeScalar:
		for i := range op.mod {
			eOut.Coeffs[i][0] = num.Neg(e.Coeffs[i][0], op.mod[i])
		}

	case e.Type == TypePoly:
		for i := range op.mod {
			vec.NegTo(eOut.Coeffs[i], e.Coeffs[i], op.mod[i])
		}
		eOut.Form = e.Form
	}
}

// Mul returns eOut = e0 * e1.
func (op *Operator) Mul(e0, e1 *Element) *Element {
	eOut := newBinaryOpOut(e0, e1)
	op.MulTo(eOut, e0, e1)
	return eOut
}

// MulTo computes eOut = e0 * e1.
func (op *Operator) MulTo(eOut, e0, e1 *Element) {
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
func (op *Operator) MulAddTo(eOut, e0, e1 *Element) {
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
func (op *Operator) MulSubTo(eOut, e0, e1 *Element) {
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
//
// Panics if p is [TypeScalar].
func (op *Operator) Aut(p *Element, idx int) *Element {
	pOut := op.NewPoly(p.Form)
	op.AutTo(pOut, p, idx)
	return pOut
}

// AutTo computes pOut = p(X^idx).
//
// Panics if p is [TypeScalar].
func (op *Operator) AutTo(pOut, p *Element, idx int) {
	if p.Type == TypeScalar {
		panic("inputs must be polynomial")
	}

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

// WithModSlice slices the modulus of [Operator] and returns the new [Operator].
func (op *Operator) WithModSlice(lo, hi int) *Operator {
	if lo == 0 && hi == len(op.mod) {
		return op
	}

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

// Append returns a new [Operator] with the moduli of op followed by those of opIn.
//
// Panics when the ring parameters do not match or the resulting moduli are not pairwise coprime.
func (op *Operator) Append(opIn *Operator) *Operator {
	if !op.params.IsEqual(opIn.params) {
		panic("inconsistent inputs")
	}

	modOut := slices.Concat(op.mod, opIn.mod)
	if !IsCoprime(modOut) {
		panic("modulus must be coprime")
	}

	opOut := &Operator{
		params: op.params,
		mod:    modOut,

		isNTTFriendly: slices.Concat(op.isNTTFriendly, opIn.isNTTFriendly),
		ntt:           slices.Concat(op.ntt, opIn.ntt),

		addSubType: op.addSubType,
		mulType:    op.mulType,
		autType:    op.autType,
	}

	switch opOut.addSubType {
	case typeTrivialAddSub:
		opOut.trivialAddSub = op.trivialAddSub.append(opIn.trivialAddSub)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSub = op.primeAutFixedAddSub.append(opIn.primeAutFixedAddSub)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMul = op.trivialMul.append(opIn.trivialMul)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMul = op.anyCyclotomicMul.append(opIn.anyCyclotomicMul)
	case typeReduceMul:
		opOut.reduceMul = op.reduceMul.append(opIn.reduceMul)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAut = op.pow2CyclotomicAut.append(opIn.pow2CyclotomicAut)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAut = op.anyCyclotomicAut.append(opIn.anyCyclotomicAut)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAut = op.pow2AutFixedAut.append(opIn.pow2AutFixedAut)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAut = op.primeAutFixedAut.append(opIn.primeAutFixedAut)
	}

	return opOut
}

// AppendTmpMod appends "temporary" modulus to the moduli chain and returns the new [Operator].
// This assumes that modulus is NTT-unfriendly, trading the appending performance with operation performance.
//
// Panics when the resulting moduli are not pairwise coprime.
func (op *Operator) AppendTmpMod(mod *num.Modulus) *Operator {
	modOut := slices.Concat(op.mod, []*num.Modulus{mod})
	if !IsCoprime(modOut) {
		panic("modulus must be coprime")
	}

	opOut := &Operator{
		params: op.params,
		mod:    modOut,

		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),
		ntt:           slices.Concat(op.ntt, []*dft.Transformer{nil}),

		addSubType: op.addSubType,
		mulType:    op.mulType,
		autType:    op.autType,
	}

	switch opOut.addSubType {
	case typeTrivialAddSub:
		opOut.trivialAddSub = op.trivialAddSub.appendTmpMod(mod)
	case typePrimeAutFixedAddSub:
		opOut.primeAutFixedAddSub = op.primeAutFixedAddSub.appendTmpMod(mod)
	}

	switch opOut.mulType {
	case typeTrivialMul:
		opOut.trivialMul = op.trivialMul.appendTmpMod(mod)
	case typeAnyCyclotomicMul:
		opOut.anyCyclotomicMul = op.anyCyclotomicMul.appendTmpMod(mod)
	case typeReduceMul:
		opOut.reduceMul = op.reduceMul.appendTmpMod(mod)
	}

	switch opOut.autType {
	case typePow2CyclotomicAut:
		opOut.pow2CyclotomicAut = op.pow2CyclotomicAut.appendTmpMod(mod)
	case typeAnyCyclotomicAut:
		opOut.anyCyclotomicAut = op.anyCyclotomicAut.appendTmpMod(mod)
	case typePow2AutFixedAut:
		opOut.pow2AutFixedAut = op.pow2AutFixedAut.appendTmpMod(mod)
	case typePrimeAutFixedAut:
		opOut.primeAutFixedAut = op.primeAutFixedAut.appendTmpMod(mod)
	}

	return opOut
}
