package rlwe

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// Element represents an element with CRT representation.
// Can be a scalar or a polynomial.
//
// A polynomial can have Standard or NTT form.
type Element struct {
	// Value is the underlying [crt.Element].
	Value *crt.Element
	// auxLen is the auxiliary modulus length.
	auxLen int
}

// NewScalar creates a new scalar [Element].
func NewScalar(params Parameters, withAux bool) *Element {
	if withAux && params.AuxModLen() == 0 {
		panic("no auxiliary modulus")
	}

	if withAux {
		return NewScalarCustom(params.BaseModLen(), params.AuxModLen())
	}
	return NewScalarCustom(params.BaseModLen(), 0)
}

// NewScalarCustom creates a new scalar [Element] with the given parameters.
func NewScalarCustom(baseLen, auxLen int) *Element {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	return &Element{
		Value:  crt.NewScalar(baseLen + auxLen),
		auxLen: auxLen,
	}
}

// NewPoly creates a new polynomial [Element].
func NewPoly(params Parameters, withAux bool, form dft.Form) *Element {
	if withAux && params.AuxModLen() == 0 {
		panic("no auxiliary modulus")
	}

	if withAux {
		return NewPolyCustom(params.Rank(), params.BaseModLen(), params.AuxModLen(), form)
	}
	return NewPolyCustom(params.Rank(), params.BaseModLen(), 0, form)
}

// NewPolyCustom creates a new polynomial [Element] with the given parameters.
func NewPolyCustom(rank, baseLen, auxLen int, form dft.Form) *Element {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	return &Element{
		Value:  crt.NewPoly(rank, baseLen+auxLen, form),
		auxLen: auxLen,
	}
}

// NewElement creates a new [Element].
// Rank and form are ignored when the type is [crt.TypeScalar].
func NewElement(params Parameters, withAux bool, elType crt.ElementType, form dft.Form) *Element {
	if elType == crt.TypeScalar {
		return NewScalar(params, withAux)
	}
	return NewPoly(params, withAux, form)
}

// NewElementCustom creates a new [Element] with the given parameters.
// Rank and form are ignored when the type is [crt.TypeScalar].
func NewElementCustom(rank, baseLen, auxLen int, elType crt.ElementType, form dft.Form) *Element {
	if elType == crt.TypeScalar {
		return NewScalarCustom(baseLen, auxLen)
	}
	return NewPolyCustom(rank, baseLen, auxLen, form)
}

// NewElementFrom creates a new [Element] from [crt.Element].
func NewElementFrom(value *crt.Element, auxLen int) *Element {
	if auxLen < 0 || auxLen > value.ModLen() {
		panic("invalid auxLen")
	}

	return &Element{
		Value:  value,
		auxLen: auxLen,
	}
}

// Rank returns the rank of e.
func (e *Element) Rank() int {
	return e.Value.Rank()
}

// Form returns the [dft.Form] of e.
func (e *Element) Form() dft.Form {
	return e.Value.Form
}

// SetForm sets the [dft.Form] of e.
// This does not change the coefficients, so it may invalidate correctness.
func (e *Element) SetForm(form dft.Form) {
	e.Value.Form = form
}

// Type returns the [crt.ElementType] of e.
func (e *Element) Type() crt.ElementType {
	return e.Value.Type
}

// FullModLen returns the full modulus length.
func (e *Element) FullModLen() int {
	return e.Value.ModLen()
}

// BaseModLen returns the base modulus length.
func (e *Element) BaseModLen() int {
	return e.Value.ModLen() - e.auxLen
}

// AuxModLen returns the auxiliary modulus length.
func (e *Element) AuxModLen() int {
	return e.auxLen
}

// Clear clears e.
func (e *Element) Clear() {
	e.Value.Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (e *Element) WithModLen(baseLen, auxLen int) *Element {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	} else if baseLen > e.BaseModLen() || auxLen > e.AuxModLen() {
		panic("out of range")
	}

	return &Element{
		Value:  e.Value.WithModSlice(e.auxLen-auxLen, e.auxLen+baseLen),
		auxLen: auxLen,
	}
}

// Resize resizes e in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is invalid.
func (e *Element) Resize(baseLen, auxLen int) {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	rank := e.Rank()
	auxLenKeep := min(e.AuxModLen(), auxLen)
	baseLenKeep := min(e.BaseModLen(), baseLen)

	if auxLen+baseLen > len(e.Value.Coeffs) {
		e.Value.Coeffs = slices.Grow(e.Value.Coeffs, auxLen+baseLen-len(e.Value.Coeffs))[:auxLen+baseLen]
	}

	for i := range e.Value.Coeffs {
		if e.Value.Coeffs[i] == nil {
			e.Value.Coeffs[i] = make([]uint64, rank)
		}
	}

	switch {
	case auxLen < e.auxLen:
		for i := -auxLenKeep; i < baseLenKeep; i++ {
			copy(e.Value.Coeffs[auxLen+i], e.Value.Coeffs[e.auxLen+i])
		}
	case auxLen > e.auxLen:
		for i := baseLenKeep - 1; i >= -auxLenKeep; i-- {
			copy(e.Value.Coeffs[auxLen+i], e.Value.Coeffs[e.auxLen+i])
		}
	}

	for i := 0; i < auxLen-auxLenKeep; i++ {
		clear(e.Value.Coeffs[i])
	}
	for i := auxLen + baseLenKeep; i < len(e.Value.Coeffs); i++ {
		clear(e.Value.Coeffs[i])
	}

	e.Value.Coeffs = e.Value.Coeffs[:auxLen+baseLen]
	e.auxLen = auxLen
}

// Copy returns a copy of e.
func (e *Element) Copy() *Element {
	return &Element{
		Value:  e.Value.Copy(),
		auxLen: e.auxLen,
	}
}

// CopyFrom copies the coefficients from eIn to e.
//
// Panics when e and eIn are not consistent.
func (e *Element) CopyFrom(eIn *Element) {
	e.Value.CopyFrom(eIn.Value)
	e.auxLen = eIn.auxLen
}

// IsEqual checks if e is equal to eIn.
func (e *Element) IsEqual(eIn *Element) bool {
	return e.Value.IsEqual(eIn.Value) && e.auxLen == eIn.auxLen
}

// IsConsistent checks if e has the same shape as eIn.
func (e *Element) IsConsistent(eIn *Element) bool {
	return e.Value.IsConsistent(eIn.Value)
}

// SecretKey is an RLWE secret key.
type SecretKey Element

// NewSecretKey creates a new [SecretKey].
func NewSecretKey(params Parameters, withAux bool, form dft.Form) *SecretKey {
	return (*SecretKey)(NewPoly(params, withAux, form))
}

// NewSecretKeyCustom creates a new [SecretKey] with the given parameters.
func NewSecretKeyCustom(rank, baseLen, auxLen int, form dft.Form) *SecretKey {
	return (*SecretKey)(NewPolyCustom(rank, baseLen, auxLen, form))
}

// Rank returns the rank of sk.
func (sk *SecretKey) Rank() int {
	return (*Element)(sk).Rank()
}

// Form returns the [dft.Form] of sk.
func (sk *SecretKey) Form() dft.Form {
	return (*Element)(sk).Form()
}

// SetForm sets the [dft.Form] of sk.
// This does not change the coefficients, so it may invalidate correctness.
func (sk *SecretKey) SetForm(form dft.Form) {
	(*Element)(sk).SetForm(form)
}

// FullModLen returns the full modulus length.
func (sk *SecretKey) FullModLen() int {
	return (*Element)(sk).FullModLen()
}

// BaseModLen returns the base modulus length.
func (sk *SecretKey) BaseModLen() int {
	return (*Element)(sk).BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (sk *SecretKey) AuxModLen() int {
	return (*Element)(sk).AuxModLen()
}

// Clear clears sk.
func (sk *SecretKey) Clear() {
	(*Element)(sk).Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (sk *SecretKey) WithModLen(baseLen, auxLen int) *SecretKey {
	return (*SecretKey)((*Element)(sk).WithModLen(baseLen, auxLen))
}

// Copy returns a copy of sk.
func (sk *SecretKey) Copy() *SecretKey {
	return (*SecretKey)((*Element)(sk).Copy())
}

// CopyFrom copies the coefficients from skIn to sk.
//
// Panics when sk and skIn are not consistent.
func (sk *SecretKey) CopyFrom(skIn *SecretKey) {
	(*Element)(sk).CopyFrom((*Element)(skIn))
}

// IsEqual checks if sk is equal to skIn.
func (sk *SecretKey) IsEqual(skIn *SecretKey) bool {
	return (*Element)(sk).IsEqual((*Element)(skIn))
}

// IsConsistent checks if sk has the same shape as skIn.
func (sk *SecretKey) IsConsistent(skIn *SecretKey) bool {
	return (*Element)(sk).IsConsistent((*Element)(skIn))
}

// Ciphertext is an RLWE ciphertext.
type Ciphertext struct {
	Body, Mask *Element
}

// NewCiphertext creates a new [Ciphertext].
func NewCiphertext(params Parameters, withAux bool, form dft.Form) *Ciphertext {
	return &Ciphertext{
		Body: NewPoly(params, withAux, form),
		Mask: NewPoly(params, withAux, form),
	}
}

// NewCiphertextCustom creates a new [Ciphertext] with the given parameters.
func NewCiphertextCustom(rank, baseLen, auxLen int, form dft.Form) *Ciphertext {
	return &Ciphertext{
		Body: NewPolyCustom(rank, baseLen, auxLen, form),
		Mask: NewPolyCustom(rank, baseLen, auxLen, form),
	}
}

// Rank returns the rank of ct.
func (ct *Ciphertext) Rank() int {
	if ct.Body.Rank() != ct.Mask.Rank() {
		panic("inconsistent rank")
	}
	return ct.Body.Rank()
}

// Form returns the [dft.Form] of ct.
func (ct *Ciphertext) Form() dft.Form {
	if ct.Body.Form() != ct.Mask.Form() {
		panic("inconsistent forms")
	}
	return ct.Body.Form()
}

// SetForm sets the [dft.Form] of ct.
// This does not change the coefficients, so it may invalidate correctness.
func (ct *Ciphertext) SetForm(form dft.Form) {
	ct.Body.SetForm(form)
	ct.Mask.SetForm(form)
}

// FullModLen returns the full modulus length.
func (ct *Ciphertext) FullModLen() int {
	if ct.Body.FullModLen() != ct.Mask.FullModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.FullModLen()
}

// BaseModLen returns the base modulus length.
func (ct *Ciphertext) BaseModLen() int {
	if ct.Body.BaseModLen() != ct.Mask.BaseModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (ct *Ciphertext) AuxModLen() int {
	if ct.Body.AuxModLen() != ct.Mask.AuxModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.AuxModLen()
}

// Clear clears ct.
func (ct *Ciphertext) Clear() {
	ct.Body.Clear()
	ct.Mask.Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (ct *Ciphertext) WithModLen(baseLen, auxLen int) *Ciphertext {
	return &Ciphertext{
		Body: ct.Body.WithModLen(baseLen, auxLen),
		Mask: ct.Mask.WithModLen(baseLen, auxLen),
	}
}

// Resize resizes ct in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is negative.
func (ct *Ciphertext) Resize(baseLen, auxLen int) {
	ct.Body.Resize(baseLen, auxLen)
	ct.Mask.Resize(baseLen, auxLen)
}

// Copy returns a copy of ct.
func (ct *Ciphertext) Copy() *Ciphertext {
	return &Ciphertext{
		Body: ct.Body.Copy(),
		Mask: ct.Mask.Copy(),
	}
}

// CopyFrom copies the coefficients from ctIn to ct.
//
// Panics when ct and ctIn are not consistent.
func (ct *Ciphertext) CopyFrom(ctIn *Ciphertext) {
	if !ct.IsConsistent(ctIn) {
		panic("inconsistent inputs")
	}
	ct.Body.CopyFrom(ctIn.Body)
	ct.Mask.CopyFrom(ctIn.Mask)
}

// IsEqual checks if ct is equal to ctIn.
func (ct *Ciphertext) IsEqual(ctIn *Ciphertext) bool {
	return ct.Body.IsEqual(ctIn.Body) && ct.Mask.IsEqual(ctIn.Mask)
}

// IsConsistent checks if ct has the same shape as ctIn.
func (ct *Ciphertext) IsConsistent(ctIn *Ciphertext) bool {
	return ct.Body.IsConsistent(ctIn.Body) && ct.Mask.IsConsistent(ctIn.Mask)
}

// PublicKey is an RLWE public key.
type PublicKey Ciphertext

// NewPublicKey creates a new [PublicKey].
func NewPublicKey(params Parameters, withAux bool, form dft.Form) *PublicKey {
	return (*PublicKey)(NewCiphertext(params, withAux, form))
}

// NewPublicKeyCustom creates a new [PublicKey] with the given parameters.
func NewPublicKeyCustom(rank, baseLen, auxLen int, form dft.Form) *PublicKey {
	return (*PublicKey)(NewCiphertextCustom(rank, baseLen, auxLen, form))
}

// Rank returns the rank of pk.
func (pk *PublicKey) Rank() int {
	return (*Ciphertext)(pk).Rank()
}

// Form returns the [dft.Form] of pk.
func (pk *PublicKey) Form() dft.Form {
	return (*Ciphertext)(pk).Form()
}

// SetForm sets the [dft.Form] of pk.
// This does not change the coefficients, so it may invalidate correctness.
func (pk *PublicKey) SetForm(form dft.Form) {
	(*Ciphertext)(pk).SetForm(form)
}

// FullModLen returns the full modulus length.
func (pk *PublicKey) FullModLen() int {
	return (*Ciphertext)(pk).FullModLen()
}

// BaseModLen returns the base modulus length.
func (pk *PublicKey) BaseModLen() int {
	return (*Ciphertext)(pk).BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (pk *PublicKey) AuxModLen() int {
	return (*Ciphertext)(pk).AuxModLen()
}

// Clear clears pk.
func (pk *PublicKey) Clear() {
	(*Ciphertext)(pk).Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (pk *PublicKey) WithModLen(baseLen, auxLen int) *PublicKey {
	return (*PublicKey)((*Ciphertext)(pk).WithModLen(baseLen, auxLen))
}

// Copy returns a copy of pk.
func (pk *PublicKey) Copy() *PublicKey {
	return (*PublicKey)((*Ciphertext)(pk).Copy())
}

// CopyFrom copies the coefficients from pkIn to pk.
//
// Panics when pk and pkIn are not consistent.
func (pk *PublicKey) CopyFrom(pkIn *PublicKey) {
	(*Ciphertext)(pk).CopyFrom((*Ciphertext)(pkIn))
}

// IsEqual checks if pk is equal to pkIn.
func (pk *PublicKey) IsEqual(pkIn *PublicKey) bool {
	return (*Ciphertext)(pk).IsEqual((*Ciphertext)(pkIn))
}

// IsConsistent checks if pk has the same shape as pkIn.
func (pk *PublicKey) IsConsistent(pkIn *PublicKey) bool {
	return (*Ciphertext)(pk).IsConsistent((*Ciphertext)(pkIn))
}

// Vector is a vector of [Element].
type Vector struct {
	Value []*Element
}

// NewVector creates a new [Vector].
func NewVector(params Parameters, length int, withAux bool, elType crt.ElementType, form dft.Form) *Vector {
	if withAux && params.AuxModLen() == 0 {
		panic("no auxiliary modulus")
	}

	if withAux {
		return NewVectorCustom(params.Rank(), params.BaseModLen(), params.AuxModLen(), length, elType, form)
	}
	return NewVectorCustom(params.Rank(), params.BaseModLen(), 0, length, elType, form)
}

// NewVectorCustom creates a new [Vector] with the given parameters.
func NewVectorCustom(rank, baseLen, auxLen, length int, elType crt.ElementType, form dft.Form) *Vector {
	value := make([]*Element, length)
	for i := range value {
		value[i] = NewElementCustom(rank, baseLen, auxLen, elType, form)
	}

	return &Vector{Value: value}
}

// Len returns the length of v.
func (v *Vector) Len() int {
	return len(v.Value)
}

// Rank returns the rank of v.
func (v *Vector) Rank() int {
	rank := v.Value[0].Rank()
	for i := 1; i < v.Len(); i++ {
		if v.Value[i].Rank() != rank {
			panic("inconsistent rank")
		}
	}
	return rank
}

// Form returns the [dft.Form] of v.
func (v *Vector) Form() dft.Form {
	form := v.Value[0].Form()
	for i := 1; i < v.Len(); i++ {
		if form != v.Value[i].Form() {
			panic("inconsistent forms")
		}
	}
	return form
}

// SetForm sets the [dft.Form] of v.
// This does not change the coefficients, so it may invalidate correctness.
func (v *Vector) SetForm(form dft.Form) {
	for i := range v.Value {
		v.Value[i].SetForm(form)
	}
}

// Type returns the [crt.ElementType] of v.
func (v *Vector) Type() crt.ElementType {
	elType := v.Value[0].Type()
	for i := 1; i < v.Len(); i++ {
		if elType != v.Value[i].Type() {
			panic("inconsistent types")
		}
	}
	return elType
}

// FullModLen returns the full modulus length.
func (v *Vector) FullModLen() int {
	fullLen := v.Value[0].FullModLen()
	for i := 1; i < v.Len(); i++ {
		if fullLen != v.Value[i].FullModLen() {
			panic("inconsistent modulus length")
		}
	}
	return fullLen
}

// BaseModLen returns the base modulus length.
func (v *Vector) BaseModLen() int {
	baseLen := v.Value[0].BaseModLen()
	for i := 1; i < v.Len(); i++ {
		if baseLen != v.Value[i].BaseModLen() {
			panic("inconsistent modulus length")
		}
	}
	return baseLen
}

// AuxModLen returns the auxiliary modulus length.
func (v *Vector) AuxModLen() int {
	auxLen := v.Value[0].AuxModLen()
	for i := 1; i < v.Len(); i++ {
		if auxLen != v.Value[i].AuxModLen() {
			panic("inconsistent modulus length")
		}
	}
	return auxLen
}

// Clear clears v.
func (v *Vector) Clear() {
	for i := range v.Value {
		v.Value[i].Clear()
	}
}

// Slice returns a shallow copy with the given element indices.
func (v *Vector) Slice(lo, hi int) *Vector {
	return &Vector{Value: v.Value[lo:hi]}
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (v *Vector) WithModLen(baseLen, auxLen int) *Vector {
	if baseLen < 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	value := make([]*Element, v.Len())
	for i := range value {
		value[i] = v.Value[i].WithModLen(baseLen, auxLen)
	}

	return &Vector{Value: value}
}

// Resize resizes v in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is invalid.
func (v *Vector) Resize(baseLen, auxLen int) {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	for i := range v.Value {
		v.Value[i].Resize(baseLen, auxLen)
	}
}

// Copy returns a copy of v.
func (v *Vector) Copy() *Vector {
	value := make([]*Element, v.Len())
	for i := 0; i < v.Len(); i++ {
		value[i] = v.Value[i].Copy()
	}

	return &Vector{Value: value}
}

// CopyFrom copies the coefficients from vIn to v.
//
// Panics when v and vIn are not consistent.
func (v *Vector) CopyFrom(vIn *Vector) {
	if !v.IsConsistent(vIn) {
		panic("inconsistent inputs")
	}

	for i := 0; i < v.Len(); i++ {
		v.Value[i].CopyFrom(vIn.Value[i])
	}
}

// IsEqual checks if v is equal to vIn.
func (v *Vector) IsEqual(vIn *Vector) bool {
	if !v.IsConsistent(vIn) {
		return false
	}

	for i := 0; i < v.Len(); i++ {
		if !v.Value[i].IsEqual(vIn.Value[i]) {
			return false
		}
	}
	return true
}

// IsConsistent checks if v has the same shape as vIn.
func (v *Vector) IsConsistent(vIn *Vector) bool {
	if v.Len() != vIn.Len() {
		return false
	}

	for i := 0; i < v.Len(); i++ {
		if !v.Value[i].IsConsistent(vIn.Value[i]) {
			return false
		}
	}
	return true
}

// GadgetEncryption is a gadget encryption.
type GadgetEncryption struct {
	Value []*Ciphertext
}

// NewGadgetEncryption creates a new [GadgetEncryption].
func NewGadgetEncryption(params Parameters, form dft.Form) *GadgetEncryption {
	return NewGadgetEncryptionCustom(params.Rank(), params.BaseModLen(), params.AuxModLen(), params.GadgetLen(), form)
}

// NewGadgetEncryptionCustom creates a new [GadgetEncryption] with the given parameters.
func NewGadgetEncryptionCustom(rank, baseLen, auxLen, gadLen int, form dft.Form) *GadgetEncryption {
	value := make([]*Ciphertext, gadLen)
	for i := range value {
		value[i] = NewCiphertextCustom(rank, baseLen, auxLen, form)
	}

	return &GadgetEncryption{Value: value}
}

// GadgetLen returns the length of the gadget decomposition.
func (ct *GadgetEncryption) GadgetLen() int {
	return len(ct.Value)
}

// Rank returns the rank of ct.
func (ct *GadgetEncryption) Rank() int {
	rank := ct.Value[0].Rank()
	for i := 1; i < ct.GadgetLen(); i++ {
		if ct.Value[i].Rank() != rank {
			panic("inconsistent rank")
		}
	}
	return rank
}

// Form returns the [dft.Form] of ct.
func (ct *GadgetEncryption) Form() dft.Form {
	form := ct.Value[0].Form()
	for i := 1; i < ct.GadgetLen(); i++ {
		if form != ct.Value[i].Form() {
			panic("inconsistent forms")
		}
	}
	return form
}

// SetForm sets the [dft.Form] of ct.
// This does not change the coefficients, so it may invalidate correctness.
func (ct *GadgetEncryption) SetForm(form dft.Form) {
	for i := range ct.Value {
		ct.Value[i].SetForm(form)
	}
}

// FullModLen returns the full modulus length.
func (ct *GadgetEncryption) FullModLen() int {
	fullLen := ct.Value[0].FullModLen()
	for i := 1; i < ct.GadgetLen(); i++ {
		if fullLen != ct.Value[i].FullModLen() {
			panic("inconsistent modulus length")
		}
	}
	return fullLen
}

// BaseModLen returns the base modulus length.
func (ct *GadgetEncryption) BaseModLen() int {
	baseLen := ct.Value[0].BaseModLen()
	for i := 1; i < ct.GadgetLen(); i++ {
		if baseLen != ct.Value[i].BaseModLen() {
			panic("inconsistent modulus length")
		}
	}
	return baseLen
}

// AuxModLen returns the auxiliary modulus length.
func (ct *GadgetEncryption) AuxModLen() int {
	auxLen := ct.Value[0].AuxModLen()
	for i := 1; i < ct.GadgetLen(); i++ {
		if auxLen != ct.Value[i].AuxModLen() {
			panic("inconsistent modulus length")
		}
	}
	return auxLen
}

// Clear clears ct.
func (ct *GadgetEncryption) Clear() {
	for i := range ct.Value {
		ct.Value[i].Clear()
	}
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (ct *GadgetEncryption) WithModLen(baseLen, auxLen int) *GadgetEncryption {
	if baseLen < 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	value := make([]*Ciphertext, ct.GadgetLen())
	for i := range value {
		value[i] = ct.Value[i].WithModLen(baseLen, auxLen)
	}

	return &GadgetEncryption{Value: value}
}

// Resize resizes ct in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is negative.
func (ct *GadgetEncryption) Resize(baseLen, auxLen int) {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	for i := range ct.Value {
		ct.Value[i].Resize(baseLen, auxLen)
	}
}

// Copy returns a copy of ct.
func (ct *GadgetEncryption) Copy() *GadgetEncryption {
	value := make([]*Ciphertext, ct.GadgetLen())
	for i := 0; i < ct.GadgetLen(); i++ {
		value[i] = ct.Value[i].Copy()
	}

	return &GadgetEncryption{Value: value}
}

// CopyFrom copies the coefficients from ctIn to ct.
//
// Panics when ct and ctIn are not consistent.
func (ct *GadgetEncryption) CopyFrom(ctIn *GadgetEncryption) {
	if !ct.IsConsistent(ctIn) {
		panic("inconsistent inputs")
	}

	for i := 0; i < ct.GadgetLen(); i++ {
		ct.Value[i].CopyFrom(ctIn.Value[i])
	}
}

// IsEqual checks if ct is equal to ctIn.
func (ct *GadgetEncryption) IsEqual(ctIn *GadgetEncryption) bool {
	if !ct.IsConsistent(ctIn) {
		return false
	}

	for i := 0; i < ct.GadgetLen(); i++ {
		if !ct.Value[i].IsEqual(ctIn.Value[i]) {
			return false
		}
	}
	return true
}

// IsConsistent checks if ct has the same shape as ctIn.
func (ct *GadgetEncryption) IsConsistent(ctIn *GadgetEncryption) bool {
	if ct.GadgetLen() != ctIn.GadgetLen() {
		return false
	}

	for i := 0; i < ct.GadgetLen(); i++ {
		if !ct.Value[i].IsConsistent(ctIn.Value[i]) {
			return false
		}
	}
	return true
}

// RGSWEncryption is a Ring-GSW encryption.
type RGSWEncryption struct {
	Body, Mask *GadgetEncryption
}

// NewRGSWEncryption creates a new [RGSWEncryption].
func NewRGSWEncryption(params Parameters, form dft.Form) *RGSWEncryption {
	return &RGSWEncryption{
		Body: NewGadgetEncryption(params, form),
		Mask: NewGadgetEncryption(params, form),
	}
}

// NewRGSWEncryptionCustom creates a new [RGSWEncryption] with the given parameters.
func NewRGSWEncryptionCustom(rank, baseLen, auxLen, gadLen int, form dft.Form) *RGSWEncryption {
	return &RGSWEncryption{
		Body: NewGadgetEncryptionCustom(rank, baseLen, auxLen, gadLen, form),
		Mask: NewGadgetEncryptionCustom(rank, baseLen, auxLen, gadLen, form),
	}
}

// GadgetLen returns the length of the gadget decomposition.
func (ct *RGSWEncryption) GadgetLen() int {
	if ct.Body.GadgetLen() != ct.Mask.GadgetLen() {
		panic("inconsistent gadget length")
	}
	return ct.Body.GadgetLen()
}

// Rank returns the rank of ct.
func (ct *RGSWEncryption) Rank() int {
	if ct.Body.Rank() != ct.Mask.Rank() {
		panic("inconsistent rank")
	}
	return ct.Body.Rank()
}

// Form returns the [dft.Form] of ct.
func (ct *RGSWEncryption) Form() dft.Form {
	if ct.Body.Form() != ct.Mask.Form() {
		panic("inconsistent forms")
	}
	return ct.Body.Form()
}

// SetForm sets the [dft.Form] of ct.
// This does not change the coefficients, so it may invalidate correctness.
func (ct *RGSWEncryption) SetForm(form dft.Form) {
	ct.Body.SetForm(form)
	ct.Mask.SetForm(form)
}

// FullModLen returns the full modulus length.
func (ct *RGSWEncryption) FullModLen() int {
	if ct.Body.FullModLen() != ct.Mask.FullModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.FullModLen()
}

// BaseModLen returns the base modulus length.
func (ct *RGSWEncryption) BaseModLen() int {
	if ct.Body.BaseModLen() != ct.Mask.BaseModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (ct *RGSWEncryption) AuxModLen() int {
	if ct.Body.AuxModLen() != ct.Mask.AuxModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.AuxModLen()
}

// Clear clears ct.
func (ct *RGSWEncryption) Clear() {
	ct.Body.Clear()
	ct.Mask.Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (ct *RGSWEncryption) WithModLen(baseLen, auxLen int) *RGSWEncryption {
	return &RGSWEncryption{
		Body: ct.Body.WithModLen(baseLen, auxLen),
		Mask: ct.Mask.WithModLen(baseLen, auxLen),
	}
}

// Resize resizes ct in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is negative.
func (ct *RGSWEncryption) Resize(baseLen, auxLen int) {
	ct.Body.Resize(baseLen, auxLen)
	ct.Mask.Resize(baseLen, auxLen)
}

// Copy returns a copy of ct.
func (ct *RGSWEncryption) Copy() *RGSWEncryption {
	return &RGSWEncryption{
		Body: ct.Body.Copy(),
		Mask: ct.Mask.Copy(),
	}
}

// CopyFrom copies the coefficients from ctIn to ct.
//
// Panics when ct and ctIn are not consistent.
func (ct *RGSWEncryption) CopyFrom(ctIn *RGSWEncryption) {
	if !ct.IsConsistent(ctIn) {
		panic("inconsistent inputs")
	}

	ct.Body.CopyFrom(ctIn.Body)
	ct.Mask.CopyFrom(ctIn.Mask)
}

// IsEqual checks if ct is equal to ctIn.
func (ct *RGSWEncryption) IsEqual(ctIn *RGSWEncryption) bool {
	return ct.Body.IsEqual(ctIn.Body) && ct.Mask.IsEqual(ctIn.Mask)
}

// IsConsistent checks if ct has the same shape as ctIn.
func (ct *RGSWEncryption) IsConsistent(ctIn *RGSWEncryption) bool {
	return ct.Body.IsConsistent(ctIn.Body) && ct.Mask.IsConsistent(ctIn.Mask)
}

// KeySwitchKey is a public key for key-switching.
type KeySwitchKey GadgetEncryption

// NewKeySwitchKey creates a new [KeySwitchKey].
func NewKeySwitchKey(params Parameters, form dft.Form) *KeySwitchKey {
	return (*KeySwitchKey)(NewGadgetEncryption(params, form))
}

// NewKeySwitchKeyCustom creates a new [KeySwitchKey] with the given parameters.
func NewKeySwitchKeyCustom(rank, baseLen, auxLen, gadLen int, form dft.Form) *KeySwitchKey {
	return (*KeySwitchKey)(NewGadgetEncryptionCustom(rank, baseLen, auxLen, gadLen, form))
}

// GadgetLen returns the length of the gadget decomposition.
func (ksk *KeySwitchKey) GadgetLen() int {
	return (*GadgetEncryption)(ksk).GadgetLen()
}

// Rank returns the rank of ksk.
func (ksk *KeySwitchKey) Rank() int {
	return (*GadgetEncryption)(ksk).Rank()
}

// Form returns the [dft.Form] of ksk.
func (ksk *KeySwitchKey) Form() dft.Form {
	return (*GadgetEncryption)(ksk).Form()
}

// SetForm sets the [dft.Form] of ksk.
// This does not change the coefficients, so it may invalidate correctness.
func (ksk *KeySwitchKey) SetForm(form dft.Form) {
	(*GadgetEncryption)(ksk).SetForm(form)
}

// FullModLen returns the full modulus length.
func (ksk *KeySwitchKey) FullModLen() int {
	return (*GadgetEncryption)(ksk).FullModLen()
}

// BaseModLen returns the base modulus length.
func (ksk *KeySwitchKey) BaseModLen() int {
	return (*GadgetEncryption)(ksk).BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (ksk *KeySwitchKey) AuxModLen() int {
	return (*GadgetEncryption)(ksk).AuxModLen()
}

// Clear clears ksk.
func (ksk *KeySwitchKey) Clear() {
	(*GadgetEncryption)(ksk).Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (ksk *KeySwitchKey) WithModLen(baseLen, auxLen int) *KeySwitchKey {
	return (*KeySwitchKey)((*GadgetEncryption)(ksk).WithModLen(baseLen, auxLen))
}

// Resize resizes ksk in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is invalid.
func (ksk *KeySwitchKey) Resize(baseLen, auxLen int) {
	(*GadgetEncryption)(ksk).Resize(baseLen, auxLen)
}

// Copy returns a copy of ksk.
func (ksk *KeySwitchKey) Copy() *KeySwitchKey {
	return (*KeySwitchKey)((*GadgetEncryption)(ksk).Copy())
}

// CopyFrom copies the coefficients from kskIn to ksk.
//
// Panics when ksk and kskIn are not consistent.
func (ksk *KeySwitchKey) CopyFrom(kskIn *KeySwitchKey) {
	(*GadgetEncryption)(ksk).CopyFrom((*GadgetEncryption)(kskIn))
}

// IsEqual checks if ksk is equal to kskIn.
func (ksk *KeySwitchKey) IsEqual(kskIn *KeySwitchKey) bool {
	return (*GadgetEncryption)(ksk).IsEqual((*GadgetEncryption)(kskIn))
}

// IsConsistent checks if ksk has the same shape as kskIn.
func (ksk *KeySwitchKey) IsConsistent(kskIn *KeySwitchKey) bool {
	return (*GadgetEncryption)(ksk).IsConsistent((*GadgetEncryption)(kskIn))
}

// RelinKey is a public key for relinearisation.
type RelinKey KeySwitchKey

// NewRelinKey creates a new [RelinKey].
func NewRelinKey(params Parameters, form dft.Form) *RelinKey {
	return (*RelinKey)(NewKeySwitchKey(params, form))
}

// NewRelinKeyCustom creates a new [RelinKey] with the given parameters.
func NewRelinKeyCustom(rank, baseLen, auxLen, gadLen int, form dft.Form) *RelinKey {
	return (*RelinKey)(NewKeySwitchKeyCustom(rank, baseLen, auxLen, gadLen, form))
}

// GadgetLen returns the length of the gadget decomposition.
func (rlk *RelinKey) GadgetLen() int {
	return (*KeySwitchKey)(rlk).GadgetLen()
}

// Rank returns the rank of rlk.
func (rlk *RelinKey) Rank() int {
	return (*KeySwitchKey)(rlk).Rank()
}

// Form returns the [dft.Form] of rlk.
func (rlk *RelinKey) Form() dft.Form {
	return (*KeySwitchKey)(rlk).Form()
}

// SetForm sets the [dft.Form] of rlk.
// This does not change the coefficients, so it may invalidate correctness.
func (rlk *RelinKey) SetForm(form dft.Form) {
	(*KeySwitchKey)(rlk).SetForm(form)
}

// FullModLen returns the full modulus length.
func (rlk *RelinKey) FullModLen() int {
	return (*KeySwitchKey)(rlk).FullModLen()
}

// BaseModLen returns the base modulus length.
func (rlk *RelinKey) BaseModLen() int {
	return (*KeySwitchKey)(rlk).BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (rlk *RelinKey) AuxModLen() int {
	return (*KeySwitchKey)(rlk).AuxModLen()
}

// Clear clears rlk.
func (rlk *RelinKey) Clear() {
	(*KeySwitchKey)(rlk).Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (rlk *RelinKey) WithModLen(baseLen, auxLen int) *RelinKey {
	return (*RelinKey)((*KeySwitchKey)(rlk).WithModLen(baseLen, auxLen))
}

// Resize resizes rlk in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is invalid.
func (rlk *RelinKey) Resize(baseLen, auxLen int) {
	(*KeySwitchKey)(rlk).Resize(baseLen, auxLen)
}

// Copy returns a copy of rlk.
func (rlk *RelinKey) Copy() *RelinKey {
	return (*RelinKey)((*KeySwitchKey)(rlk).Copy())
}

// CopyFrom copies the coefficients from rlkIn to rlk.
//
// Panics when rlk and rlkIn are not consistent.
func (rlk *RelinKey) CopyFrom(rlkIn *RelinKey) {
	(*KeySwitchKey)(rlk).CopyFrom((*KeySwitchKey)(rlkIn))
}

// IsEqual checks if rlk is equal to rlkIn.
func (rlk *RelinKey) IsEqual(rlkIn *RelinKey) bool {
	return (*KeySwitchKey)(rlk).IsEqual((*KeySwitchKey)(rlkIn))
}

// IsConsistent checks if rlk has the same shape as rlkIn.
func (rlk *RelinKey) IsConsistent(rlkIn *RelinKey) bool {
	return (*KeySwitchKey)(rlk).IsConsistent((*KeySwitchKey)(rlkIn))
}

// AutKey is a public key for automorphism.
type AutKey struct {
	Value *KeySwitchKey
	// Idx is the automorphism index.
	Idx int
}

// NewAutKey creates a new [AutKey] with the given automorphism index.
func NewAutKey(params Parameters, idx int, form dft.Form) *AutKey {
	return &AutKey{
		Value: NewKeySwitchKey(params, form),
		Idx:   idx,
	}
}

// NewAutKeyCustom creates a new [AutKey] with the given parameters.
func NewAutKeyCustom(rank, baseLen, auxLen, gadLen, idx int, form dft.Form) *AutKey {
	return &AutKey{
		Value: NewKeySwitchKeyCustom(rank, baseLen, auxLen, gadLen, form),
		Idx:   idx,
	}
}

// GadgetLen returns the length of the gadget decomposition.
func (atk *AutKey) GadgetLen() int {
	return atk.Value.GadgetLen()
}

// Rank returns the rank of atk.
func (atk *AutKey) Rank() int {
	return atk.Value.Rank()
}

// Form returns the [dft.Form] of atk.
func (atk *AutKey) Form() dft.Form {
	return atk.Value.Form()
}

// SetForm sets the [dft.Form] of atk.
// This does not change the coefficients, so it may invalidate correctness.
func (atk *AutKey) SetForm(form dft.Form) {
	atk.Value.SetForm(form)
}

// FullModLen returns the full modulus length.
func (atk *AutKey) FullModLen() int {
	return atk.Value.FullModLen()
}

// BaseModLen returns the base modulus length.
func (atk *AutKey) BaseModLen() int {
	return atk.Value.BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (atk *AutKey) AuxModLen() int {
	return atk.Value.AuxModLen()
}

// Clear clears the coefficients of atk, preserving its automorphism index.
func (atk *AutKey) Clear() {
	atk.Value.Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (atk *AutKey) WithModLen(baseLen, auxLen int) *AutKey {
	return &AutKey{
		Value: atk.Value.WithModLen(baseLen, auxLen),
		Idx:   atk.Idx,
	}
}

// Resize resizes atk in place to the given base and auxiliary modulus lengths.
//
// Panics when baseLen or auxLen is invalid.
func (atk *AutKey) Resize(baseLen, auxLen int) {
	atk.Value.Resize(baseLen, auxLen)
}

// Copy returns a copy of atk.
func (atk *AutKey) Copy() *AutKey {
	return &AutKey{
		Value: atk.Value.Copy(),
		Idx:   atk.Idx,
	}
}

// CopyFrom copies the coefficients and automorphism index from atkIn to atk.
//
// Panics when atk and atkIn are not consistent.
func (atk *AutKey) CopyFrom(atkIn *AutKey) {
	atk.Value.CopyFrom(atkIn.Value)
	atk.Idx = atkIn.Idx
}

// IsEqual checks if atk is equal to atkIn, including its automorphism index.
func (atk *AutKey) IsEqual(atkIn *AutKey) bool {
	return atk.Value.IsEqual(atkIn.Value) && atk.Idx == atkIn.Idx
}

// IsConsistent checks if atk has the same shape as atkIn.
func (atk *AutKey) IsConsistent(atkIn *AutKey) bool {
	return atk.Value.IsConsistent(atkIn.Value)
}
