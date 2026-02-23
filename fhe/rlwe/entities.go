package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
)

type Element struct {
	Value *crt.Element
	// hasAux indicates whether Value has an auxiliary modulus.
	// TODO: We should consider saving AuxModLen instead of HasAux.
	hasAux bool
}

// TODO: By adding AuxLen, we should add/change the following:
// - All functions currently accepting (modLen, hasAux) should accept (baseLen, auxLen)
// - In `Copy`, `CopyFrom`, additional checks for AuxLen should be made
// - We should remove `ModLen`, and add `BaseModLen`, `AuxModLen` and `FullModLen`.

// NewScalar creates a new scalar [Element].
func NewScalar(params Parameters, hasAux bool) *Element {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	if !hasAux {
		return NewScalarCustom(len(params.baseMod), hasAux)
	}
	return NewScalarCustom(len(params.fullMod), hasAux)
}

// NewScalarCustom creates a new scalar [Element] with the given parameters.
func NewScalarCustom(modLen int, hasAux bool) *Element {
	return &Element{
		Value:  crt.NewScalar(modLen),
		hasAux: hasAux,
	}
}

// NewPoly creates a new polynomial [Element].
func NewPoly(params Parameters, hasAux, isNTT bool) *Element {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	if !hasAux {
		return NewPolyCustom(params.RingParams().Rank(), len(params.baseMod), hasAux, isNTT)
	}
	return NewPolyCustom(params.RingParams().Rank(), len(params.fullMod), hasAux, isNTT)
}

// NewPolyCustom creates a new polynomial [Element] with the given parameters.
func NewPolyCustom(rank, modLen int, hasAux, isNTT bool) *Element {
	return &Element{
		Value:  crt.NewPolyCustom(rank, modLen, isNTT),
		hasAux: hasAux,
	}
}

// Rank returns the rank.
func (e *Element) Rank() int {
	return e.Value.Rank()
}

// ModLen returns the modulus length.
func (e *Element) ModLen() int {
	return e.Value.ModLen()
}

// IsNTT returns whether value is in NTT form.
func (e *Element) IsNTT() bool {
	return e.Value.IsNTT
}

// HasAuxModulus returns whether auxiliary modulus is used.
func (e *Element) HasAuxModulus() bool {
	return e.hasAux
}

// Clear clears value.
func (e *Element) Clear() {
	e.Value.Clear()
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
//
// TODO: This is not robust, as `idx` also spans auxillary modulus,
// thus we don't know whether hasAux flag is valid.
func (e *Element) WithModIdx(idx ...int) *Element {
	return &Element{
		Value:  e.Value.WithModIdx(idx...),
		hasAux: e.hasAux,
	}
}

// Copy returns a copy.
func (e *Element) Copy() *Element {
	return &Element{
		Value:  e.Value.Copy(),
		hasAux: e.hasAux,
	}
}

// CopyFrom copies the given value.
//
// Panics when not consistent.
func (e *Element) CopyFrom(eIn *Element) {
	if !e.IsConsistent(eIn) {
		panic("inconsistent input(s)")
	}
	e.Value.CopyFrom(eIn.Value)
	e.hasAux = eIn.hasAux
}

// IsEqual checks if two values are equal.
func (e *Element) IsEqual(e0 *Element) bool {
	if !e.IsConsistent(e0) {
		return false
	}
	return e.Value.IsEqual(e0.Value) && e.hasAux == e0.hasAux
}

// IsConsistent checks if two values have the same shape.
func (e *Element) IsConsistent(e0 *Element) bool {
	return e.Value.IsConsistent(e0.Value)
}

// SecretKey is an RLWE secret key.
type SecretKey Element

// NewSecretKey creates a new [SecretKey].
func NewSecretKey(params Parameters, hasAux bool) *SecretKey {
	return (*SecretKey)(NewPoly(params, hasAux, false))
}

// NewSecretKeyCustom creates a new [SecretKey] with the given parameters.
func NewSecretKeyCustom(rank, modLen int, hasAux, isNTT bool) *SecretKey {
	return (*SecretKey)(NewPolyCustom(rank, modLen, hasAux, isNTT))
}

// Rank returns the rank.
func (sk *SecretKey) Rank() int {
	return (*Element)(sk).Rank()
}

// ModLen returns the modulus length.
func (sk *SecretKey) ModLen() int {
	return (*Element)(sk).ModLen()
}

// IsNTT returns whether value is in NTT form.
func (sk *SecretKey) IsNTT() bool {
	return (*Element)(sk).IsNTT()
}

// HasAuxModulus returns whether auxiliary modulus is used.
func (sk *SecretKey) HasAuxModulus() bool {
	return (*Element)(sk).HasAuxModulus()
}

// Clear clears value.
func (sk *SecretKey) Clear() {
	(*Element)(sk).Clear()
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
func (sk *SecretKey) WithModIdx(idx ...int) *SecretKey {
	return (*SecretKey)((*Element)(sk).WithModIdx(idx...))
}

// Copy returns a copy.
func (sk *SecretKey) Copy() *SecretKey {
	return (*SecretKey)((*Element)(sk).Copy())
}

// CopyFrom copies the given value.
//
// Panics when not consistent.
func (sk *SecretKey) CopyFrom(skIn *SecretKey) {
	(*Element)(sk).CopyFrom((*Element)(skIn))
}

// IsEqual checks if two values are equal.
func (sk *SecretKey) IsEqual(e0 *Element) bool {
	return (*Element)(sk).IsEqual(e0)
}

// IsConsistent checks if two values have the same shape.
func (sk *SecretKey) IsConsistent(e0 *Element) bool {
	return (*Element)(sk).IsConsistent(e0)
}

// Ciphertext is an RLWE ciphertext.
type Ciphertext struct {
	Body, Mask *Element
}

// NewCiphertext creates a new [Ciphertext].
func NewCiphertext(params Parameters, hasAux bool, isNTT bool) *Ciphertext {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	if !hasAux {
		return NewCiphertextCustom(params.RingParams().Rank(), len(params.baseMod), hasAux, isNTT)
	}
	return NewCiphertextCustom(params.RingParams().Rank(), len(params.fullMod), hasAux, isNTT)
}

// NewCiphertextCustom creates a new [Ciphertext] with the given parameters.
func NewCiphertextCustom(rank, modLen int, hasAux, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Body: NewPolyCustom(rank, modLen, hasAux, isNTT),
		Mask: NewPolyCustom(rank, modLen, hasAux, isNTT),
	}
}

// Rank returns the rank.
func (ct *Ciphertext) Rank() int {
	if ct.Body.Rank() != ct.Mask.Rank() {
		panic("inconsistent rank")
	}
	return ct.Body.Rank()
}

// ModLen returns the modulus length.
func (ct *Ciphertext) ModLen() int {
	if ct.Body.ModLen() != ct.Mask.ModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.ModLen()
}

// IsNTT returns whether value is in NTT form.
func (ct *Ciphertext) IsNTT() bool {
	if ct.Body.IsNTT() != ct.Mask.IsNTT() {
		panic("inconsistent NTT form")
	}
	return ct.Body.IsNTT()
}

// HasAuxModulus returns whether ct has an auxiliary modulus.
func (ct *Ciphertext) HasAuxModulus() bool {
	if ct.Body.HasAuxModulus() != ct.Mask.HasAuxModulus() {
		panic("inconsistent auxiliary modulus")
	}
	return ct.Body.HasAuxModulus()
}

// Clear clears value.
func (ct *Ciphertext) Clear() {
	ct.Body.Clear()
	ct.Mask.Clear()
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
func (ct *Ciphertext) WithModIdx(idx ...int) *Ciphertext {
	return &Ciphertext{
		Body: ct.Body.WithModIdx(idx...),
		Mask: ct.Mask.WithModIdx(idx...),
	}
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
		panic("inconsistent input(s)")
	}
	ct.Body.CopyFrom(ctIn.Body)
	ct.Mask.CopyFrom(ctIn.Mask)
}

// IsEqual checks if two values are equal.
func (ct *Ciphertext) IsEqual(ct0 *Ciphertext) bool {
	return ct.Body.IsEqual(ct0.Body) && ct.Mask.IsEqual(ct0.Mask)
}

// IsConsistent checks if two values have the same shape.
func (ct *Ciphertext) IsConsistent(ct0 *Ciphertext) bool {
	return ct.Body.IsConsistent(ct0.Body) && ct.Mask.IsConsistent(ct0.Mask)
}

// PublicKey is an RLWE public key.
type PublicKey Ciphertext

// NewPublicKey creates a new [PublicKey].
func NewPublicKey(params Parameters, hasAux bool, isNTT bool) *PublicKey {
	return (*PublicKey)(NewCiphertext(params, hasAux, isNTT))
}

// NewPublicKeyCustom creates a new [PublicKey] with the given parameters.
func NewPublicKeyCustom(rank, modLen int, hasAux, isNTT bool) *PublicKey {
	return (*PublicKey)(NewCiphertextCustom(rank, modLen, hasAux, isNTT))
}

// Rank returns the rank.
func (pk *PublicKey) Rank() int {
	return (*Ciphertext)(pk).Rank()
}

// ModLen returns the modulus length.
func (pk *PublicKey) ModLen() int {
	return (*Ciphertext)(pk).ModLen()
}

// HasAux returns whether auxiliary modulus is used.
func (pk *PublicKey) HasAux() bool {
	return (*Ciphertext)(pk).HasAuxModulus()
}

// Clear clears value.
func (pk *PublicKey) Clear() {
	(*Ciphertext)(pk).Clear()
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
func (pk *PublicKey) WithModIdx(idx ...int) *PublicKey {
	return (*PublicKey)((*Ciphertext)(pk).WithModIdx(idx...))

}

// Copy returns a copy.
func (pk *PublicKey) Copy() *PublicKey {
	return (*PublicKey)((*Ciphertext)(pk).Copy())
}

// CopyFrom copies the given value.
//
// Panics when not consistent.
func (pk *PublicKey) CopyFrom(pkIn *PublicKey) {
	(*Ciphertext)(pk).CopyFrom((*Ciphertext)(pkIn))

}

// IsEqual checks if two values are equal.
func (pk *PublicKey) IsEqual(pk0 *PublicKey) bool {
	return (*Ciphertext)(pk).IsEqual((*Ciphertext)(pk0))
}

// IsConsistent checks if two values have the same shape.
func (pk *PublicKey) IsConsistent(pk0 *PublicKey) bool {
	return (*Ciphertext)(pk).IsConsistent((*Ciphertext)(pk0))
}

// Vector is a vector of RLWE elements.
type Vector struct {
	Value []*Element
}

// NewVector creates a new [Vector].
func NewVector(params Parameters, length int, hasAux, isNTT bool) *Vector {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	if !hasAux {
		return NewVectorCustom(params.RingParams().Rank(), len(params.baseMod), hasAux, length, isNTT)
	}
	return NewVectorCustom(params.RingParams().Rank(), len(params.fullMod), hasAux, length, isNTT)
}

// NewVectorCustom creates a new [Vector] with the given parameters.
func NewVectorCustom(rank, modLen int, hasAux bool, length int, isNTT bool) *Vector {
	value := make([]*Element, length)
	for i := range value {
		value[i] = NewPolyCustom(rank, modLen, hasAux, isNTT)
	}

	return &Vector{Value: value}
}

// Len returns the length.
func (v *Vector) Len() int {
	return len(v.Value)
}

// Rank returns the rank.
func (v *Vector) Rank() int {
	rank := v.Value[0].Rank()
	for i := 1; i < v.Len(); i++ {
		if v.Value[i].Rank() != rank {
			panic("inconsistent rank")
		}
	}
	return rank
}

// ModLen returns the modulus length of v.
func (v *Vector) ModLen() int {
	modLen := v.Value[0].ModLen()
	for i := 1; i < v.Len(); i++ {
		if modLen != v.Value[i].ModLen() {
			panic("inconsistent modulus length")
		}
	}
	return modLen
}

// IsNTT returns whether value is in NTT form.
func (v *Vector) IsNTT() bool {
	isNTT := v.Value[0].IsNTT()
	for i := 1; i < v.Len(); i++ {
		if v.Value[i].IsNTT() != isNTT {
			panic("inconsistent NTT form")
		}
	}
	return isNTT
}

// HasAuxModulus returns whether auxiliary modulus is used.
func (v *Vector) HasAuxModulus() bool {
	hasAux := v.Value[0].HasAuxModulus()
	for i := 1; i < v.Len(); i++ {
		if v.Value[i].HasAuxModulus() != hasAux {
			panic("inconsistent auxiliary modulus")
		}
	}
	return hasAux
}

// Clear clears value.
func (v *Vector) Clear() {
	for _, p := range v.Value {
		p.Clear()
	}
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
func (v *Vector) WithLenModIdx(length int, idx ...int) *Vector {
	if length < 1 || length > v.Len() {
		panic("length out of range")
	}

	value := make([]*Element, length)
	for i := 0; i < length; i++ {
		value[i] = v.Value[i].WithModIdx(idx...)
	}

	return &Vector{Value: value}
}

// Copy returns a copy.
func (v *Vector) Copy() *Vector {
	value := make([]*Element, v.Len())
	for i := 0; i < v.Len(); i++ {
		value[i] = v.Value[i].Copy()
	}

	return &Vector{Value: value}
}

// CopyFrom copies the given value.
//
// Panics when not consistent.
func (v *Vector) CopyFrom(vIn *Vector) {
	if !v.IsConsistent(vIn) {
		panic("inconsistent input(s)")
	}

	for i := 0; i < v.Len(); i++ {
		v.Value[i].CopyFrom(vIn.Value[i])
	}
}

// IsEqual checks if two values are equal.
func (v *Vector) IsEqual(v0 *Vector) bool {
	if !v.IsConsistent(v0) {
		return false
	}

	for i := 0; i < v.Len(); i++ {
		if !v.Value[i].IsEqual(v0.Value[i]) {
			return false
		}
	}
	return true
}

// IsConsistent checks if two values have the same shape.
func (v *Vector) IsConsistent(v0 *Vector) bool {
	if v.Len() != v0.Len() {
		return false
	}

	for i := 0; i < v.Len(); i++ {
		if !v.Value[i].IsConsistent(v0.Value[i]) {
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
func NewGadgetEncryption(params Parameters, isNTT bool) *GadgetEncryption {
	return NewGadgetEncryptionCustom(params.RingParams().Rank(), len(params.fullMod), params.HasAuxModulus(), params.GadgetLen(), isNTT)
}

// NewGadgetEncryptionCustom creates a new [GadgetEncryption] with the given parameters.
func NewGadgetEncryptionCustom(rank, modLen int, hasAux bool, gadLen int, isNTT bool) *GadgetEncryption {
	value := make([]*Ciphertext, gadLen)
	for i := range value {
		value[i] = NewCiphertextCustom(rank, modLen, hasAux, isNTT)
	}
	return &GadgetEncryption{Value: value}
}

// Rank returns the rank.
func (ct *GadgetEncryption) Rank() int {
	rank := ct.Value[0].Rank()
	for i := 1; i < len(ct.Value); i++ {
		if ct.Value[i].Rank() != rank {
			panic("inconsistent rank")
		}
	}
	return rank
}

// ModLen returns the modulus length.
func (ct *GadgetEncryption) ModLen() int {
	modLen := ct.Value[0].ModLen()
	for i := 1; i < len(ct.Value); i++ {
		if modLen != ct.Value[i].ModLen() {
			panic("inconsistent modulus length")
		}
	}
	return modLen
}

// IsNTT returns whether value is in NTT form.
func (ct *GadgetEncryption) IsNTT() bool {
	isNTT := ct.Value[0].IsNTT()
	for i := 1; i < len(ct.Value); i++ {
		if ct.Value[i].IsNTT() != isNTT {
			panic("inconsistent NTT form")
		}
	}
	return isNTT
}

// HasAuxModulus returns whether auxiliary modulus is used.
func (ct *GadgetEncryption) HasAuxModulus() bool {
	hasAux := ct.Value[0].HasAuxModulus()
	for i := 1; i < len(ct.Value); i++ {
		if ct.Value[i].HasAuxModulus() != hasAux {
			panic("inconsistent auxiliary modulus")
		}
	}
	return hasAux
}

// GadgetLen returns the length of the gadget.
func (ct *GadgetEncryption) GadgetLen() int {
	return len(ct.Value)
}

// Clear clears value.
func (ct *GadgetEncryption) Clear() {
	for i := range ct.Value {
		ct.Value[i].Clear()
	}
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
func (ct *GadgetEncryption) WithModIdx(idx ...int) *GadgetEncryption {
	value := make([]*Ciphertext, len(ct.Value))
	for i := 0; i < len(ct.Value); i++ {
		value[i] = ct.Value[i].WithModIdx(idx...)
	}
	return &GadgetEncryption{Value: value}
}

// Copy returns a copy.
func (ct *GadgetEncryption) Copy() *GadgetEncryption {
	value := make([]*Ciphertext, len(ct.Value))
	for i := 0; i < len(ct.Value); i++ {
		value[i] = ct.Value[i].Copy()
	}
	return &GadgetEncryption{Value: value}
}

// CopyFrom copies the given value.
//
// Panics when not consistent.
func (ct *GadgetEncryption) CopyFrom(ctIn *GadgetEncryption) {
	if !ct.IsConsistent(ctIn) {
		panic("inconsistent input(s)")
	}
	for i := 0; i < len(ct.Value); i++ {
		ct.Value[i].CopyFrom(ctIn.Value[i])
	}
}

// IsEqual checks if two values are equal.
func (ct *GadgetEncryption) IsEqual(ct0 *GadgetEncryption) bool {
	if !ct.IsConsistent(ct0) {
		return false
	}

	for i := 0; i < len(ct.Value); i++ {
		if !ct.Value[i].IsEqual(ct0.Value[i]) {
			return false
		}
	}
	return true
}

// IsConsistent checks if two values have the same shape.
func (ct *GadgetEncryption) IsConsistent(ct0 *GadgetEncryption) bool {
	if len(ct.Value) != len(ct0.Value) {
		return false
	}

	for i := 0; i < len(ct.Value); i++ {
		if !ct.Value[i].IsConsistent(ct0.Value[i]) {
			return false
		}
	}
	return true
}

// RGSW is a Ring-GSW ciphertext.
type RGSW struct {
	// TODO: Is this a good name?
	Body, Mask *GadgetEncryption
}

// NewRGSW creates a new [RGSW].
func NewRGSW(params Parameters, isNTT bool) *RGSW {
	return &RGSW{
		Body: NewGadgetEncryption(params, isNTT),
		Mask: NewGadgetEncryption(params, isNTT),
	}
}

// NewRGSWCustom creates a new [RGSW] with the given parameters.
func NewRGSWCustom(rank, modLen int, hasAux bool, gadLen int, isNTT bool) *RGSW {
	return &RGSW{
		Body: NewGadgetEncryptionCustom(rank, modLen, hasAux, gadLen, isNTT),
		Mask: NewGadgetEncryptionCustom(rank, modLen, hasAux, gadLen, isNTT),
	}
}

// Rank returns the rank.
func (ct *RGSW) Rank() int {
	if ct.Body.Rank() != ct.Mask.Rank() {
		panic("inconsistent rank")
	}
	return ct.Body.Rank()
}

// ModLen returns the modulus length.
func (ct *RGSW) ModLen() int {
	if ct.Body.ModLen() != ct.Mask.ModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.ModLen()
}

// IsNTT returns whether value is in NTT form.
func (ct *RGSW) IsNTT() bool {
	if ct.Body.IsNTT() != ct.Mask.IsNTT() {
		panic("inconsistent NTT form")
	}
	return ct.Body.IsNTT()
}

// HasAuxModulus returns whether auxiliary modulus is used.
func (ct *RGSW) HasAuxModulus() bool {
	if ct.Body.HasAuxModulus() != ct.Mask.HasAuxModulus() {
		panic("inconsistent auxiliary modulus")
	}
	return ct.Body.HasAuxModulus()
}

// GadgetLen returns the length of the gadget.
func (ct *RGSW) GadgetLen() int {
	if ct.Body.GadgetLen() != ct.Mask.GadgetLen() {
		panic("inconsistent gadget length")
	}
	return ct.Body.GadgetLen()
}

// Clear clears value.
func (ct *RGSW) Clear() {
	ct.Body.Clear()
	ct.Mask.Clear()
}

// WithModIdx returns a shallow copy with the given modulus indices.
// HasAux is preserved.
//
// Panics when idx is out of range.
func (ct *RGSW) WithModIdx(idx ...int) *RGSW {
	return &RGSW{
		Body: ct.Body.WithModIdx(idx...),
		Mask: ct.Mask.WithModIdx(idx...),
	}
}

// Copy returns a copy of ct.
func (ct *RGSW) Copy() *RGSW {
	return &RGSW{
		Body: ct.Body.Copy(),
		Mask: ct.Mask.Copy(),
	}
}

// CopyFrom copies the coefficients from ctIn to ct.
func (ct *RGSW) CopyFrom(ctIn *RGSW) {
	if !ct.IsConsistent(ctIn) {
		panic("inconsistent input(s)")
	}
	ct.Body.CopyFrom(ctIn.Body)
	ct.Mask.CopyFrom(ctIn.Mask)
}

// IsEqual checks if two values are equal.
func (ct *RGSW) IsEqual(ct0 *RGSW) bool {
	return ct.Body.IsEqual(ct0.Body) && ct.Mask.IsEqual(ct0.Mask)
}

// IsConsistent checks if two values have the same shape.
func (ct *RGSW) IsConsistent(ct0 *RGSW) bool {
	return ct.Body.IsConsistent(ct0.Body) && ct.Mask.IsConsistent(ct0.Mask)
}

// // RelinKey is a relinearisation key.
// type RelinKey GadgetEncryption

// // KeySwitchKey is a key switch key.
// type KeySwitchKey GadgetEncryption

// // AutomorphismKey is an automorphism key.
// type AutomorphismKey struct {
// 	// value is the gadget encryption.
// 	Value []*Ciphertext
// 	// idx is the automorphism index.
// 	Idx int
// }
