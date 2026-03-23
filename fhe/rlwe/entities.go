package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/vec"
)

type Element struct {
	Value  *crt.Element
	auxLen int
}

// NewScalar creates a new scalar [Element].
func NewScalar(params Parameters, hasAux bool) *Element {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	auxLen := len(params.auxMod)
	if !hasAux {
		auxLen = 0
	}

	return NewElement(1, len(params.baseMod), auxLen, false)
}

// NewPoly creates a new polynomial [Element].
func NewPoly(params Parameters, hasAux bool, isNTT bool) *Element {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	auxLen := len(params.auxMod)
	if !hasAux {
		auxLen = 0
	}

	return NewElement(params.RingParams().Rank(), len(params.baseMod), auxLen, isNTT)
}

// NewElement creates a new [Element] with the given parameters.
func NewElement(rank, baseLen, auxLen int, isNTT bool) *Element {
	return &Element{
		Value:  crt.NewPolyCustom(rank, baseLen+auxLen, isNTT),
		auxLen: auxLen,
	}
}

// NewElementFrom creates a new [Element] from a [crt.Element].
func NewElementFrom(value *crt.Element, auxLen int) *Element {
	return &Element{
		Value:  value,
		auxLen: auxLen,
	}
}

// Rank returns the rank.
func (e *Element) Rank() int {
	return e.Value.Rank()
}

// IsNTT returns whether value is in NTT form.
func (e *Element) IsNTT() bool {
	return e.Value.IsNTT
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

// Clear clears value.
func (e *Element) Clear() {
	e.Value.Clear()
}

// Resize resizes the element to the given length.
func (e *Element) Resize(baseLen, auxLen int) {
	switch {
	case e.auxLen >= auxLen && e.BaseModLen() >= baseLen:
		e.Value.Coeffs = e.Value.Coeffs[e.auxLen-auxLen : e.auxLen+baseLen]
	case e.auxLen >= auxLen && e.BaseModLen() < baseLen:
		extraBase := make([][]uint64, baseLen-e.BaseModLen())
		for i := range extraBase {
			extraBase[i] = make([]uint64, e.Value.Rank())
		}
		e.Value.Coeffs = append(e.Value.Coeffs, extraBase...)
	case e.auxLen < auxLen && e.BaseModLen() >= baseLen:
		extraAux := make([][]uint64, auxLen-e.auxLen)
		for i := range extraAux {
			extraAux[i] = make([]uint64, e.Value.Rank())
		}
		e.Value.Coeffs = append(extraAux, e.Value.Coeffs[e.auxLen-auxLen:e.auxLen+baseLen]...)
	case e.auxLen < auxLen && e.BaseModLen() < baseLen:
		extraAux := make([][]uint64, auxLen-e.auxLen)
		for i := range extraAux {
			extraAux[i] = make([]uint64, e.Value.Rank())
		}
		extraBase := make([][]uint64, baseLen-e.BaseModLen())
		for i := range extraBase {
			extraBase[i] = make([]uint64, e.Value.Rank())
		}

		e.Value.Coeffs = append(extraAux, extraBase...)
		e.Value.Coeffs = append(e.Value.Coeffs, extraBase...)
	}
	e.auxLen = auxLen
}

// Type returns the type of e.
func (e *Element) Type() crt.ElementType {
	return e.Value.Type()
}

// WithModIdx returns a shallow copy with the given modulus indices.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (e *Element) WithModLen(baseLen, auxLen int) *Element {
	if baseLen > e.BaseModLen() || auxLen > e.AuxModLen() {
		panic("out of range")
	}

	return &Element{
		Value:  e.Value.WithModIdx(vec.Range(e.auxLen-auxLen, e.auxLen+baseLen)...),
		auxLen: auxLen,
	}
}

// Copy returns a copy.
func (e *Element) Copy() *Element {
	return &Element{
		Value:  e.Value.Copy(),
		auxLen: e.auxLen,
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
	e.auxLen = eIn.auxLen
}

// IsEqual checks if two values are equal.
func (e *Element) IsEqual(e0 *Element) bool {
	if !e.IsConsistent(e0) {
		return false
	}
	return e.Value.IsEqual(e0.Value) && e.auxLen == e0.auxLen
}

// IsConsistent checks if two values have the same shape.
func (e *Element) IsConsistent(e0 *Element) bool {
	return e.Value.IsConsistent(e0.Value)
}

// SecretKey is an RLWE secret key.
type SecretKey Element

// NewSecretKey creates a new [SecretKey].
func NewSecretKey(params Parameters, hasAux bool) *SecretKey {
	if hasAux && len(params.auxMod) == 0 {
		panic("no auxiliary modulus")
	}

	return (*SecretKey)(NewPoly(params, hasAux, false))
}

// NewSecretKeyCustom creates a new [SecretKey] with the given parameters.
func NewSecretKeyCustom(rank, baseLen, auxLen int, isNTT bool) *SecretKey {
	return (*SecretKey)(NewElement(rank, baseLen, auxLen, isNTT))
}

// Rank returns the rank.
func (sk *SecretKey) Rank() int {
	return (*Element)(sk).Rank()
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

// IsNTT returns whether value is in NTT form.
func (sk *SecretKey) IsNTT() bool {
	return (*Element)(sk).IsNTT()
}

// Clear clears value.
func (sk *SecretKey) Clear() {
	(*Element)(sk).Clear()
}

// WithModIdx returns a shallow copy with the given modulus indices.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (sk *SecretKey) WithModLen(baseLen, auxLen int) *SecretKey {
	return (*SecretKey)((*Element)(sk).WithModLen(baseLen, auxLen))
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

	auxLen := len(params.auxMod)
	if !hasAux {
		auxLen = 0
	}

	return NewCiphertextCustom(params.RingParams().Rank(), len(params.baseMod), auxLen, isNTT)
}

// NewCiphertextCustom creates a new [Ciphertext] with the given parameters.
func NewCiphertextCustom(rank, baseLen, auxLen int, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Body: NewElement(rank, baseLen, auxLen, isNTT),
		Mask: NewElement(rank, baseLen, auxLen, isNTT),
	}
}

// Rank returns the rank.
func (ct *Ciphertext) Rank() int {
	if ct.Body.Rank() != ct.Mask.Rank() {
		panic("inconsistent rank")
	}
	return ct.Body.Rank()
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

// IsNTT returns whether value is in NTT form.
func (ct *Ciphertext) IsNTT() bool {
	if ct.Body.IsNTT() != ct.Mask.IsNTT() {
		panic("inconsistent NTT form")
	}
	return ct.Body.IsNTT()
}

// Resize resizes the ciphertext to the given length.
func (ct *Ciphertext) Resize(baseLen, auxLen int) {
	ct.Body.Resize(baseLen, auxLen)
	ct.Mask.Resize(baseLen, auxLen)
}

// Clear clears value.
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
func NewPublicKeyCustom(rank, baseLen, auxLen int, isNTT bool) *PublicKey {
	return (*PublicKey)(NewCiphertextCustom(rank, baseLen, auxLen, isNTT))
}

// Rank returns the rank.
func (pk *PublicKey) Rank() int {
	return (*Ciphertext)(pk).Rank()
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

// Clear clears value.
func (pk *PublicKey) Clear() {
	(*Ciphertext)(pk).Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (pk *PublicKey) WithModLen(baseLen, auxLen int) *PublicKey {
	return (*PublicKey)((*Ciphertext)(pk).WithModLen(baseLen, auxLen))

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

	auxLen := len(params.auxMod)
	if !hasAux {
		auxLen = 0
	}

	return NewVectorCustom(params.RingParams().Rank(), len(params.baseMod), auxLen, length, isNTT)
}

// NewVectorCustom creates a new [Vector] with the given parameters.
func NewVectorCustom(rank, baseLen, auxLen int, length int, isNTT bool) *Vector {
	value := make([]*Element, length)
	for i := range value {
		value[i] = NewElement(rank, baseLen, auxLen, isNTT)
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

// FullModLen returns the full modulus length of v.
func (v *Vector) FullModLen() int {
	fullLen := v.Value[0].FullModLen()
	for i := 1; i < v.Len(); i++ {
		if fullLen != v.Value[i].FullModLen() {
			panic("inconsistent modulus length")
		}
	}
	return fullLen
}

// BaseModLen returns the base modulus length of v.
func (v *Vector) BaseModLen() int {
	baseLen := v.Value[0].BaseModLen()
	for i := 1; i < v.Len(); i++ {
		if baseLen != v.Value[i].BaseModLen() {
			panic("inconsistent modulus length")
		}
	}
	return baseLen
}

// AuxModLen returns the auxiliary modulus length of v.
func (v *Vector) AuxModLen() int {
	auxLen := v.Value[0].AuxModLen()
	for i := 1; i < v.Len(); i++ {
		if auxLen != v.Value[i].AuxModLen() {
			panic("inconsistent modulus length")
		}
	}
	return auxLen
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

// Resize resizes the vector to the given length.
func (v *Vector) Resize(baseLen, auxLen int) {
	for _, e := range v.Value {
		e.Resize(baseLen, auxLen)
	}
}

// Clear clears value.
func (v *Vector) Clear() {
	for _, p := range v.Value {
		p.Clear()
	}
}

// Slice returns a shallow copy of the elements in the range idx.
//
// Panics when idx is out of range.
func (v *Vector) Slice(idx ...int) *Vector {
	length := len(idx)
	value := make([]*Element, length)
	for i := range value {
		if idx[i] < 0 || idx[i] >= v.Len() {
			panic("index out of range")
		}
		value[i] = v.Value[idx[i]]
	}

	return &Vector{Value: value}
}

// WithModLen returns a shallow copy with the given length, base modulus length and auxiliary modulus length.
//
// Panics when length is out of range.
func (v *Vector) WithModLen(baseLen, auxLen int) *Vector {
	value := make([]*Element, v.Len())
	for i := range value {
		value[i] = v.Value[i].WithModLen(baseLen, auxLen)
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
	return NewGadgetEncryptionCustom(params.RingParams().Rank(), len(params.baseMod), len(params.auxMod), params.GadgetLen(), isNTT)
}

// NewGadgetEncryptionCustom creates a new [GadgetEncryption] with the given parameters.
func NewGadgetEncryptionCustom(rank, baseLen, auxLen int, gadLen int, isNTT bool) *GadgetEncryption {
	value := make([]*Ciphertext, gadLen)
	for i := range value {
		value[i] = NewCiphertextCustom(rank, baseLen, auxLen, isNTT)
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

// FullModLen returns the full modulus length.
func (ct *GadgetEncryption) FullModLen() int {
	fullLen := ct.Value[0].FullModLen()
	for i := 1; i < len(ct.Value); i++ {
		if fullLen != ct.Value[i].FullModLen() {
			panic("inconsistent modulus length")
		}
	}
	return fullLen
}

// BaseModLen returns the base modulus length.
func (ct *GadgetEncryption) BaseModLen() int {
	baseLen := ct.Value[0].BaseModLen()
	for i := 1; i < len(ct.Value); i++ {
		if baseLen != ct.Value[i].BaseModLen() {
			panic("inconsistent modulus length")
		}
	}
	return baseLen
}

// AuxModLen returns the auxiliary modulus length.
func (ct *GadgetEncryption) AuxModLen() int {
	auxLen := ct.Value[0].AuxModLen()
	for i := 1; i < len(ct.Value); i++ {
		if auxLen != ct.Value[i].AuxModLen() {
			panic("inconsistent modulus length")
		}
	}
	return auxLen
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
func (ct *GadgetEncryption) WithModLen(baseLen, auxLen int) *GadgetEncryption {
	value := make([]*Ciphertext, len(ct.Value))
	for i := range value {
		value[i] = ct.Value[i].WithModLen(baseLen, auxLen)
	}
	return &GadgetEncryption{Value: value}
}

// Copy returns a copy.
func (ct *GadgetEncryption) Copy() *GadgetEncryption {
	value := make([]*Ciphertext, len(ct.Value))
	for i := range value {
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
	for i := range ct.Value {
		ct.Value[i].CopyFrom(ctIn.Value[i])
	}
}

// IsEqual checks if two values are equal.
func (ct *GadgetEncryption) IsEqual(ct0 *GadgetEncryption) bool {
	if !ct.IsConsistent(ct0) {
		return false
	}

	for i := range ct.Value {
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

	for i := range ct.Value {
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
func NewRGSWCustom(rank, baseLen, auxLen int, gadLen int, isNTT bool) *RGSW {
	return &RGSW{
		Body: NewGadgetEncryptionCustom(rank, baseLen, auxLen, gadLen, isNTT),
		Mask: NewGadgetEncryptionCustom(rank, baseLen, auxLen, gadLen, isNTT),
	}
}

// Rank returns the rank.
func (ct *RGSW) Rank() int {
	if ct.Body.Rank() != ct.Mask.Rank() {
		panic("inconsistent rank")
	}
	return ct.Body.Rank()
}

// FullModLen returns the full modulus length.
func (ct *RGSW) FullModLen() int {
	if ct.Body.FullModLen() != ct.Mask.FullModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.FullModLen()
}

// BaseModLen returns the base modulus length.
func (ct *RGSW) BaseModLen() int {
	if ct.Body.BaseModLen() != ct.Mask.BaseModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (ct *RGSW) AuxModLen() int {
	if ct.Body.AuxModLen() != ct.Mask.AuxModLen() {
		panic("inconsistent modulus length")
	}
	return ct.Body.AuxModLen()
}

// IsNTT returns whether value is in NTT form.
func (ct *RGSW) IsNTT() bool {
	if ct.Body.IsNTT() != ct.Mask.IsNTT() {
		panic("inconsistent NTT form")
	}
	return ct.Body.IsNTT()
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

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (ct *RGSW) WithModLen(baseLen, auxLen int) *RGSW {
	return &RGSW{
		Body: ct.Body.WithModLen(baseLen, auxLen),
		Mask: ct.Mask.WithModLen(baseLen, auxLen),
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

// RelinKey is a relinearisation key.
type RelinKey GadgetEncryption

// KeySwitchKey is a key switch key.
type KeySwitchKey GadgetEncryption

// AutomorphismKey is an automorphism key.
type AutomorphismKey struct {
	// value is the gadget encryption.
	Value []*Ciphertext
	// idx is the automorphism index.
	Idx int
}
