package rlwe

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// Scalar is a scalar in CRT representation.
type Scalar struct {
	// Value is the underlying [crt.Scalar].
	Value *crt.Scalar
	// auxLen is the auxiliary modulus length.
	auxLen int
}

// NewScalar creates a new [Scalar].
func NewScalar(params Parameters, withAux bool) *Scalar {
	if withAux && !params.HasAuxMod() {
		panic("no auxiliary modulus")
	}

	if withAux {
		return &Scalar{
			Value:  crt.NewScalar(len(params.fullMod)),
			auxLen: len(params.auxMod),
		}
	}
	return &Scalar{
		Value:  crt.NewScalar(len(params.baseMod)),
		auxLen: 0,
	}
}

// NewScalarCustom creates a new [Scalar] with the given parameters.
func NewScalarCustom(baseLen, auxLen int) *Scalar {
	return &Scalar{
		Value:  crt.NewScalar(baseLen + auxLen),
		auxLen: auxLen,
	}
}

// NewScalarFrom creates a new [Scalar] from [crt.Scalar].
func NewScalarFrom(value *crt.Scalar, auxLen int) *Scalar {
	return &Scalar{
		Value:  value,
		auxLen: auxLen,
	}
}

// FullModLen returns the full modulus length.
func (c *Scalar) FullModLen() int {
	return c.Value.ModLen()
}

// BaseModLen returns the base modulus length.
func (c *Scalar) BaseModLen() int {
	return c.Value.ModLen() - c.auxLen
}

// AuxModLen returns the auxiliary modulus length.
func (c *Scalar) AuxModLen() int {
	return c.auxLen
}

// Clear clears c.
func (c *Scalar) Clear() {
	c.Value.Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (c *Scalar) WithModLen(baseLen, auxLen int) *Scalar {
	if baseLen < 0 || auxLen < 0 {
		panic("negative modulus length")
	} else if baseLen > c.BaseModLen() || auxLen > c.AuxModLen() {
		panic("out of range")
	}

	return &Scalar{
		Value:  c.Value.Slice(c.auxLen-auxLen, c.auxLen+baseLen),
		auxLen: auxLen,
	}
}

// Resize resizes c in place to the given base and auxiliary modulus lengths.
// It preserves the auxiliary suffix and base prefix, and zeroes new coefficients.
//
// Panics when baseLen or auxLen is invalid.
func (c *Scalar) Resize(baseLen, auxLen int) {
	if baseLen <= 0 || auxLen < 0 {
		panic("invalid modulus length")
	}

	auxLenKeep := min(c.AuxModLen(), auxLen)
	baseLenKeep := min(c.BaseModLen(), baseLen)

	if auxLen+baseLen > len(c.Value.Value) {
		c.Value.Value = slices.Grow(c.Value.Value, auxLen+baseLen-len(c.Value.Value))[:auxLen+baseLen]
	}

	copy(c.Value.Value[auxLen-auxLenKeep:auxLen+baseLenKeep], c.Value.Value[c.auxLen-auxLenKeep:c.auxLen+baseLenKeep])
	clear(c.Value.Value[:auxLen-auxLenKeep])
	clear(c.Value.Value[auxLen+baseLenKeep:])

	c.Value.Value = c.Value.Value[:auxLen+baseLen]
	c.auxLen = auxLen
}

// Copy returns a copy of c.
func (c *Scalar) Copy() *Scalar {
	return &Scalar{
		Value:  c.Value.Copy(),
		auxLen: c.auxLen,
	}
}

// CopyFrom copies the coefficients from c0 to c.
//
// Panics when c and c0 are not consistent.
func (c *Scalar) CopyFrom(c0 *Scalar) {
	c.Value.CopyFrom(c0.Value)
	c.auxLen = c0.auxLen
}

// IsEqual checks if c is equal to c0.
func (c *Scalar) IsEqual(c0 *Scalar) bool {
	return c.Value.IsEqual(c0.Value) && c.auxLen == c0.auxLen
}

// IsConsistent checks if c has the same shape as c0.
func (c *Scalar) IsConsistent(c0 *Scalar) bool {
	return c.Value.IsConsistent(c0.Value)
}

// Poly is a polynomial in CRT representation.
type Poly struct {
	// Value is the underlying [crt.Poly].
	Value *crt.Poly
	// auxLen is the auxiliary modulus length.
	auxLen int
}

// NewPoly creates a new [Poly].
func NewPoly(params Parameters, withAux bool, form dft.Form) *Poly {
	if withAux && !params.HasAuxMod() {
		panic("no auxiliary modulus")
	}

	if withAux {
		return NewPolyCustom(params.Rank(), len(params.baseMod), len(params.auxMod), form)
	}
	return NewPolyCustom(params.Rank(), len(params.baseMod), 0, form)
}

// NewPolyCustom creates a new [Poly] with the given parameters.
func NewPolyCustom(rank, baseLen, auxLen int, form dft.Form) *Poly {
	return &Poly{
		Value:  crt.NewPoly(rank, baseLen+auxLen, form),
		auxLen: auxLen,
	}
}

// NewPolyFrom creates a new [Poly] from [crt.Poly].
func NewPolyFrom(value *crt.Poly, auxLen int) *Poly {
	return &Poly{
		Value:  value,
		auxLen: auxLen,
	}
}

// Rank returns the rank of p.
func (p *Poly) Rank() int {
	return p.Value.Rank()
}

// Form returns the [dft.Form] of p.
func (p *Poly) Form() dft.Form {
	return p.Value.Form
}

// FullModLen returns the full modulus length.
func (p *Poly) FullModLen() int {
	return p.Value.ModLen()
}

// BaseModLen returns the base modulus length.
func (p *Poly) BaseModLen() int {
	return p.Value.ModLen() - p.auxLen
}

// AuxModLen returns the auxiliary modulus length.
func (p *Poly) AuxModLen() int {
	return p.auxLen
}

// Clear clears p.
func (p *Poly) Clear() {
	p.Value.Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (p *Poly) WithModLen(baseLen, auxLen int) *Poly {
	if baseLen < 0 || auxLen < 0 {
		panic("negative modulus length")
	} else if baseLen > p.BaseModLen() || auxLen > p.AuxModLen() {
		panic("out of range")
	}

	return &Poly{
		Value:  p.Value.Slice(p.auxLen-auxLen, p.auxLen+baseLen),
		auxLen: auxLen,
	}
}

// Resize resizes p in place to the given base and auxiliary modulus lengths.
// It preserves the auxiliary suffix and base prefix, and zeroes new coefficients.
//
// Panics when baseLen or auxLen is invalid.
func (p *Poly) Resize(baseLen, auxLen int) {
	if baseLen <= 0 || auxLen < 0 {
		panic("negative modulus length")
	}

	rank := p.Rank()
	auxLenKeep := min(p.AuxModLen(), auxLen)
	baseLenKeep := min(p.BaseModLen(), baseLen)

	if auxLen+baseLen > len(p.Value.Coeffs) {
		p.Value.Coeffs = slices.Grow(p.Value.Coeffs, auxLen+baseLen-len(p.Value.Coeffs))[:auxLen+baseLen]
	}

	copy(p.Value.Coeffs[auxLen-auxLenKeep:auxLen+baseLenKeep], p.Value.Coeffs[p.auxLen-auxLenKeep:p.auxLen+baseLenKeep])
	for i := 0; i < auxLen-auxLenKeep; i++ {
		p.Value.Coeffs[i] = make([]uint64, rank)
	}
	for i := auxLen + baseLenKeep; i < auxLen+baseLen; i++ {
		p.Value.Coeffs[i] = make([]uint64, rank)
	}
	clear(p.Value.Coeffs[auxLen+baseLen:])

	p.Value.Coeffs = p.Value.Coeffs[:auxLen+baseLen]
	p.auxLen = auxLen
}

// Copy returns a copy of p.
func (p *Poly) Copy() *Poly {
	return &Poly{
		Value:  p.Value.Copy(),
		auxLen: p.auxLen,
	}
}

// CopyFrom copies the coefficients from p0 to p.
//
// Panics when p and p0 are not consistent.
func (p *Poly) CopyFrom(p0 *Poly) {
	p.Value.CopyFrom(p0.Value)
	p.auxLen = p0.auxLen
}

// IsEqual checks if p is equal to p0.
func (p *Poly) IsEqual(p0 *Poly) bool {
	return p.Value.IsEqual(p0.Value) && p.auxLen == p0.auxLen
}

// IsConsistent checks if p has the same shape as p0.
func (p *Poly) IsConsistent(p0 *Poly) bool {
	return p.Value.IsConsistent(p0.Value)
}

// SecretKey is an RLWE secret key.
type SecretKey Poly

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
	return (*Poly)(sk).Rank()
}

// Form returns the [dft.Form] of sk.
func (sk *SecretKey) Form() dft.Form {
	return (*Poly)(sk).Form()
}

// FullModLen returns the full modulus length.
func (sk *SecretKey) FullModLen() int {
	return (*Poly)(sk).FullModLen()
}

// BaseModLen returns the base modulus length.
func (sk *SecretKey) BaseModLen() int {
	return (*Poly)(sk).BaseModLen()
}

// AuxModLen returns the auxiliary modulus length.
func (sk *SecretKey) AuxModLen() int {
	return (*Poly)(sk).AuxModLen()
}

// Clear clears sk.
func (sk *SecretKey) Clear() {
	(*Poly)(sk).Clear()
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (sk *SecretKey) WithModLen(baseLen, auxLen int) *SecretKey {
	return (*SecretKey)((*Poly)(sk).WithModLen(baseLen, auxLen))
}

// Copy returns a copy of sk.
func (sk *SecretKey) Copy() *SecretKey {
	return (*SecretKey)((*Poly)(sk).Copy())
}

// CopyFrom copies the coefficients from sk0 to sk.
//
// Panics when sk and sk0 are not consistent.
func (sk *SecretKey) CopyFrom(sk0 *SecretKey) {
	(*Poly)(sk).CopyFrom((*Poly)(sk0))
}

// IsEqual checks if sk is equal to sk0.
func (sk *SecretKey) IsEqual(sk0 *SecretKey) bool {
	return (*Poly)(sk).IsEqual((*Poly)(sk0))
}

// IsConsistent checks if sk has the same shape as sk0.
func (sk *SecretKey) IsConsistent(sk0 *SecretKey) bool {
	return (*Poly)(sk).IsConsistent((*Poly)(sk0))
}

// Ciphertext is an RLWE ciphertext.
type Ciphertext struct {
	Body, Mask *Poly
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
		panic("inconsistent form")
	}
	return ct.Body.Form()
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
// It preserves the auxiliary suffix and base prefix, and zeroes new coefficients.
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

// CopyFrom copies the coefficients from ct0 to ct.
//
// Panics when ct and ct0 are not consistent.
func (ct *Ciphertext) CopyFrom(ct0 *Ciphertext) {
	if !ct.IsConsistent(ct0) {
		panic("inconsistent input(s)")
	}
	ct.Body.CopyFrom(ct0.Body)
	ct.Mask.CopyFrom(ct0.Mask)
}

// IsEqual checks if ct is equal to ct0.
func (ct *Ciphertext) IsEqual(ct0 *Ciphertext) bool {
	return ct.Body.IsEqual(ct0.Body) && ct.Mask.IsEqual(ct0.Mask)
}

// IsConsistent checks if ct has the same shape as ct0.
func (ct *Ciphertext) IsConsistent(ct0 *Ciphertext) bool {
	return ct.Body.IsConsistent(ct0.Body) && ct.Mask.IsConsistent(ct0.Mask)
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

// CopyFrom copies the coefficients from pk0 to pk.
//
// Panics when pk and pk0 are not consistent.
func (pk *PublicKey) CopyFrom(pk0 *PublicKey) {
	(*Ciphertext)(pk).CopyFrom((*Ciphertext)(pk0))
}

// IsEqual checks if pk is equal to pk0.
func (pk *PublicKey) IsEqual(pk0 *PublicKey) bool {
	return (*Ciphertext)(pk).IsEqual((*Ciphertext)(pk0))
}

// IsConsistent checks if pk has the same shape as pk0.
func (pk *PublicKey) IsConsistent(pk0 *PublicKey) bool {
	return (*Ciphertext)(pk).IsConsistent((*Ciphertext)(pk0))
}

// Vector is a vector of [Poly].
type Vector struct {
	Value []*Poly
}

// NewVector creates a new [Vector].
func NewVector(params Parameters, length int, withAux bool, form dft.Form) *Vector {
	if withAux && !params.HasAuxMod() {
		panic("no auxiliary modulus")
	}

	if withAux {
		return NewVectorCustom(params.Rank(), len(params.baseMod), len(params.auxMod), length, form)
	}
	return NewVectorCustom(params.Rank(), len(params.baseMod), 0, length, form)
}

// NewVectorCustom creates a new [Vector] with the given parameters.
func NewVectorCustom(rank, baseLen, auxLen, length int, form dft.Form) *Vector {
	value := make([]*Poly, length)
	for i := range value {
		value[i] = NewPolyCustom(rank, baseLen, auxLen, form)
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
			panic("inconsistent form")
		}
	}
	return form
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
//
// Panics when idx is out of range.
func (v *Vector) Slice(idx ...int) *Vector {
	length := len(idx)
	value := make([]*Poly, length)
	for i := range value {
		if idx[i] < 0 || idx[i] >= v.Len() {
			panic("index out of range")
		}
		value[i] = v.Value[idx[i]]
	}

	return &Vector{Value: value}
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (v *Vector) WithModLen(baseLen, auxLen int) *Vector {
	if baseLen < 0 || auxLen < 0 {
		panic("negative modulus length")
	}

	value := make([]*Poly, v.Len())
	for i := range value {
		value[i] = v.Value[i].WithModLen(baseLen, auxLen)
	}

	return &Vector{Value: value}
}

// Resize resizes v in place to the given base and auxiliary modulus lengths.
// It preserves the auxiliary suffix and base prefix, and zeroes new coefficients.
//
// Panics when baseLen or auxLen is invalid.
func (v *Vector) Resize(baseLen, auxLen int) {
	if baseLen <= 0 || auxLen < 0 {
		panic("negative modulus length")
	}

	for i := range v.Value {
		v.Value[i].Resize(baseLen, auxLen)
	}
}

// Copy returns a copy of v.
func (v *Vector) Copy() *Vector {
	value := make([]*Poly, v.Len())
	for i := 0; i < v.Len(); i++ {
		value[i] = v.Value[i].Copy()
	}

	return &Vector{Value: value}
}

// CopyFrom copies the coefficients from v0 to v.
//
// Panics when v and v0 are not consistent.
func (v *Vector) CopyFrom(v0 *Vector) {
	if !v.IsConsistent(v0) {
		panic("inconsistent input(s)")
	}

	for i := 0; i < v.Len(); i++ {
		v.Value[i].CopyFrom(v0.Value[i])
	}
}

// IsEqual checks if v is equal to v0.
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

// IsConsistent checks if v has the same shape as v0.
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
func NewGadgetEncryption(params Parameters, form dft.Form) *GadgetEncryption {
	return NewGadgetEncryptionCustom(params.Rank(), len(params.baseMod), len(params.auxMod), params.GadgetLen(), form)
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
			panic("inconsistent form")
		}
	}
	return form
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

// Slice returns a shallow copy with the given element indices.
//
// Panics when idx is out of range.
func (ct *GadgetEncryption) Slice(idx ...int) *GadgetEncryption {
	length := len(idx)
	value := make([]*Ciphertext, length)
	for i := range value {
		if idx[i] < 0 || idx[i] >= ct.GadgetLen() {
			panic("index out of range")
		}
		value[i] = ct.Value[idx[i]]
	}

	return &GadgetEncryption{Value: value}
}

// WithModLen returns a shallow copy with the given modulus lengths.
//
// Panics when baseLen or auxLen is larger than the current base or auxiliary modulus lengths.
func (ct *GadgetEncryption) WithModLen(baseLen, auxLen int) *GadgetEncryption {
	if baseLen < 0 || auxLen < 0 {
		panic("negative modulus length")
	}

	value := make([]*Ciphertext, ct.GadgetLen())
	for i := range value {
		value[i] = ct.Value[i].WithModLen(baseLen, auxLen)
	}

	return &GadgetEncryption{Value: value}
}

// Resize resizes ct in place to the given base and auxiliary modulus lengths.
// It preserves the auxiliary suffix and base prefix, and zeroes new coefficients.
//
// Panics when baseLen or auxLen is negative.
func (ct *GadgetEncryption) Resize(baseLen, auxLen int) {
	if baseLen < 0 || auxLen < 0 {
		panic("negative modulus length")
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

// CopyFrom copies the coefficients from ct0 to ct.
//
// Panics when ct and ct0 are not consistent.
func (ct *GadgetEncryption) CopyFrom(ct0 *GadgetEncryption) {
	if !ct.IsConsistent(ct0) {
		panic("inconsistent input(s)")
	}

	for i := 0; i < ct.GadgetLen(); i++ {
		ct.Value[i].CopyFrom(ct0.Value[i])
	}
}

// IsEqual checks if ct is equal to ct0.
func (ct *GadgetEncryption) IsEqual(ct0 *GadgetEncryption) bool {
	if !ct.IsConsistent(ct0) {
		return false
	}

	for i := 0; i < ct.GadgetLen(); i++ {
		if !ct.Value[i].IsEqual(ct0.Value[i]) {
			return false
		}
	}
	return true
}

// IsConsistent checks if ct has the same shape as ct0.
func (ct *GadgetEncryption) IsConsistent(ct0 *GadgetEncryption) bool {
	if ct.GadgetLen() != ct0.GadgetLen() {
		return false
	}

	for i := 0; i < ct.GadgetLen(); i++ {
		if !ct.Value[i].IsConsistent(ct0.Value[i]) {
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
		panic("inconsistent form")
	}
	return ct.Body.Form()
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

// Slice returns a shallow copy with the given element indices.
//
// Panics when idx is out of range.
func (ct *RGSWEncryption) Slice(idx ...int) *RGSWEncryption {
	return &RGSWEncryption{
		Body: ct.Body.Slice(idx...),
		Mask: ct.Mask.Slice(idx...),
	}
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
// It preserves the auxiliary suffix and base prefix, and zeroes new coefficients.
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

// CopyFrom copies the coefficients from ct0 to ct.
//
// Panics when ct and ct0 are not consistent.
func (ct *RGSWEncryption) CopyFrom(ct0 *RGSWEncryption) {
	if !ct.IsConsistent(ct0) {
		panic("inconsistent input(s)")
	}

	ct.Body.CopyFrom(ct0.Body)
	ct.Mask.CopyFrom(ct0.Mask)
}

// IsEqual checks if ct is equal to ct0.
func (ct *RGSWEncryption) IsEqual(ct0 *RGSWEncryption) bool {
	return ct.Body.IsEqual(ct0.Body) && ct.Mask.IsEqual(ct0.Mask)
}

// IsConsistent checks if ct has the same shape as ct0.
func (ct *RGSWEncryption) IsConsistent(ct0 *RGSWEncryption) bool {
	return ct.Body.IsConsistent(ct0.Body) && ct.Mask.IsConsistent(ct0.Mask)
}
