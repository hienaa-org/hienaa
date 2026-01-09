package rlwe

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/crt"
)

// SecretKey is a RLWE secret key.
type SecretKey PlainPoly

// PublicKey is a RLWE public key.
type PublicKey Ciphertext

// Plaintext is a RLWE plaintext.
type Plaintext interface {
	*PlainScalar | *PlainPoly
}

// PlainScalar is a RLWE plaintext scalar.
type PlainScalar struct {
	Value  crt.Scalar
	HasAux bool
}

// NewPlainScalar creates a new [PlainScalar] with the given parameters, auxiliary flag, and NTT flag.
func NewPlainScalar(params Parameters, hasAux bool) *PlainScalar {
	modLen := len(params.modulus)
	auxLen := len(params.auxModulus)

	return NewPlainScalarCustom(modLen, auxLen)
}

// NewPlainScalarCustom creates a new [PlainScalar] with the given modulus length and auxiliary modulus length.
func NewPlainScalarCustom(modLen int, auxLen int) *PlainScalar {
	return &PlainScalar{
		Value:  make([]uint64, modLen+auxLen),
		HasAux: auxLen > 0,
	}
}

// Clear clears s.
func (s *PlainScalar) Clear() {
	clear(s.Value)
}

// WithModIdx returns a copy of s with the given modulus indices.
//
// Panics when idx is out of range.
func (s *PlainScalar) WithModIdx(idx ...int) *PlainScalar {
	for _, idxi := range idx {
		if idxi < 0 || idxi >= len(s.Value) {
			panic("WithModIdx: index out of range")
		}
	}

	value := make([]uint64, len(idx))
	for i, idxi := range idx {
		value[i] = s.Value[idxi]
	}

	return &PlainScalar{
		Value:  value,
		HasAux: s.HasAux,
	}
}

// Copy returns a copy of s.
func (s *PlainScalar) Copy() *PlainScalar {
	value := make([]uint64, len(s.Value))
	copy(value, s.Value)

	return &PlainScalar{
		Value:  value,
		HasAux: s.HasAux,
	}
}

// CopyFrom copies the coefficients from sIn to s.
func (s *PlainScalar) CopyFrom(sIn *PlainScalar) {
	if !s.IsConsistent(sIn) {
		panic("CopyFrom: inconsistent plain scalars")
	}

	copy(s.Value, sIn.Value)
}

// ModLen returns the modulus length of s.
func (s *PlainScalar) ModLen() int {
	return len(s.Value)
}

// IsEqual checks if s is equal to s0.
func (s *PlainScalar) IsEqual(s0 *PlainScalar) bool {
	if !s.IsConsistent(s0) {
		return false
	}

	return slices.Equal(s.Value, s0.Value)
}

// IsConsistent checks if s has the same shape as s0.
func (s *PlainScalar) IsConsistent(s0 *PlainScalar) bool {
	return s.HasAux == s0.HasAux && len(s.Value) == len(s0.Value)
}

// PlainPoly is a RLWE plaintext polynomial.
type PlainPoly struct {
	Value  *crt.Poly
	HasAux bool
}

// NewPlainPoly creates a new [PlainPoly] with the given parameters, auxiliary flag, and NTT flag.
func NewPlainPoly(params Parameters, hasAux bool, isNTT bool) *PlainPoly {
	rank := params.ringParams.Rank()
	modLen := len(params.modulus)
	auxLen := len(params.auxModulus)

	return NewPlainPolyCustom(rank, modLen, auxLen, isNTT)
}

// NewPlainPolyCustom creates a new [PlainPoly] with the given rank, modulus length, auxiliary modulus length, and NTT flag.
func NewPlainPolyCustom(rank int, modLen int, auxLen int, isNTT bool) *PlainPoly {
	return &PlainPoly{
		Value:  crt.NewPolyCustom(rank, modLen+auxLen, isNTT),
		HasAux: auxLen > 0,
	}
}

// Clear clears pt.
func (pt *PlainPoly) Clear() {
	pt.Value.Clear()
}

// WithModIdx returns a copy of pt with the given modulus indices.
func (pt *PlainPoly) WithModIdx(idx ...int) *PlainPoly {
	return &PlainPoly{
		Value:  pt.Value.WithModIdx(idx...),
		HasAux: pt.HasAux,
	}
}

// Copy returns a copy of pt.
func (pt *PlainPoly) Copy() *PlainPoly {
	return &PlainPoly{
		Value:  pt.Value.Copy(),
		HasAux: pt.HasAux,
	}
}

// CopyFrom copies the coefficients from ptIn to pt.
//
// Panics when pt and ptIn are not consistent.
func (pt *PlainPoly) CopyFrom(ptIn *PlainPoly) {
	if !pt.IsConsistent(ptIn) {
		panic("CopyFrom: inconsistent plaintexts")
	}
	pt.Value.CopyFrom(ptIn.Value)
}

// ModLen returns the modulus length of pt.
func (pt *PlainPoly) ModLen() int {
	return pt.Value.ModLen()
}

// IsEqual checks if pt is equal to pt0.
func (pt *PlainPoly) IsEqual(pt0 *PlainPoly) bool {
	if !pt.IsConsistent(pt0) {
		return false
	}
	return pt.Value.IsEqual(pt0.Value)
}

// IsConsistent checks if pt has the same shape as pt0.
func (pt *PlainPoly) IsConsistent(pt0 *PlainPoly) bool {
	return pt.Value.IsConsistent(pt0.Value)
}

// Ciphertext is a RLWE ciphertext.
type Ciphertext struct {
	Body   *crt.Poly
	Mask   *crt.Poly
	HasAux bool
}

// NewCiphertext creates a new [Ciphertext] with the given parameters, auxiliary flag, and NTT flag.
func NewCiphertext(params Parameters, hasAux bool, isNTT bool) *Ciphertext {
	rank := params.ringParams.Rank()
	modLen := len(params.modulus)

	var auxLen int
	if hasAux {
		if params.auxModulus == nil {
			panic("auxiliary modulus is not set")
		}
		auxLen = len(params.auxModulus)
	}

	return NewCiphertextCustom(rank, modLen, auxLen, isNTT)
}

// NewCiphertextCustom creates a new [Ciphertext] with the given rank, modulus length, auxiliary modulus length, and NTT flag.
func NewCiphertextCustom(rank int, modLen int, auxLen int, isNTT bool) *Ciphertext {
	return &Ciphertext{
		Body:   crt.NewPolyCustom(rank, modLen+auxLen, isNTT),
		Mask:   crt.NewPolyCustom(rank, modLen+auxLen, isNTT),
		HasAux: auxLen > 0,
	}
}

// Clear clears c.
func (c *Ciphertext) Clear() {
	c.Body.Clear()
	c.Mask.Clear()
}

// WithModIdx returns a copy of c with the given modulus indices.
func (c *Ciphertext) WithModIdx(idx ...int) *Ciphertext {
	return &Ciphertext{
		Body:   c.Body.WithModIdx(idx...),
		Mask:   c.Mask.WithModIdx(idx...),
		HasAux: c.HasAux,
	}
}

// Copy returns a copy of c.
func (c *Ciphertext) Copy() *Ciphertext {
	return &Ciphertext{
		Body:   c.Body.Copy(),
		Mask:   c.Mask.Copy(),
		HasAux: c.HasAux,
	}
}

// CopyFrom copies the coefficients from cIn to c.
//
// Panics when c and cIn are not consistent.
func (c *Ciphertext) CopyFrom(cIn *Ciphertext) {
	if !c.IsConsistent(cIn) {
		panic("CopyFrom: inconsistent ciphertexts")
	}
	c.Body.CopyFrom(cIn.Body)
	c.Mask.CopyFrom(cIn.Mask)
}

// ModLen returns the modulus length of c.
func (c *Ciphertext) ModLen() int {
	modLen := len(c.Body.Coeffs)
	if modLen != len(c.Mask.Coeffs) {
		panic("ModLen: inconsistent modulus length")
	}
	return modLen
}

// IsEqual checks if c is equal to c0.
func (c *Ciphertext) IsEqual(c0 *Ciphertext) bool {
	if !c.IsConsistent(c0) {
		return false
	}
	return c.Body.IsEqual(c0.Body) && c.Mask.IsEqual(c0.Mask)
}

// IsConsistent checks if c has the same shape as c0.
func (c *Ciphertext) IsConsistent(c0 *Ciphertext) bool {
	return c.Body.IsConsistent(c0.Body) && c.Mask.IsConsistent(c0.Mask) && c.HasAux == c0.HasAux
}

// IsCompatible checks if c is compatible with p0.
func (c *Ciphertext) IsCompatible(p0 *PlainPoly) bool {
	return c.Body.IsConsistent(p0.Value) && c.Mask.IsConsistent(p0.Value) && c.HasAux == p0.HasAux
}

func (c *Ciphertext) IsCompatibleScalar(s *PlainScalar) bool {
	return c.Body.ModLen() == s.ModLen() && c.HasAux == s.HasAux
}

// CheckSanity checks if c is sane.
func (c *Ciphertext) CheckSanity() bool {
	if c.Body.Rank() != c.Mask.Rank() {
		return false
	} else if c.Body.ModLen() != c.Mask.ModLen() {
		return false
	} else if c.Body.IsNTT != c.Mask.IsNTT {
		return false
	}
	return true
}

// Tensor is a vector of polynomials.
type Tensor struct {
	Value  []*crt.Poly
	HasAux bool
}

// NewTensor creates a new [Tensor] with the given parameters, auxiliary flag, degree, and NTT flag.
func NewTensor(params Parameters, hasAux bool, degree int, isNTT bool) *Tensor {
	rank := params.ringParams.Rank()
	modLen := len(params.modulus)
	auxLen := len(params.auxModulus)

	return NewTensorCustom(rank, modLen, auxLen, degree, isNTT)
}

// NewTensorCustom creates a new [Tensor] with the given rank, modulus length, auxiliary modulus length, degree, and NTT flag.
func NewTensorCustom(rank int, modLen int, auxLen int, degree int, isNTT bool) *Tensor {
	value := make([]*crt.Poly, degree)
	for i := 0; i < degree; i++ {
		value[i] = crt.NewPolyCustom(rank, modLen+auxLen, isNTT)
	}
	return &Tensor{Value: value, HasAux: auxLen > 0}
}

// Degree returns the degree of t.
func (t *Tensor) Degree() int {
	return len(t.Value)
}

// Clear clears t.
func (t *Tensor) Clear() {
	for _, p := range t.Value {
		p.Clear()
	}
}

// WithModIdx returns a copy of t with the given modulus indices.
func (t *Tensor) WithDegreeAndModIdx(degree int, idx ...int) *Tensor {
	if degree > t.Degree() || degree < 1 {
		panic("Degree out of range")
	}

	value := make([]*crt.Poly, degree)
	for i := 0; i < degree; i++ {
		value[i] = t.Value[i].WithModIdx(idx...)
	}
	return &Tensor{Value: value, HasAux: t.HasAux}
}

// Copy returns a copy of t.
func (t *Tensor) Copy() *Tensor {
	value := make([]*crt.Poly, t.Degree())
	for i := 0; i < t.Degree(); i++ {
		value[i] = t.Value[i].Copy()
	}
	return &Tensor{Value: value, HasAux: t.HasAux}
}

// CopyFrom copies the coefficients from tIn to t.
//
// Panics when t and tIn are not consistent.
func (t *Tensor) CopyFrom(tIn *Tensor) {
	if !t.IsConsistent(tIn) {
		panic("CopyFrom: inconsistent tensor")
	}
}

// ModLen returns the modulus length of t.
func (t *Tensor) ModLen() int {
	modLen := len(t.Value[0].Coeffs)
	for i := 1; i < t.Degree(); i++ {
		if modLen != len(t.Value[i].Coeffs) {
			panic("ModLen: inconsistent modulus length")
		}
	}
	return modLen
}

// IsEqual checks if t is equal to t0.
func (t *Tensor) IsEqual(t0 *Tensor) bool {
	if !t.IsConsistent(t0) {
		return false
	}
	b := true
	for i := 0; i < t.Degree(); i++ {
		b = b && t.Value[i].IsEqual(t0.Value[i])
	}
	return b
}

// IsConsistent checks if t has the same shape as t0.
func (t *Tensor) IsConsistent(t0 *Tensor) bool {
	if t.Degree() != t0.Degree() {
		return false
	}

	b := true
	for i := 0; i < t.Degree(); i++ {
		b = b && t.Value[i].IsConsistent(t0.Value[i])
	}
	return b
}

// GadgetEncryption is a gadget encryption.
type GadgetEncryption struct {
	Value []*Ciphertext
}

// NewGadgetEncryption creates a new [GadgetEncryption] with the given parameters, degree, and NTT flag.
func NewGadgetEncryption(params Parameters, isNTT bool) *GadgetEncryption {
	gadLen := maxGadgetLen(params)
	value := make([]*Ciphertext, gadLen)
	for i := 0; i < gadLen; i++ {
		value[i] = NewCiphertext(params, true, isNTT)
	}
	return &GadgetEncryption{Value: value}
}

// Clear clears g.
func (g *GadgetEncryption) Clear() {
	for _, c := range g.Value {
		c.Clear()
	}
}

// Copy returns a copy of g.
func (g *GadgetEncryption) Copy() *GadgetEncryption {
	value := make([]*Ciphertext, len(g.Value))
	for i := 0; i < len(g.Value); i++ {
		value[i] = g.Value[i].Copy()
	}
	return &GadgetEncryption{Value: value}
}

// CopyFrom copies the coefficients from gIn to g.
//
// Panics when g and gIn are not consistent.
func (g *GadgetEncryption) CopyFrom(gIn *GadgetEncryption) {
	if !g.IsConsistent(gIn) {
		panic("CopyFrom: inconsistent gadget encryption")
	}
	for i := 0; i < len(g.Value); i++ {
		g.Value[i].CopyFrom(gIn.Value[i])
	}
}

// ModLen returns the modulus length of g.
func (g *GadgetEncryption) ModLen() int {
	modLen := g.Value[0].ModLen()
	for i := 1; i < len(g.Value); i++ {
		if modLen != g.Value[i].ModLen() {
			panic("ModLen: inconsistent modulus length")
		}
	}
	return modLen
}

// GadgetLen returns the length of the gadget.
func (g *GadgetEncryption) GadgetLen() int {
	return len(g.Value)
}

// IsEqual checks if g is equal to g0.
func (g *GadgetEncryption) IsEqual(g0 *GadgetEncryption) bool {
	if !g.IsConsistent(g0) {
		return false
	}
	b := true
	for i := 0; i < len(g.Value); i++ {
		b = b && g.Value[i].IsEqual(g0.Value[i])
	}
	return b
}

// IsConsistent checks if g has the same shape as g0.
func (g *GadgetEncryption) IsConsistent(g0 *GadgetEncryption) bool {
	if len(g.Value) != len(g0.Value) {
		return false
	}
	b := true
	for i := 0; i < len(g.Value); i++ {
		b = b && g.Value[i].IsConsistent(g0.Value[i])
	}
	return b
}

// RGSW is a Ring-GSW ciphertext.
type RGSW struct {
	Body *GadgetEncryption
	Mask *GadgetEncryption
}

// NewRGSW creates a new [RGSW] with the given parameters and NTT flag.
func NewRGSW(params Parameters, isNTT bool) *RGSW {
	return &RGSW{
		Body: NewGadgetEncryption(params, isNTT),
		Mask: NewGadgetEncryption(params, isNTT),
	}
}

// Clear clears r.
func (r *RGSW) Clear() {
	r.Body.Clear()
	r.Mask.Clear()
}

// Copy returns a copy of r.
func (g *RGSW) Copy() *RGSW {
	return &RGSW{
		Body: g.Body.Copy(),
		Mask: g.Mask.Copy(),
	}
}

// CopyFrom copies the coefficients from rIn to r.
func (r *RGSW) CopyFrom(rIn *RGSW) {
	if !r.IsConsistent(rIn) {
		panic("CopyFrom: inconsistent RGSW")
	}
	r.Body.CopyFrom(rIn.Body)
	r.Mask.CopyFrom(rIn.Mask)
}

// ModLen returns the modulus length of r.
func (r *RGSW) ModLen() int {
	return r.Body.ModLen()
}

// GadgetLen returns the length of the gadget.
func (r *RGSW) GadgetLen() int {
	return r.Body.GadgetLen()
}

// IsEqual checks if r is equal to r0.
func (r *RGSW) IsEqual(r0 *RGSW) bool {
	return r.Body.IsEqual(r0.Body) && r.Mask.IsEqual(r0.Mask)
}

// IsConsistent checks if r has the same shape as r0.
func (r *RGSW) IsConsistent(r0 *RGSW) bool {
	return r.Body.IsConsistent(r0.Body) && r.Mask.IsConsistent(r0.Mask)
}
