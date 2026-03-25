package bfv

import (
	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encoder encodes/decodes [*crt.Element] into/from [*Plaintext].
type Encoder struct {
	intEcd *heint.Encoder
}

// NewEncoder creates a new [Encoder].
func NewEncoder(params rlwe.Parameters, msgMod *num.Modulus) *Encoder {
	return &Encoder{
		intEcd: heint.NewEncoder(params, msgMod),
	}
}

// Encode encodes a [*Plaintext] into a [*rlwe.Element].
func (ecd *Encoder) Encode(eIn Plaintext, hasAux, isNTT bool) *rlwe.Element {
	return ecd.intEcd.Encode(eIn, hasAux, isNTT)
}

// EncodeCustom encodes a [*Plaintext] into a [*rlwe.Element] with custom parameters.
func (ecd *Encoder) EncodeCustom(eIn Plaintext, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	return ecd.intEcd.EncodeCustom(eIn, baseLen, auxLen, isNTT)
}

// EncodeTo encodes a [*Plaintext] into a [*rlwe.Element].
func (ecd *Encoder) EncodeTo(eOut *rlwe.Element, eIn Plaintext, isNTT bool) {
	ecd.intEcd.EncodeTo(eOut, eIn, isNTT)
}

// Decode decodes a [*rlwe.Element] into a [*Plaintext].
func (ecd *Encoder) Decode(e *rlwe.Element) Plaintext {
	return ecd.intEcd.Decode(e)
}

// DecodeTo decodes a [*rlwe.Element] into a [*Plaintext].
func (ecd *Encoder) DecodeTo(eOut Plaintext, e *rlwe.Element) {
	ecd.intEcd.DecodeTo(eOut, e)
}
