package bgv

import (
	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encoder encodes/decodes [*crt.Element] into/from [*Plaintext].
type Encoder struct {
	encoder *heint.Encoder
}

// NewEncoder creates a new [Encoder].
func NewEncoder(params rlwe.Parameters, msgMod *num.Modulus) *Encoder {
	return &Encoder{
		encoder: heint.NewEncoder(params, msgMod),
	}
}

// Encode encodes a [*Plaintext] into a [*rlwe.Element].
func (enc *Encoder) Encode(eIn *Plaintext, hasAux, isNTT bool) *rlwe.Element {
	return enc.encoder.Encode((*crt.Element)(eIn), hasAux, isNTT)
}

// EncodeCustom encodes a [*Plaintext] into a [*rlwe.Element] with custom parameters.
func (enc *Encoder) EncodeCustom(eIn *Plaintext, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	return enc.encoder.EncodeCustom((*crt.Element)(eIn), baseLen, auxLen, isNTT)
}

// EncodeTo encodes a [*Plaintext] into a [*rlwe.Element].
func (enc *Encoder) EncodeTo(eOut *rlwe.Element, eIn *Plaintext, isNTT bool) {
	enc.encoder.EncodeTo(eOut, (*crt.Element)(eIn), isNTT)
}

// Decode decodes a [*rlwe.Element] into a [*Plaintext].
func (enc *Encoder) Decode(e *rlwe.Element) *Plaintext {
	return (*Plaintext)(enc.encoder.Decode(e))
}

// DecodeTo decodes a [*rlwe.Element] into a [*Plaintext].
func (enc *Encoder) DecodeTo(eOut *Plaintext, e *rlwe.Element) {
	enc.encoder.DecodeTo((*crt.Element)(eOut), e)
}
