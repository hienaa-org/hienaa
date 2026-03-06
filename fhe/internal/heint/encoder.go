package heint

import (
	"sync"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encoder encodes/decodes [*crt.Element] into/from [*rlwe.Element].
type Encoder struct {
	rlweParams rlwe.Parameters
	msgMod     *num.Modulus

	pOp *rlwe.PlainOperator

	pPool *sync.Pool

	msgToFull []*crt.Embedder
	baseToMsg []*crt.Embedder
}

// NewEncoder creates a new [Encoder].
func NewEncoder(rlweParams rlwe.Parameters, msgMod *num.Modulus) *Encoder {
	pOp := rlwe.NewPlainOperator(rlweParams)

	msgToFull := make([]*crt.Embedder, 1+len(rlweParams.AuxModulus()))
	paramAuxLen := len(rlweParams.AuxModulus())
	fullMod := rlweParams.FullModulus()
	for i := range msgToFull {
		msgToFull[i] = crt.NewEmbedder(fullMod[paramAuxLen-i:], []*num.Modulus{msgMod})
	}

	baseToMsg := make([]*crt.Embedder, len(rlweParams.BaseModulus()))
	for i := range baseToMsg {
		baseToMsg[i] = crt.NewEmbedder([]*num.Modulus{msgMod}, rlweParams.BaseModulus()[:i+1])
	}

	return &Encoder{
		rlweParams: rlweParams,
		msgMod:     msgMod,

		pOp: pOp,

		pPool: &sync.Pool{
			New: func() any {
				return rlwe.NewPoly(rlweParams, rlweParams.HasAuxModulus(), true)
			},
		},

		msgToFull: msgToFull,
		baseToMsg: baseToMsg,
	}
}

// Encode encodes a [*crt.Element] into a [*rlwe.Element].
func (enc *Encoder) Encode(eIn *crt.Element, hasAux, isNTT bool) *rlwe.Element {
	auxLen := len(enc.rlweParams.AuxModulus())
	if !hasAux {
		auxLen = 0
	}

	eOut := rlwe.NewElement(eIn.Rank(), len(enc.rlweParams.BaseModulus()), auxLen, isNTT)
	enc.EncodeTo(eOut, eIn, isNTT)
	return eOut
}

// EncodeCustom encodes a [*crt.Element] into a [*rlwe.Element] with custom parameters.
func (enc *Encoder) EncodeCustom(eIn *crt.Element, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	eOut := rlwe.NewElement(eIn.Rank(), baseLen, auxLen, isNTT)
	enc.EncodeTo(eOut, eIn, isNTT)
	return eOut
}

// EncodeTo encodes a [*crt.Element] into a [*rlwe.Element].
func (enc *Encoder) EncodeTo(eOut *rlwe.Element, eIn *crt.Element, isNTT bool) {
	if eOut.Rank() != eIn.Rank() {
		panic("inconsistent input(s)")
	}

	auxLen := eOut.AuxModLen()
	enc.msgToFull[auxLen].EmbedTo(eOut.Value, (*crt.Element)(eIn))
	eOut.Value.IsNTT = false

	if isNTT && eOut.Value.Type() == crt.TypePoly {
		enc.pOp.FwdNTTTo((*rlwe.Element)(eOut), (*rlwe.Element)(eOut))
	}
}

// Decode decodes a [*rlwe.Element] into a [*crt.Element].
func (enc *Encoder) Decode(e *rlwe.Element) *crt.Element {
	pt := crt.NewPoly(e.Rank(), 1)
	enc.DecodeTo(pt, e)
	return pt
}

// DecodeTo decodes a [*rlwe.Element] into a [*crt.Element].
func (enc *Encoder) DecodeTo(eOut *crt.Element, e *rlwe.Element) {
	if eOut.Rank() != e.Rank() {
		panic("inconsistent input(s)")
	}

	baseLen := e.BaseModLen()
	pNTT := enc.pPool.Get().(*rlwe.Element)
	defer enc.pPool.Put(pNTT)
	pNTT = pNTT.WithModLen(baseLen, 0)

	eBase := (*rlwe.Element)(e.WithModLen(baseLen, 0))
	pNTT.CopyFrom(eBase)

	if e.IsNTT() {
		enc.pOp.InvNTTTo(pNTT, pNTT)
	}

	enc.baseToMsg[baseLen-1].EmbedTo((*crt.Element)(eOut), pNTT.Value)
}
