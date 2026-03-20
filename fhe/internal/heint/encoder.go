package heint

import (
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Encoder encodes/decodes [*crt.Element] into/from [*rlwe.Element].
type Encoder struct {
	rlweParams rlwe.Parameters
	msgMod     *num.Modulus

	pOp *rlwe.PlainOperator

	ePool *rlwe.ElementPool

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

		ePool: rlwe.NewElementPool(rlweParams, true, true),

		msgToFull: msgToFull,
		baseToMsg: baseToMsg,
	}
}

// Encode encodes a []uint64 into a [*rlwe.Element].
func (enc *Encoder) Encode(eIn []uint64, hasAux, isNTT bool) *rlwe.Element {
	auxLen := len(enc.rlweParams.AuxModulus())
	if !hasAux {
		auxLen = 0
	}

	eOut := rlwe.NewElement(len(eIn), len(enc.rlweParams.BaseModulus()), auxLen, isNTT)
	enc.EncodeTo(eOut, eIn, isNTT)
	return eOut
}

// EncodeCustom encodes a []uint64 into a [*rlwe.Element] with custom parameters.
func (enc *Encoder) EncodeCustom(eIn []uint64, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	eOut := rlwe.NewElement(len(eIn), baseLen, auxLen, isNTT)
	enc.EncodeTo(eOut, eIn, isNTT)
	return eOut
}

// EncodeTo encodes a []uint64 into a [*rlwe.Element].
func (enc *Encoder) EncodeTo(eOut *rlwe.Element, eIn []uint64, isNTT bool) {
	if len(eIn) != eOut.Rank() {
		panic("inconsistent input(s)")
	}

	var buf *rlwe.Element
	switch eOut.Rank() {
	case 1:
		buf = enc.ePool.Get(crt.TypeScalar)
	case enc.rlweParams.Rank():
		buf = enc.ePool.Get(crt.TypePoly)
	default:
		panic("invalid input length")
	}
	defer enc.ePool.Put(buf)
	buf = buf.WithModLen(1, 0)
	buf.Value.IsNTT = false
	vec.ReduceTo(buf.Value.Coeffs[0], eIn, enc.msgMod)

	auxLen := eOut.AuxModLen()
	enc.msgToFull[auxLen].EmbedTo(eOut.Value, buf.Value)
	eOut.Value.IsNTT = false

	if isNTT && eOut.Value.Type() == crt.TypePoly {
		enc.pOp.FwdNTTTo(eOut, eOut)
	}
}

// Decode decodes a [*rlwe.Element] into a []uint64.
func (enc *Encoder) Decode(e *rlwe.Element) []uint64 {
	eOut := make([]uint64, e.Rank())
	enc.DecodeTo(eOut, e)
	return eOut
}

// DecodeTo decodes a [*rlwe.Element] into a []uint64.
func (enc *Encoder) DecodeTo(eOut []uint64, e *rlwe.Element) {
	if e.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	} else if len(eOut) != e.Rank() {
		panic("invalid output length")
	}

	baseLen := e.BaseModLen()

	buf := enc.ePool.Get(e.Type())
	defer enc.ePool.Put(buf)
	buf = buf.WithModLen(baseLen, 0)
	bufOut := buf.WithModLen(1, 0)

	if e.IsNTT() {
		enc.pOp.InvNTTTo(buf, e)
		enc.baseToMsg[baseLen-1].EmbedTo(bufOut.Value, buf.Value)
	} else {
		enc.baseToMsg[baseLen-1].EmbedTo(bufOut.Value, e.Value)
	}

	copy(eOut, bufOut.Value.Coeffs[0][:len(eOut)])
}
