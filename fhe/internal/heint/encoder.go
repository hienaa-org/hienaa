package heint

import (
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Encoder encodes/decodes []uint64 into/from [*rlwe.Element].
type Encoder struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	pOp       *rlwe.PlainOperator
	msgToFull []*crt.Embedder
	baseToMsg []*crt.Embedder

	ePool *rlwe.ElementPool
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
		params: rlweParams,
		msgMod: msgMod,

		pOp:       pOp,
		msgToFull: msgToFull,
		baseToMsg: baseToMsg,

		ePool: rlwe.NewElementPool(rlweParams, true, true),
	}
}

// Encode encodes a []uint64 into a [*rlwe.Element].
func (ecd *Encoder) Encode(eIn []uint64, hasAux, isNTT bool) *rlwe.Element {
	auxLen := len(ecd.params.AuxModulus())
	if !hasAux {
		auxLen = 0
	}

	eOut := rlwe.NewElement(len(eIn), len(ecd.params.BaseModulus()), auxLen, isNTT)
	ecd.EncodeTo(eOut, eIn, isNTT)
	return eOut
}

// EncodeCustom encodes a []uint64 into a [*rlwe.Element] with custom parameters.
func (ecd *Encoder) EncodeCustom(eIn []uint64, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	eOut := rlwe.NewElement(len(eIn), baseLen, auxLen, isNTT)
	ecd.EncodeTo(eOut, eIn, isNTT)
	return eOut
}

// EncodeTo encodes a []uint64 into a [*rlwe.Element].
func (ecd *Encoder) EncodeTo(eOut *rlwe.Element, eIn []uint64, isNTT bool) {
	if len(eIn) != eOut.Rank() {
		panic("inconsistent input(s)")
	}

	if eOut.Rank() != 1 && eOut.Rank() != ecd.params.Rank() {
		panic("invalid output length")
	}

	eBase := eOut.WithModLen(1, 0)
	eBase.Value.IsNTT = false
	vec.ReduceTo(eBase.Value.Coeffs[0], eIn, ecd.msgMod)

	auxLen := eOut.AuxModLen()
	ecd.msgToFull[auxLen].EmbedTo(eOut.Value, eBase.Value)
	eOut.Value.IsNTT = false

	if isNTT && eOut.Value.Type() == crt.TypePoly {
		ecd.pOp.FwdNTTTo(eOut, eOut)
	}
}

// Decode decodes a [*rlwe.Element] into a []uint64.
func (ecd *Encoder) Decode(e *rlwe.Element) []uint64 {
	eOut := make([]uint64, e.Rank())
	ecd.DecodeTo(eOut, e)
	return eOut
}

// DecodeTo decodes a [*rlwe.Element] into a []uint64.
func (ecd *Encoder) DecodeTo(eOut []uint64, e *rlwe.Element) {
	if e.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	} else if len(eOut) != e.Rank() {
		panic("invalid output length")
	}

	baseLen := e.BaseModLen()

	buf := ecd.ePool.Get(e.Type())
	defer ecd.ePool.Put(buf)
	buf = buf.WithModLen(baseLen, 0)
	bufOut := buf.WithModLen(1, 0)

	if e.IsNTT() {
		ecd.pOp.InvNTTTo(buf, e)
		ecd.baseToMsg[baseLen-1].EmbedTo(bufOut.Value, buf.Value)
	} else {
		ecd.baseToMsg[baseLen-1].EmbedTo(bufOut.Value, e.Value)
	}

	copy(eOut, bufOut.Value.Coeffs[0][:len(eOut)])
}
