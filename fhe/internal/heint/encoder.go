package heint

import (
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Encoder encodes/decodes []uint64 into/from [*rlwe.Element].
type Encoder struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	pOp       *rlwe.PlainOperator
	embEncode []*crt.Embedder
	embDecode []*crt.Embedder
	scEncode  []*crt.Scaler
	scDecode  []*crt.Scaler

	ePool *rlwe.ElementPool
}

// NewEncoder creates a new [Encoder].
func NewEncoder(rlweParams rlwe.Parameters, msgMod *num.Modulus) *Encoder {
	pOp := rlwe.NewPlainOperator(rlweParams)
	msgOp := crt.NewOperator(rlweParams.RingParams(), []*num.Modulus{msgMod})

	embEncode := make([]*crt.Embedder, 1+len(rlweParams.AuxModulus()))
	auxLen := len(rlweParams.AuxModulus())
	modLen := len(rlweParams.BaseModulus())
	embPool := pool.NewPool(func() *[]uint64 {
		v := make([]uint64, rlweParams.Rank())
		return &v
	})
	for i := range embEncode {
		modOp := rlweParams.Operator().WithModIdx(vec.Range(auxLen-i, auxLen+modLen)...)
		embEncode[i] = crt.NewEmbedder(modOp, msgOp).WithPool(embPool)
	}

	embDecode := make([]*crt.Embedder, len(rlweParams.BaseModulus()))
	for i := range embDecode {
		modOp := rlweParams.Operator().WithModIdx(vec.Range(auxLen, auxLen+i+1)...)
		embDecode[i] = crt.NewEmbedder(msgOp, modOp).WithPool(embPool)
	}

	scEncode := make([]*crt.Scaler, len(rlweParams.BaseModulus()))
	scDecode := make([]*crt.Scaler, len(rlweParams.BaseModulus()))
	for i := range scEncode {
		modOp := rlweParams.Operator().WithModIdx(vec.Range(auxLen, auxLen+i+1)...)
		scEncode[i] = crt.NewScaler(modOp, msgOp).WithPool(embPool)
		scDecode[i] = crt.NewScaler(msgOp, modOp).WithPool(embPool)
	}

	return &Encoder{
		params: rlweParams,
		msgMod: msgMod,

		pOp:       pOp,
		embEncode: embEncode,
		embDecode: embDecode,
		scEncode:  scEncode,
		scDecode:  scDecode,

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
	ecd.embEncode[auxLen].EmbedTo(eOut.Value, eBase.Value, isNTT)
}

// ScaleEncode encodes a []uint64 into a [*rlwe.Element] while scaling.
func (ecd *Encoder) ScaleEncode(eIn []uint64, isNTT bool) *rlwe.Element {
	eOut := rlwe.NewElement(len(eIn), len(ecd.params.BaseModulus()), 0, isNTT)
	ecd.ScaleEncodeTo(eOut, eIn, isNTT)
	return eOut
}

// ScaleEncodeCustom encodes a []uint64 into a [*rlwe.Element] with custom parameters while scaling.
func (ecd *Encoder) ScaleEncodeCustom(eIn []uint64, baseLen int, isNTT bool) *rlwe.Element {
	eOut := rlwe.NewElement(len(eIn), baseLen, 0, isNTT)
	ecd.ScaleEncodeTo(eOut, eIn, isNTT)
	return eOut
}

// ScaleEncodeTo encodes a []uint64 into a [*rlwe.Element] while scaling.
func (ecd *Encoder) ScaleEncodeTo(eOut *rlwe.Element, eIn []uint64, isNTT bool) {
	if len(eIn) != eOut.Rank() {
		panic("inconsistent input(s)")
	}

	if eOut.Rank() != 1 && eOut.Rank() != ecd.params.Rank() {
		panic("invalid output length")
	}

	if eOut.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	}

	eBase := eOut.WithModLen(1, 0)
	eBase.Value.IsNTT = false
	vec.ReduceTo(eBase.Value.Coeffs[0], eIn, ecd.msgMod)

	baseLen := eOut.BaseModLen()
	ecd.scEncode[baseLen-1].ScaleTo(eOut.Value, eBase.Value, isNTT)
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
		ecd.embDecode[baseLen-1].EmbedTo(bufOut.Value, buf.Value, false)
	} else {
		ecd.embDecode[baseLen-1].EmbedTo(bufOut.Value, e.Value, false)
	}

	copy(eOut, bufOut.Value.Coeffs[0][:len(eOut)])
}

// ScaleDecode decodes a [*rlwe.Element] into a []uint64 while scaling.
func (ecd *Encoder) ScaleDecode(e *rlwe.Element) []uint64 {
	eOut := make([]uint64, e.Rank())
	ecd.ScaleDecodeTo(eOut, e)
	return eOut
}

// ScaleDecodeTo decodes a [*rlwe.Element] into a []uint64 while scaling.
func (ecd *Encoder) ScaleDecodeTo(eOut []uint64, e *rlwe.Element) {
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
		ecd.scDecode[baseLen-1].ScaleTo(bufOut.Value, buf.Value, false)
	} else {
		ecd.scDecode[baseLen-1].ScaleTo(bufOut.Value, e.Value, false)
	}

	copy(eOut, bufOut.Value.Coeffs[0][:len(eOut)])
}
