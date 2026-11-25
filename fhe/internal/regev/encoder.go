package regev

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/fhe/internal/rlwecheck"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Encoder encodes and decodes [rlwe.Element] values.
// Encode, EncodeCustom, and ScaleEncode take a uint64 for a scalar or a []uint64
// with one coefficient per ring dimension for a polynomial, even at rank one.
type Encoder struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	op *rlwe.PlainOperator

	embEcd    []*crt.Embedder
	embDcd    []*crt.Embedder
	embAuxDcd []*crt.Embedder

	scEcd []*crt.Scaler
	scDcd []*crt.Scaler

	msgEmbModEcd   *num.Modulus
	msgEmbModDcd   []*num.Modulus
	scale          []*rlwe.Element
	scaleRem       []uint64
	embModToExtDcd []*crt.VecEmbedder
	scExtToEmbDcd  []*crt.VecScaler
	embEmbToMsgDcd *crt.VecEmbedder

	poolEmb    *pool.Pool[*[]uint64]
	poolScalar *pool.Pool[*rlwe.Element]
	poolPoly   *pool.Pool[*rlwe.Element]
}

// NewEncoder creates a new [Encoder].
func NewEncoder(params rlwe.Parameters, msgMod *num.Modulus) *Encoder {
	op := rlwe.NewPlainOperator(params)
	msgOp := crt.NewOperator(params.RingParams(), []*num.Modulus{msgMod})

	poolEmb := crt.NewEmbedderPool(params.Rank(), params.FullModLen())
	embEcd := make([]*crt.Embedder, 1+params.AuxModLen())
	for i := range embEcd {
		modOp := params.Op().WithModSlice(params.AuxModLen()-i, params.FullModLen())
		embEcd[i] = crt.NewEmbedderUnsafe(modOp, msgOp, poolEmb)
	}

	embDcd := make([]*crt.Embedder, params.BaseModLen())
	embAuxDcd := make([]*crt.Embedder, 1+params.BaseModLen())
	for i := range embDcd {
		modOp := params.Op().WithModSlice(params.AuxModLen(), params.AuxModLen()+i+1)
		embDcd[i] = crt.NewEmbedderUnsafe(msgOp, modOp, poolEmb)
		modAuxOp := params.Op().WithModSlice(0, params.AuxModLen()+i+1)
		embAuxDcd[i] = crt.NewEmbedderUnsafe(msgOp, modAuxOp, poolEmb)
	}

	isOverlapping := false
	isCoprime := true
	for _, q := range params.BaseMod() {
		if msgMod.Value() == q.Value() {
			isOverlapping = true
		}
		if num.GCD(msgMod.Value(), q.Value()) != 1 {
			isCoprime = false
		}
	}

	var scEcd, scDcd []*crt.Scaler
	var msgEmbModEcd *num.Modulus
	var msgEmbModDcd []*num.Modulus
	var scale []*rlwe.Element
	var scaleRem []uint64
	var embModToExtDcd []*crt.VecEmbedder
	var scExtToEmbDcd []*crt.VecScaler
	var embEmbToMsgDcd *crt.VecEmbedder
	if isOverlapping || isCoprime {
		scEcd = make([]*crt.Scaler, params.BaseModLen())
		scDcd = make([]*crt.Scaler, params.BaseModLen())
		for i := 0; i < params.BaseModLen(); i++ {
			modOp := params.Op().WithModSlice(params.AuxModLen(), params.AuxModLen()+i+1)
			scEcd[i] = crt.NewScalerUnsafe(modOp, msgOp, poolEmb)
			scDcd[i] = crt.NewScalerUnsafe(msgOp, modOp, poolEmb)
		}
	} else {
		if msgMod.Value() == num.MaxModulus {
			panic("msgMod too large")
		}
		msgEmbModEcd = num.NewModulus(msgMod.Value() + 1)

		msgEmbModDcd = make([]*num.Modulus, 0)
		msgModLeft := msgMod.Value()
		for currMod := uint64(num.MaxModulus); currMod >= 2 && msgModLeft != 0; currMod-- {
			isCoprime := true
			for _, q := range params.BaseMod() {
				if num.GCD(currMod, q.Value()) != 1 {
					isCoprime = false
					break
				}
			}
			for _, q := range msgEmbModDcd {
				if num.GCD(currMod, q.Value()) != 1 {
					isCoprime = false
					break
				}
			}
			if !isCoprime {
				continue
			}

			msgEmbModDcd = append(msgEmbModDcd, num.NewModulus(currMod))
			msgModLeft /= currMod
		}

		scale = make([]*rlwe.Element, params.BaseModLen())
		scaleRem = make([]uint64, params.BaseModLen())
		remCurr := uint64(1)
		for i, q := range params.BaseMod() {
			remHi, remLo := bits.Mul64(remCurr, q.Value())
			remCarry, remNext := bits.Div64(remHi, remLo, msgMod.Value())

			scale[i] = rlwe.NewScalarCustom(i+1, 0)
			scale[i].Value.Coeffs[i][0] = remCarry
			scaleRem[i] = remNext
			for j, qj := range params.BaseMod()[:i] {
				scaleQuo := num.Mul(scale[i-1].Value.Coeffs[j][0], num.Reduce(q.Value(), qj), qj)
				scale[i].Value.Coeffs[j][0] = num.Add(scaleQuo, num.Reduce(remCarry, qj), qj)
			}
			remCurr = remNext
		}

		embModToExtDcd = make([]*crt.VecEmbedder, params.BaseModLen())
		scExtToEmbDcd = make([]*crt.VecScaler, params.BaseModLen())
		for i := 0; i < params.BaseModLen(); i++ {
			modExt := slices.Concat(params.BaseMod()[:i+1], msgEmbModDcd)
			embModToExtDcd[i] = crt.NewVecEmbedderUnsafe(modExt, params.BaseMod()[:i+1])
			scExtToEmbDcd[i] = crt.NewVecScalerUnsafe(msgEmbModDcd, modExt)
		}
		embEmbToMsgDcd = crt.NewVecEmbedderUnsafe([]*num.Modulus{msgMod}, msgEmbModDcd)
	}

	return &Encoder{
		params: params,
		msgMod: msgMod,

		op: op,

		embEcd:    embEcd,
		embDcd:    embDcd,
		embAuxDcd: embAuxDcd,

		scEcd: scEcd,
		scDcd: scDcd,

		msgEmbModEcd:   msgEmbModEcd,
		msgEmbModDcd:   msgEmbModDcd,
		scale:          scale,
		scaleRem:       scaleRem,
		embModToExtDcd: embModToExtDcd,
		scExtToEmbDcd:  scExtToEmbDcd,
		embEmbToMsgDcd: embEmbToMsgDcd,

		poolEmb: poolEmb,
		poolScalar: pool.NewPool(func() *rlwe.Element {
			baseLen := params.BaseModLen()
			if msgEmbModDcd != nil {
				baseLen += len(msgEmbModDcd)
			}
			return rlwe.NewScalarCustom(baseLen, params.AuxModLen())
		}),
		poolPoly: pool.NewPool(func() *rlwe.Element {
			baseLen := params.BaseModLen()
			if msgEmbModDcd != nil {
				baseLen += len(msgEmbModDcd)
			}
			return rlwe.NewPolyCustom(params.Rank(), baseLen, params.AuxModLen(), dft.FormCoeff)
		}),
	}
}

func (ecd *Encoder) encodingInput[T uint64 | []uint64](v T) (crt.ElementType, []uint64) {
	if x, ok := any(v).(uint64); ok {
		return crt.TypeScalar, []uint64{x}
	}
	coeffs := any(v).([]uint64)
	if len(coeffs) != ecd.params.Rank() {
		panic("inconsistent inputs")
	}
	return crt.TypePoly, coeffs
}

// Encode returns the encoding of v.
func (ecd *Encoder) Encode[T uint64 | []uint64](v T, withAux bool, outForm dft.Form) *rlwe.Element {
	return ecd.EncodeCustom(v, ecd.params.BaseModLen(), withAux, outForm)
}

// EncodeCustom returns the encoding of v with custom parameters.
func (ecd *Encoder) EncodeCustom[T uint64 | []uint64](v T, baseLen int, withAux bool, outForm dft.Form) *rlwe.Element {
	v64, vIsUint64 := any(v).(uint64)
	vVec, _ := any(v).([]uint64)
	if vIsUint64 {
		vVec = []uint64{v64}
	}

	if withAux && ecd.params.AuxModLen() == 0 {
		panic("no auxiliary modulus")
	}

	auxLen := ecd.params.AuxModLen()
	if !withAux {
		auxLen = 0
	}

	var eOut *rlwe.Element
	switch {
	case vIsUint64:
		eOut = rlwe.NewScalarCustom(baseLen, auxLen)
	case !vIsUint64:
		eOut = rlwe.NewPolyCustom(ecd.params.Rank(), baseLen, auxLen, dft.FormCoeff)
	}

	ecd.EncodeTo(eOut, vVec, outForm)
	return eOut
}

// EncodeTo encodes v to eOut.
func (ecd *Encoder) EncodeTo(eOut *rlwe.Element, v []uint64, outForm dft.Form) {
	rlwecheck.CheckElementShape(ecd.params, eOut)
	if len(v) != eOut.Rank() {
		panic("inconsistent inputs")
	}

	eBase := eOut.WithModLen(1, 0)
	vec.ReduceTo(eBase.Value.Coeffs[0], v, ecd.msgMod)
	eBase.SetForm(dft.FormCoeff)

	ecd.embEcd[eOut.AuxModLen()].EmbedTo(eOut.Value, eBase.Value, outForm)
}

// ScaleEncode returns the scaled encoding of v at the full base modulus, without auxiliary limbs.
func (ecd *Encoder) ScaleEncode[T uint64 | []uint64](v T, outForm dft.Form) *rlwe.Element {
	eOutType, coeffs := ecd.encodingInput(v)
	eOut := rlwe.NewElement(ecd.params, false, eOutType, dft.FormCoeff)
	ecd.ScaleEncodeTo(eOut, coeffs, outForm)
	return eOut
}

// ScaleEncodeTo scales and encodes v to eOut.
func (ecd *Encoder) ScaleEncodeTo(eOut *rlwe.Element, v []uint64, outForm dft.Form) {
	rlwecheck.CheckElementShape(ecd.params, eOut)
	if eOut.AuxModLen() > 0 || len(v) != eOut.Rank() {
		panic("inconsistent inputs")
	}

	if ecd.scEcd != nil {
		eBase := eOut.WithModLen(1, 0)
		vec.ReduceTo(eBase.Value.Coeffs[0], v, ecd.msgMod)
		eBase.SetForm(dft.FormCoeff)

		ecd.scEcd[eOut.BaseModLen()-1].ScaleTo(eOut.Value, eBase.Value, outForm)
		return
	}

	bias := (ecd.msgMod.Value() - 1) >> 1

	vBufPtr := ecd.poolEmb.Get()
	vBuf := (*vBufPtr)[:eOut.Rank()]
	defer ecd.poolEmb.Put(vBufPtr)

	vec.ReduceTo(vBuf, v, ecd.msgMod)

	remMsgModPtr := ecd.poolEmb.Get()
	remMsgMod := (*remMsgModPtr)[:eOut.Rank()]
	defer ecd.poolEmb.Put(remMsgModPtr)

	remMsgModEmbPtr := ecd.poolEmb.Get()
	remMsgModEmb := (*remMsgModEmbPtr)[:eOut.Rank()]
	defer ecd.poolEmb.Put(remMsgModEmbPtr)

	vec.MulTo(remMsgMod, vBuf, ecd.scaleRem[eOut.BaseModLen()-1], ecd.msgMod)
	vec.AddTo(remMsgMod, remMsgMod, bias, ecd.msgMod)

	vec.MulTo(remMsgModEmb, vBuf, ecd.scaleRem[eOut.BaseModLen()-1], ecd.msgEmbModEcd)
	vec.AddTo(remMsgModEmb, remMsgModEmb, bias, ecd.msgEmbModEcd)

	vec.SubTo(remMsgMod, remMsgMod, remMsgModEmb, ecd.msgEmbModEcd)

	for i, q := range ecd.params.BaseMod()[:eOut.BaseModLen()] {
		vec.ReduceTo(remMsgModEmb, vBuf, q)
		vec.ReduceTo(eOut.Value.Coeffs[i], remMsgMod, q)
		vec.MulAddTo(eOut.Value.Coeffs[i], remMsgModEmb, ecd.scale[eOut.BaseModLen()-1].Value.Coeffs[i][0], q)
	}

	if eOut.Type() == crt.TypePoly {
		eOut.SetForm(dft.FormCoeff)
		if outForm == dft.FormNTT {
			ecd.op.FwdNTTTo(eOut, eOut)
		}
	}
}

// Decode returns the decoding of e.
func (ecd *Encoder) Decode(e *rlwe.Element) []uint64 {
	vOut := make([]uint64, e.Rank())
	ecd.DecodeTo(vOut, e)
	return vOut
}

// DecodeTo decodes e to vOut.
func (ecd *Encoder) DecodeTo(vOut []uint64, e *rlwe.Element) {
	rlwecheck.CheckElementShape(ecd.params, e)
	if e.Rank() != len(vOut) {
		panic("inconsistent inputs")
	}

	var eBufFull *rlwe.Element
	switch {
	case e.Type() == crt.TypeScalar:
		eBufFull = ecd.poolScalar.Get()
		defer ecd.poolScalar.Put(eBufFull)
	case e.Type() == crt.TypePoly:
		eBufFull = ecd.poolPoly.Get()
		defer ecd.poolPoly.Put(eBufFull)
	}
	eBase := eBufFull.WithModLen(1, 0)

	if e.AuxModLen() > 0 {
		ecd.embAuxDcd[e.BaseModLen()-1].EmbedTo(eBase.Value, e.Value, dft.FormCoeff)
	} else {
		ecd.embDcd[e.BaseModLen()-1].EmbedTo(eBase.Value, e.Value, dft.FormCoeff)
	}
	copy(vOut, eBase.Value.Coeffs[0])
}

// ScaleDecode returns the scaled decoding of e at its active base modulus.
func (ecd *Encoder) ScaleDecode(e *rlwe.Element) []uint64 {
	vOut := make([]uint64, e.Rank())
	ecd.ScaleDecodeTo(vOut, e)
	return vOut
}

// ScaleDecodeTo scales and decodes e to vOut.
func (ecd *Encoder) ScaleDecodeTo(vOut []uint64, e *rlwe.Element) {
	rlwecheck.CheckElementShape(ecd.params, e)
	if e.AuxModLen() > 0 || e.Rank() != len(vOut) {
		panic("inconsistent inputs")
	}

	var eInvNTTFull *rlwe.Element
	switch {
	case e.Type() == crt.TypeScalar:
		eInvNTTFull = ecd.poolScalar.Get()
		defer ecd.poolScalar.Put(eInvNTTFull)
	case e.Type() == crt.TypePoly:
		eInvNTTFull = ecd.poolPoly.Get()
		defer ecd.poolPoly.Put(eInvNTTFull)
	}

	eInvNTT := eInvNTTFull.WithModLen(e.BaseModLen(), 0)
	eInvNTT.CopyFrom(e)
	if e.Type() == crt.TypePoly && e.Form() == dft.FormNTT {
		ecd.op.InvNTTTo(eInvNTT, eInvNTT)
	}

	if ecd.scDcd != nil {
		eInvNTTBase := eInvNTTFull.WithModLen(1, 0)
		eInvNTTBase.SetForm(dft.FormCoeff)

		ecd.scDcd[e.BaseModLen()-1].ScaleTo(eInvNTTBase.Value, eInvNTT.Value, dft.FormCoeff)
		copy(vOut, eInvNTTBase.Value.Coeffs[0])
		return
	}

	eInvNTTExt := eInvNTTFull.WithModLen(e.BaseModLen()+len(ecd.msgEmbModDcd), 0)
	eInvNTTExt.SetForm(dft.FormCoeff)

	ecd.embModToExtDcd[e.BaseModLen()-1].EmbedTo(eInvNTTExt.Value.Coeffs, eInvNTT.Value.Coeffs)
	for i, q := range ecd.params.BaseMod()[:e.BaseModLen()] {
		vec.MulTo(eInvNTTExt.Value.Coeffs[i], eInvNTTExt.Value.Coeffs[i], num.Reduce(ecd.msgMod.Value(), q), q)
	}
	for i, q := range ecd.msgEmbModDcd {
		vec.MulTo(eInvNTTExt.Value.Coeffs[e.BaseModLen()+i], eInvNTTExt.Value.Coeffs[e.BaseModLen()+i], num.Reduce(ecd.msgMod.Value(), q), q)
	}

	eInvNTTEmb := eInvNTTFull.WithModLen(len(ecd.msgEmbModDcd), 0)
	eInvNTTEmb.SetForm(dft.FormCoeff)

	ecd.scExtToEmbDcd[e.BaseModLen()-1].ScaleTo(eInvNTTEmb.Value.Coeffs, eInvNTTExt.Value.Coeffs)
	ecd.embEmbToMsgDcd.EmbedTo([][]uint64{vOut}, eInvNTTEmb.Value.Coeffs)
}
