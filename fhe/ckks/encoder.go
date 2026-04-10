package ckks

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
)

// Encoder encodes/decodes []complex128 into/from [*rlwe.Element].
type Encoder struct {
	params rlwe.Parameters

	pOp *rlwe.PlainOperator

	fPool   *pool.Pool[*big.Float]
	bigPool *pool.Pool[*big.Int]
	vPool   *pool.Pool[[]*big.Int]
	ePool   *rlwe.ElementPool
}

// NewEncoder creates a new [Encoder].
func NewEncoder(params rlwe.Parameters) *Encoder {
	pOp := rlwe.NewPlainOperator(params)

	return &Encoder{
		params: params,
		pOp:    pOp,
		fPool: pool.NewPool(func() *big.Float {
			return big.NewFloat(0).SetPrec(52 + 52 + 64)
		}),
		bigPool: pool.NewPool(func() *big.Int {
			return new(big.Int)
		}),
		vPool: pool.NewPool(func() []*big.Int {
			res := make([]*big.Int, params.Rank())
			for i := range res {
				res[i] = new(big.Int)
			}
			return res
		}),
		ePool: rlwe.NewElementPool(params, true, true),
	}
}

// Encode encodes a []float64 into a [*rlwe.Element].
func (ecd *Encoder) Encode(eIn *Plaintext, scFac float64, hasAux, isNTT bool) *rlwe.Element {
	auxLen := len(ecd.params.AuxModulus())
	if !hasAux {
		auxLen = 0
	}

	eOut := rlwe.NewElement(eIn.Rank(), len(ecd.params.BaseModulus()), auxLen, isNTT)
	ecd.EncodeTo(eOut, eIn, scFac, isNTT)
	return eOut
}

// EncodeCustom encodes a []float64 into a [*rlwe.Element] with custom parameters.
func (ecd *Encoder) EncodeCustom(eIn *Plaintext, scFac float64, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	eOut := rlwe.NewElement(eIn.Rank(), baseLen, auxLen, isNTT)
	ecd.EncodeTo(eOut, eIn, scFac, isNTT)
	return eOut
}

// EncodeTo encodes a []float64 into a [*rlwe.Element].
func (ecd *Encoder) EncodeTo(eOut *rlwe.Element, eIn *Plaintext, scFac float64, isNTT bool) {
	if eIn.Rank() != eOut.Rank() {
		panic("inconsistent input(s)")
	}

	if eOut.Rank() != 1 && eOut.Rank() != ecd.params.Rank() {
		panic("invalid output length")
	}

	tmpFloat := ecd.fPool.Get()
	scBig := ecd.fPool.Get()
	defer ecd.fPool.Put(tmpFloat)
	defer ecd.fPool.Put(scBig)

	v := ecd.vPool.Get()
	defer ecd.vPool.Put(v)
	v = v[:eOut.Rank()]

	scBig.SetFloat64(scFac)
	for i := range v {
		tmpFloat.SetFloat64(eIn.Value[i])
		tmpFloat.Mul(tmpFloat, scBig)
		tmpFloat.Add(tmpFloat, big.NewFloat(0.5))
		tmpFloat.Int(v[i])
	}

	tmpInt := ecd.bigPool.Get()
	modBig := ecd.bigPool.Get()
	zero := ecd.bigPool.Get()
	defer ecd.bigPool.Put(tmpInt)
	defer ecd.bigPool.Put(modBig)
	defer ecd.bigPool.Put(zero)

	zero.SetUint64(0)
	auxLen := len(ecd.params.AuxModulus())
	for i := 0; i < eOut.BaseModLen(); i++ {
		modBig.SetUint64(ecd.params.BaseModulus()[i].Value())
		for j := range v {
			if v[j].Cmp(zero) >= 0 {
				tmpInt.Rem(v[j], modBig)
				eOut.Value.Coeffs[eOut.AuxModLen()+i][j] = tmpInt.Uint64()
			} else {
				tmpInt.Neg(v[j])
				tmpInt.Rem(tmpInt, modBig)
				eOut.Value.Coeffs[eOut.AuxModLen()+i][j] = ecd.params.BaseModulus()[i].Value() - tmpInt.Uint64()
			}
		}
	}
	for i := 0; i < eOut.AuxModLen(); i++ {
		modBig.SetUint64(ecd.params.AuxModulus()[auxLen-eOut.AuxModLen()+i].Value())
		for j := range v {
			if v[j].Cmp(zero) >= 0 {
				tmpInt.Rem(v[j], modBig)
				eOut.Value.Coeffs[i][j] = tmpInt.Uint64()
			} else {
				tmpInt.Neg(v[j])
				tmpInt.Rem(tmpInt, modBig)
				eOut.Value.Coeffs[i][j] = ecd.params.AuxModulus()[auxLen-eOut.AuxModLen()+i].Value() - tmpInt.Uint64()
			}
		}
	}

	eOut.Value.IsNTT = false
	if isNTT && eOut.Value.Type() == crt.TypePoly {
		ecd.pOp.FwdNTTTo(eOut, eOut)
	}
}

// Decode decodes a [*rlwe.Element] into a []float64.
func (ecd *Encoder) Decode(e *rlwe.Element, scFac float64) *Plaintext {
	eOut := NewPoly(e.Rank(), TypeReal)
	ecd.DecodeTo(eOut, e, scFac)
	return eOut
}

// DecodeTo decodes a [*rlwe.Element] into a [*Plaintext].
func (ecd *Encoder) DecodeTo(eOut *Plaintext, e *rlwe.Element, scFac float64) {
	if e.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	} else if eOut.Rank() != e.Rank() {
		panic("invalid output length")
	}

	buf := ecd.ePool.Get(e.Type())
	defer ecd.ePool.Put(buf)
	buf = buf.WithModLen(e.BaseModLen(), 0)

	// TODO: refrain from heavy bigInt vector allocation.
	var out []*big.Int
	if e.IsNTT() {
		ecd.pOp.InvNTTTo(buf, e)
		out = ecd.pOp.AsBig(buf)
	} else {
		out = ecd.pOp.AsBig(e)
	}

	for i := range eOut.Value {
		val, _ := out[i].Float64()
		eOut.Value[i] = val / scFac
	}
}
