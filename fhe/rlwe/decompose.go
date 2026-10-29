package rlwe

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Decomposer gadget decomposes polynomial [Element] values.
type Decomposer struct {
	params          Parameters
	rnsDecomposer   *rnsDecomposer
	digitDecomposer *digitDecomposer
}

// NewDecomposer creates a new [Decomposer].
func NewDecomposer(params Parameters) *Decomposer {
	dcmp := &Decomposer{
		params: params,
	}

	switch params.gadgetParams.(type) {
	case RNSGadgetParameters:
		dcmp.rnsDecomposer = newRNSDecomposer(params)
	case DigitGadgetParameters:
		dcmp.digitDecomposer = newDigitDecomposer(params)
	}

	return dcmp
}

// Params returns the [Parameters] of the decomposer.
func (d *Decomposer) Params() Parameters {
	return d.params
}

// GadgetParams returns the [GadgetParameters] of the decomposer.
func (d *Decomposer) GadgetParams() GadgetParameters {
	return d.params.gadgetParams
}

// DecomposeLen returns the length of the decomposition of given base modulus length.
func (d *Decomposer) DecomposeLen(baseLen int) int {
	switch d.params.gadgetParams.(type) {
	case RNSGadgetParameters:
		return d.rnsDecomposer.decomposeLen(baseLen)
	case DigitGadgetParameters:
		return d.digitDecomposer.decomposeLen(baseLen)
	}
	return 0
}

// Decompose returns the decomposition of p.
func (d *Decomposer) Decompose(p *Element, outForm dft.Form) *Vector {
	pOut := NewVectorCustom(d.params.Rank(), p.BaseModLen(), d.params.AuxModLen(), d.DecomposeLen(p.BaseModLen()), p.Type(), outForm)
	d.DecomposeTo(pOut, p, outForm)
	return pOut
}

// DecomposeTo decomposes p to pOut.
func (d *Decomposer) DecomposeTo(pOut *Vector, p *Element, outForm dft.Form) {
	checkElementShape(d.params, pOut)
	checkElementShape(d.params, p)
	if p.AuxModLen() > 0 {
		panic("inputs must not have auxiliary modulus")
	} else if pOut.Len() != d.DecomposeLen(p.BaseModLen()) || pOut.AuxModLen() != d.params.AuxModLen() || pOut.BaseModLen() != p.BaseModLen() {
		panic("inconsistent inputs")
	}

	switch d.params.gadgetParams.(type) {
	case RNSGadgetParameters:
		d.rnsDecomposer.decomposeTo(pOut, p, outForm)
	case DigitGadgetParameters:
		d.digitDecomposer.decomposeTo(pOut, p, outForm)
	}
}

// rnsDecomposer computes RNS decomposition.
type rnsDecomposer struct {
	gadgetParams RNSGadgetParameters
	embedders    [][]*crt.Embedder
}

// newRNSDecomposer creates a new [rnsDecomposer].
func newRNSDecomposer(params Parameters) *rnsDecomposer {
	chunkSize := params.gadgetParams.(RNSGadgetParameters).chunkSize

	embPool := crt.NewEmbedderPool(params.Rank(), params.FullModLen())

	embedders := make([][]*crt.Embedder, params.GadgetLen())
	for i := range embedders {
		start := i * chunkSize
		end := min((i+1)*chunkSize, params.BaseModLen())
		embedders[i] = make([]*crt.Embedder, end-start)
		for j := range embedders[i] {
			modOp := params.Op().WithModSlice(params.AuxModLen()+start, params.AuxModLen()+start+j+1)
			embedders[i][j] = crt.NewEmbedder(params.Op(), modOp).WithPool(embPool)
		}
	}

	return &rnsDecomposer{
		gadgetParams: params.gadgetParams.(RNSGadgetParameters),
		embedders:    embedders,
	}
}

func (d *rnsDecomposer) decomposeLen(baseLen int) int {
	return num.DivCeil(baseLen, d.gadgetParams.chunkSize)
}

func (d *rnsDecomposer) decomposeTo(pOut *Vector, p *Element, outForm dft.Form) {
	for i := 0; i < pOut.Len(); i++ {
		start := i * d.gadgetParams.chunkSize
		end := min((i+1)*d.gadgetParams.chunkSize, p.BaseModLen())
		d.embedders[i][end-start-1].EmbedTo(pOut.Value[i].Value, p.Value.WithModSlice(start, end), outForm)
	}
}

// digitDecomposer computes digit decomposition.
type digitDecomposer struct {
	params       Parameters
	gadgetParams DigitGadgetParameters

	digitEmbedder []*crt.Embedder
	lastEmbedder  []*crt.Embedder
	modEmbedder   *crt.Embedder

	baseInv []uint64

	poolScalar *pool.Pool[*Element]
	poolPoly   *pool.Pool[*Element]
}

// newDigitDecomposer creates a new [digitDecomposer].
func newDigitDecomposer(params Parameters) *digitDecomposer {
	logDigitBase := params.gadgetParams.(DigitGadgetParameters).logDigitBase
	baseMod := num.NewModulus(uint64(1) << logDigitBase)

	embPool := crt.NewEmbedderPool(params.Rank(), params.FullModLen())

	digitEmbedder := make([]*crt.Embedder, params.BaseModLen())
	lastEmbedder := make([]*crt.Embedder, params.BaseModLen())
	baseModOp := crt.NewOperator(params.RingParams(), []*num.Modulus{baseMod})
	for i := range digitEmbedder {
		modOp := params.Op().WithModSlice(params.AuxModLen(), params.AuxModLen()+i+1)
		digitEmbedder[i] = crt.NewEmbedder(baseModOp, modOp).WithPool(embPool)
		lastEmbedder[i] = crt.NewEmbedder(params.Op(), modOp).WithPool(embPool)
	}
	modEmbedder := crt.NewEmbedder(params.Op(), baseModOp).WithPool(embPool)

	baseInv := make([]uint64, params.BaseModLen())
	for i := range baseInv {
		baseInv[i] = num.Inv(uint64(1)<<logDigitBase, params.baseMod[i])
	}

	return &digitDecomposer{
		params:       params,
		gadgetParams: params.gadgetParams.(DigitGadgetParameters),

		digitEmbedder: digitEmbedder,
		lastEmbedder:  lastEmbedder,
		modEmbedder:   modEmbedder,

		baseInv: baseInv,

		poolScalar: pool.NewPool(func() *Element {
			return NewScalar(params, false)
		}),
		poolPoly: pool.NewPool(func() *Element {
			return NewPoly(params, false, dft.FormCoeff)
		}),
	}
}

func (d *digitDecomposer) decomposeLen(baseLen int) int {
	modBitLen := 0
	for i := 0; i < baseLen; i++ {
		modBitLen += bits.Len64(d.params.baseMod[i].Value() - 1)
	}
	return num.DivCeil(modBitLen, d.gadgetParams.logDigitBase)
}

func (d *digitDecomposer) decomposeTo(pOut *Vector, p *Element, outForm dft.Form) {
	var pBufFull *Element
	switch {
	case p.Type() == crt.TypeScalar:
		pBufFull = d.poolScalar.Get()
		defer d.poolScalar.Put(pBufFull)
	case p.Type() == crt.TypePoly:
		pBufFull = d.poolPoly.Get()
		defer d.poolPoly.Put(pBufFull)
	}
	pBuf := pBufFull.WithModLen(p.BaseModLen(), 0)

	switch {
	case p.Form() == dft.FormCoeff:
		pBuf.CopyFrom(p)
	case p.Form() == dft.FormNTT:
		d.params.crtOp.WithModSlice(d.params.AuxModLen(), d.params.AuxModLen()+p.BaseModLen()).InvNTTTo(pBuf.Value, p.Value)
	}

	outOp := d.params.crtOp.WithModSlice(0, pOut.FullModLen())
	for i := 0; i < pOut.Len()-1; i++ {
		pOutDigit := pOut.Value[i].Value.WithModSlice(0, 1)
		d.digitEmbedder[pBuf.BaseModLen()-1].EmbedTo(pOutDigit, pBuf.Value, dft.FormCoeff)
		d.modEmbedder.EmbedTo(pOut.Value[i].Value, pOutDigit, dft.FormCoeff)
		for j := 0; j < pBuf.BaseModLen(); j++ {
			vec.SubTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], pOut.Value[i].Value.Coeffs[pOut.AuxModLen()+j], d.params.baseMod[j])
			vec.MulTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], d.baseInv[j], d.params.baseMod[j])
		}

		if p.Type() == crt.TypePoly && outForm == dft.FormNTT {
			outOp.FwdNTTTo(pOut.Value[i].Value, pOut.Value[i].Value)
		}
	}

	d.lastEmbedder[pBuf.BaseModLen()-1].EmbedTo(pOut.Value[pOut.Len()-1].Value, pBuf.Value, outForm)
}
