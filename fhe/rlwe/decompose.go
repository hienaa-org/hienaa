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

// DecomposeLen returns the length of the decomposition of p.
func (d *Decomposer) DecomposeLen(p *Element) int {
	checkPolyShape(d.params, p)
	switch d.params.gadgetParams.(type) {
	case RNSGadgetParameters:
		return d.rnsDecomposer.decomposeLen(p.BaseModLen())
	case DigitGadgetParameters:
		return d.digitDecomposer.decomposeLen(p.BaseModLen())
	}
	return 0
}

// Decompose returns the decomposition of p.
//
// Panics when p is in NTT form.
func (d *Decomposer) Decompose(p *Element, outForm dft.Form) *Vector {
	pOut := NewVectorCustom(d.params.Rank(), p.BaseModLen(), len(d.params.auxMod), d.DecomposeLen(p), outForm)
	d.DecomposeTo(pOut, p, outForm)
	return pOut
}

// DecomposeTo decomposes p to pOut.
//
// Panics when p is in NTT form.
func (d *Decomposer) DecomposeTo(pOut *Vector, p *Element, outForm dft.Form) {
	checkPolyShape(d.params, p)
	if p.Form() == dft.FormNTT {
		panic("inputs must be in standard form")
	} else if p.AuxModLen() > 0 {
		panic("inputs must not have auxiliary modulus")
	} else if pOut.Len() != d.DecomposeLen(p) || pOut.AuxModLen() != len(d.params.auxMod) || pOut.BaseModLen() != p.BaseModLen() {
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

	embPool := crt.NewEmbedderPool(params.Rank(), len(params.fullMod))

	embedders := make([][]*crt.Embedder, params.GadgetLen())
	for i := range embedders {
		start := i * chunkSize
		end := min((i+1)*chunkSize, len(params.baseMod))
		embedders[i] = make([]*crt.Embedder, end-start)
		for j := range embedders[i] {
			modOp := params.Operator().Slice(len(params.auxMod)+start, len(params.auxMod)+start+j+1)
			embedders[i][j] = crt.NewEmbedder(params.Operator(), modOp).WithPool(embPool)
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
		d.embedders[i][end-start-1].EmbedTo(pOut.Value[i].Value, p.Value.Slice(start, end), outForm)
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

	pool *pool.Pool[*Element]
}

// newDigitDecomposer creates a new [digitDecomposer].
func newDigitDecomposer(params Parameters) *digitDecomposer {
	logDigitBase := params.gadgetParams.(DigitGadgetParameters).logDigitBase
	baseMod := num.NewModulus(uint64(1) << logDigitBase)

	embPool := crt.NewEmbedderPool(params.Rank(), len(params.fullMod))

	digitEmbedder := make([]*crt.Embedder, len(params.baseMod))
	lastEmbedder := make([]*crt.Embedder, len(params.baseMod))
	baseModOp := crt.NewOperator(params.RingParams(), []*num.Modulus{baseMod})
	for i := range digitEmbedder {
		modOp := params.Operator().Slice(len(params.auxMod), len(params.auxMod)+i+1)
		digitEmbedder[i] = crt.NewEmbedder(baseModOp, modOp).WithPool(embPool)
		lastEmbedder[i] = crt.NewEmbedder(params.Operator(), modOp).WithPool(embPool)
	}
	modEmbedder := crt.NewEmbedder(params.Operator(), baseModOp).WithPool(embPool)

	baseInv := make([]uint64, len(params.baseMod))
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

		pool: pool.NewPool(func() *Element {
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
	pBufFull := d.pool.Get()
	defer d.pool.Put(pBufFull)

	pBuf := pBufFull.WithModLen(p.BaseModLen(), 0)
	pBuf.CopyFrom(p)

	outOp := d.params.Operator().Slice(0, pOut.FullModLen())
	for i := 0; i < pOut.Len()-1; i++ {
		pOutDigit := pOut.Value[i].Value.Slice(0, 1)
		d.digitEmbedder[p.BaseModLen()-1].EmbedTo(pOutDigit, pBuf.Value, dft.FormCoeff)
		d.modEmbedder.EmbedTo(pOut.Value[i].Value, pOutDigit, dft.FormCoeff)
		for j := 0; j < p.BaseModLen(); j++ {
			vec.SubTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], pOut.Value[i].Value.Coeffs[pOut.AuxModLen()+j], d.params.baseMod[j])
			vec.MulTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], d.baseInv[j], d.params.baseMod[j])
		}

		if outForm == dft.FormNTT {
			outOp.FwdNTTTo(pOut.Value[i].Value, pOut.Value[i].Value)
		}
	}

	d.lastEmbedder[p.BaseModLen()-1].EmbedTo(pOut.Value[pOut.Len()-1].Value, pBuf.Value, outForm)
}
