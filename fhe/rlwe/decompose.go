package rlwe

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Decomposer gadget decomposes [Poly].
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
func (d *Decomposer) DecomposeLen(p *Poly) int {
	switch d.params.gadgetParams.(type) {
	case RNSGadgetParameters:
		return d.rnsDecomposer.decomposeLen(p.BaseModLen())
	case DigitGadgetParameters:
		return d.digitDecomposer.decomposeLen(p.BaseModLen())
	}
	return 0
}

// Decompose decomposes p.
func (d *Decomposer) Decompose(p *Poly, outForm dft.Form) *Vector {
	pOut := NewVectorCustom(d.params.Rank(), p.BaseModLen(), len(d.params.auxMod), d.DecomposeLen(p), outForm)
	d.DecomposeTo(pOut, p, outForm)
	return pOut
}

// DecomposeTo decomposes p to pOut.
func (d *Decomposer) DecomposeTo(pOut *Vector, p *Poly, outForm dft.Form) {
	if p.Form() == dft.FormNTT {
		panic("input(s) must be in standard form")
	} else if p.AuxModLen() > 0 {
		panic("input(s) must not have auxiliary modulus")
	} else if pOut.Len() != d.DecomposeLen(p) || pOut.BaseModLen() != p.BaseModLen() || pOut.AuxModLen() != len(d.params.auxMod) {
		panic("input(s) not consistent")
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
	params    Parameters
	gadparams RNSGadgetParameters
	embedders [][]*crt.Embedder
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
		params:    params,
		gadparams: params.gadgetParams.(RNSGadgetParameters),
		embedders: embedders,
	}
}

func (d *rnsDecomposer) decomposeLen(baseLen int) int {
	return num.DivCeil(baseLen, d.gadparams.chunkSize)
}

func (d *rnsDecomposer) decomposeTo(pOut *Vector, p *Poly, outForm dft.Form) {
	for i := 0; i < pOut.Len(); i++ {
		start := i * d.gadparams.chunkSize
		end := min((i+1)*d.gadparams.chunkSize, p.BaseModLen())
		d.embedders[i][end-start-1].EmbedTo(pOut.Value[i].Value, p.Value.Slice(start, end), outForm)
	}
}

// digitDecomposer computes Digit Decomposition.
type digitDecomposer struct {
	params    Parameters
	gadParams DigitGadgetParameters

	digitEmbedder []*crt.Embedder
	lastEmbedder  []*crt.Embedder
	modEmbedder   *crt.Embedder

	pool *pool.Pool[*Poly]
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

	return &digitDecomposer{
		params:    params,
		gadParams: params.gadgetParams.(DigitGadgetParameters),

		digitEmbedder: digitEmbedder,
		lastEmbedder:  lastEmbedder,
		modEmbedder:   modEmbedder,

		pool: pool.NewPool(func() *Poly {
			return NewPoly(params, false, dft.FormCoeff)
		}),
	}
}

func (d *digitDecomposer) decomposeLen(baseLen int) int {
	modBitLen := 0
	for i := 0; i < baseLen; i++ {
		modBitLen += bits.Len64(d.params.baseMod[i].Value() - 1)
	}
	return num.DivCeil(modBitLen, d.gadParams.logDigitBase)
}

func (d *digitDecomposer) decomposeTo(pOut *Vector, p *Poly, outForm dft.Form) {
	pBufFull := d.pool.Get()
	defer d.pool.Put(pBufFull)

	pBuf := pBufFull.WithModLen(p.BaseModLen(), 0)
	pBuf.CopyFrom(p)

	base := uint64(1 << d.gadParams.logDigitBase)
	outOp := d.params.Operator().Slice(0, pOut.FullModLen())
	for i := 0; i < pOut.Len()-1; i++ {
		pOutDigit := pOut.Value[i].Value.WithModIdx(0)
		d.digitEmbedder[p.BaseModLen()-1].EmbedTo(pOutDigit, pBuf.Value, dft.FormCoeff)
		d.modEmbedder.EmbedTo(pOut.Value[i].Value, pOutDigit, dft.FormCoeff)
		for j := 0; j < p.BaseModLen(); j++ {
			vec.SubTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], pOut.Value[i].Value.Coeffs[pOut.AuxModLen()+j], d.params.baseMod[j])
			vec.MulTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], num.Inv(base, d.params.baseMod[j]), d.params.baseMod[j])
		}

		if outForm == dft.FormNTT {
			outOp.FwdNTTTo(pOut.Value[i].Value, pOut.Value[i].Value)
		}
	}

	d.lastEmbedder[p.BaseModLen()-1].EmbedTo(pOut.Value[pOut.Len()-1].Value, pBuf.Value, outForm)
}
