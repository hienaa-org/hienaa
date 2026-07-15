package rlwe

import (
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Decomposer is a struct that decomposes a polynomial into a tensor.
type Decomposer interface {
	// Params returns the parameters.
	Params() Parameters
	// GadgetParams returns the gadget parameters.
	GadgetParams() GadgetParameters
	// GadgetVector returns the gadget vector.
	GadgetVector() []*Element
	// DecomposeLen returns the length of the decomposition for a given base modulus length.
	DecomposeLen(baseLen int) int
	// AuxModLen returns the auxiliary modulus length for the given base modulus length.
	AuxModLen(baseLen int) int
	// Decompose decomposes p.
	Decompose(p *Element, isNTT bool) *Vector
	// DecomposeTo decomposes p into pOut.
	DecomposeTo(pOut *Vector, p *Element, isNTT bool)
}

// NewDecomposer creates a new [Decomposer] for the given parameters.
func NewDecomposer(params Parameters) Decomposer {
	switch params.gadgetParams.(type) {
	case RNSGadgetParameters:
		return newRNSDecomposer(params)
	case DigitGadgetParameters:
		return newDigitDecomposer(params)
	}
	panic("unsupported parameters")
}

// rnsDecomposer computes RNS decomposition.
type rnsDecomposer struct {
	params    Parameters
	gadparams RNSGadgetParameters

	gadVec    []*Element
	embedders [][]*crt.Embedder
}

// newRNSDecomposer creates a new [rnsDecomposer].
func newRNSDecomposer(params Parameters) Decomposer {
	chunkSize := params.gadgetParams.(RNSGadgetParameters).chunkSize
	gadLen := params.GadgetLen()

	gadVec := make([]*Element, gadLen)
	embedders := make([][]*crt.Embedder, gadLen)

	auxModBig := big.NewInt(1)
	for _, q := range params.auxMod {
		auxModBig.Mul(auxModBig, new(big.Int).SetUint64(q.Value()))
	}

	auxLen := len(params.auxMod)
	embPool := pool.NewPool(func() *[]uint64 {
		v := make([]uint64, params.Rank())
		return &v
	})
	for i := 0; i < gadLen; i++ {
		start := i * chunkSize
		end := min((i+1)*chunkSize, len(params.baseMod))

		embedders[i] = make([]*crt.Embedder, end-start)
		fullOp := params.Operator()
		for j := range embedders[i] {
			modOp := fullOp.Slice(auxLen+start, auxLen+start+j+1)
			embedders[i][j] = crt.NewEmbedder(fullOp, modOp).WithPool(embPool)
		}

		g := big.NewInt(1)
		chunk := big.NewInt(1)
		for j, q := range params.baseMod {
			if start <= j && j < end {
				chunk.Mul(chunk, new(big.Int).SetUint64(q.Value()))
			} else {
				g.Mul(g, new(big.Int).SetUint64(q.Value()))
			}
		}
		g.Mul(g, new(big.Int).ModInverse(g, chunk))
		g.Mul(g, auxModBig)

		gadVec[i] = &Element{
			Value:  crt.NewScalarFrom(g, params.fullMod),
			auxLen: len(params.auxMod),
		}
	}

	return &rnsDecomposer{
		params:    params,
		gadparams: params.gadgetParams.(RNSGadgetParameters),

		gadVec:    gadVec,
		embedders: embedders,
	}
}

// Params returns the parameters of the decomposer.
func (d *rnsDecomposer) Params() Parameters {
	return d.params
}

// GadgetParams returns the gadget parameters.
func (d *rnsDecomposer) GadgetParams() GadgetParameters {
	return d.gadparams
}

// GadgetVector returns the gadget vector.
func (d *rnsDecomposer) GadgetVector() []*Element {
	return d.gadVec
}

// DecomposeLen returns the length of the decomposition.
func (d *rnsDecomposer) DecomposeLen(baseLen int) int {
	return int(math.Ceil(float64(baseLen) / float64(d.gadparams.chunkSize)))
}

// TODO: Consider the edge cases where the auxiliary modulus can be smaller.
// AuxModLen returns the auxiliary modulus length for the given base modulus length.
func (d *rnsDecomposer) AuxModLen(baseLen int) int {
	return len(d.params.auxMod)
}

// Decompose decomposes a polynomial into a tensor using the RNS gadget.
func (d *rnsDecomposer) Decompose(p *Element, isNTT bool) *Vector {
	auxLen := d.AuxModLen(p.BaseModLen())

	pOut := NewVectorCustom(
		p.Rank(), p.BaseModLen(), auxLen, d.DecomposeLen(p.BaseModLen()), false,
	)
	d.DecomposeTo(pOut, p, isNTT)
	return pOut
}

// DecomposeTo decomposes p into pOut using the RNS gadget.
func (d *rnsDecomposer) DecomposeTo(pOut *Vector, p *Element, isNTT bool) {
	if p.IsNTT() {
		panic("input(s) must be in standard form")
	} else if p.AuxModLen() > 0 {
		panic("input(s) must not have auxiliary modulus")
	} else if pOut.Len() != d.DecomposeLen(p.BaseModLen()) || pOut.BaseModLen() != p.BaseModLen() || pOut.AuxModLen() != d.AuxModLen(p.BaseModLen()) {
		panic("output not consistent")
	}

	for i := 0; i < pOut.Len(); i++ {
		start := i * d.gadparams.chunkSize
		end := min((i+1)*d.gadparams.chunkSize, p.BaseModLen())
		d.embedders[i][end-start-1].EmbedTo(pOut.Value[i].Value, p.Value.Slice(start, end), isNTT)
	}
}

// digitDecomposer computes digit decomposition.
type digitDecomposer struct {
	params    Parameters
	gadParams DigitGadgetParameters

	gadVec []*Element

	digitEmbedder []*crt.Embedder
	modEmbedder   *crt.Embedder

	pool *pool.Pool[*Element]
}

// newDigitDecomposer creates a new [digitDecomposer].
func newDigitDecomposer(params Parameters) Decomposer {
	logDigitBase := params.gadgetParams.(DigitGadgetParameters).logDigitBase
	baseMod := num.NewModulus(1 << logDigitBase)

	gadLen := params.GadgetLen()
	gadVec := make([]*Element, gadLen)

	g := big.NewInt(1)
	for _, q := range params.auxMod {
		g.Mul(g, new(big.Int).SetUint64(q.Value()))
	}

	for i := 0; i < gadLen; i++ {
		gadVec[i] = &Element{
			Value:  crt.NewScalarFrom(g, params.fullMod),
			auxLen: len(params.auxMod),
		}
		g.Lsh(g, uint(logDigitBase))
	}

	digitEmbedder := make([]*crt.Embedder, len(params.baseMod))
	baseModOp := crt.NewOperator(params.RingParams(), []*num.Modulus{baseMod})
	auxLen := len(params.auxMod)
	embPool := pool.NewPool(func() *[]uint64 {
		v := make([]uint64, params.Rank())
		return &v
	})
	for i := range params.baseMod {
		modOp := params.Operator().Slice(auxLen, auxLen+i+1)
		digitEmbedder[i] = crt.NewEmbedder(baseModOp, modOp).WithPool(embPool)
	}
	modEmbedder := crt.NewEmbedder(params.Operator(), baseModOp).WithPool(embPool)

	return &digitDecomposer{
		params:    params,
		gadParams: params.gadgetParams.(DigitGadgetParameters),

		gadVec: gadVec,

		digitEmbedder: digitEmbedder,
		modEmbedder:   modEmbedder,

		pool: pool.NewPool(func() *Element {
			return NewPoly(params, false, false)
		}),
	}
}

// Params returns the parameters of the decomposer.
func (d *digitDecomposer) Params() Parameters {
	return d.params
}

// GadgetParams returns the gadget parameters.
func (d *digitDecomposer) GadgetParams() GadgetParameters {
	return d.gadParams
}

// GadgetVector returns the gadget vector.
func (d *digitDecomposer) GadgetVector() []*Element {
	return d.gadVec
}

// DecomposeLen returns the length of the decomposition.
func (d *digitDecomposer) DecomposeLen(baseLen int) int {
	modBitLen := 0.0
	for i := 0; i < baseLen; i++ {
		modBitLen += num.Log2(d.params.baseMod[i].Value())
	}
	return int(math.Ceil(modBitLen / float64(d.gadParams.logDigitBase)))
}

// AuxModLen returns the auxiliary modulus length for the given base modulus length.
func (d *digitDecomposer) AuxModLen(baseLen int) int {
	return len(d.params.auxMod)
}

// Decompose outputs the decomposition of p using the digit gadget.
func (d *digitDecomposer) Decompose(p *Element, isNTT bool) *Vector {
	auxLen := d.AuxModLen(p.BaseModLen())
	pOut := NewVectorCustom(
		p.Rank(), p.BaseModLen(), auxLen, d.DecomposeLen(p.BaseModLen()), false,
	)
	d.DecomposeTo(pOut, p, isNTT)
	return pOut
}

// DecomposeTo decomposes p into pOut using the digit gadget.
func (d *digitDecomposer) DecomposeTo(pOut *Vector, p *Element, isNTT bool) {
	baseLen, auxLen, dcmpLen := p.BaseModLen(), d.AuxModLen(p.BaseModLen()), d.DecomposeLen(p.BaseModLen())

	if p.IsNTT() {
		panic("input(s) must be in standard form")
	} else if p.auxLen > 0 {
		panic("input(s) must not have auxiliary modulus")
	} else if pOut.Len() != dcmpLen || pOut.BaseModLen() != baseLen || pOut.AuxModLen() != auxLen {
		panic("output not consistent")
	}

	pBuf := d.pool.Get()
	defer d.pool.Put(pBuf)

	pBuf = pBuf.WithModLen(baseLen, 0)
	pBuf.CopyFrom(p)

	base := uint64(1 << d.gadParams.logDigitBase)
	for i := 0; i < dcmpLen; i++ {
		d.digitEmbedder[baseLen-1].EmbedTo(pOut.Value[i].Value.WithModIdx(0), pBuf.Value, false)
		d.modEmbedder.EmbedTo(pOut.Value[i].Value, pOut.Value[i].Value.WithModIdx(0), isNTT)

		for j := 0; j < baseLen; j++ {
			vec.SubTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], pOut.Value[i].Value.Coeffs[auxLen+j], d.params.baseMod[j])
			vec.MulScalarTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], num.Inv(base, d.params.baseMod[j]), d.params.baseMod[j])
		}
	}
}
