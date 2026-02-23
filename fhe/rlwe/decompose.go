package rlwe

import (
	"math"
	"math/big"
	"sync"

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
	// DecomposeLen returns the length of the decomposition for a given modulus length.
	DecomposeLen(modLen int) int
	// Decompose decomposes p.
	Decompose(p *Element) *Vector
	// DecomposeTo decomposes p into pOut.
	DecomposeTo(pOut *Vector, p *Element)
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

	for i := 0; i < gadLen; i++ {
		start := i * chunkSize
		end := min((i+1)*chunkSize, len(params.baseMod))

		embedders[i] = make([]*crt.Embedder, end-start)
		for j := range embedders[i] {
			embedders[i][j] = crt.NewEmbedder(params.fullMod, params.baseMod[start:start+j+1])
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
			hasAux: params.HasAuxModulus(),
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
func (d *rnsDecomposer) DecomposeLen(modLen int) int {
	return int(math.Ceil(float64(modLen) / float64(d.gadparams.chunkSize)))
}

// Decompose decomposes a polynomial into a tensor using the RNS gadget.
func (d *rnsDecomposer) Decompose(p *Element) *Vector {
	pOut := NewVectorCustom(
		p.Rank(), p.ModLen()+len(d.params.auxMod), d.params.HasAuxModulus(), d.DecomposeLen(p.ModLen()), false,
	)
	d.DecomposeTo(pOut, p)
	return pOut
}

// DecomposeTo decomposes p into pOut using the RNS gadget.
func (d *rnsDecomposer) DecomposeTo(pOut *Vector, p *Element) {
	if p.IsNTT() {
		panic("input(s) must be in standard form")
	} else if p.HasAuxModulus() {
		panic("input(s) must not have auxiliary modulus")
	} else if pOut.Len() != d.DecomposeLen(p.ModLen()) || pOut.ModLen() != p.ModLen()+len(d.params.auxMod) {
		panic("output not consistent")
	}

	for i := 0; i < pOut.Len(); i++ {
		start := i * d.gadparams.chunkSize
		end := min((i+1)*d.gadparams.chunkSize, p.ModLen())
		d.embedders[i][end-start-1].EmbedTo(pOut.Value[i].Value, p.Value.WithModIdx(vec.Range(start, end)...))
	}
}

// digitDecomposer computes digit decomposition.
type digitDecomposer struct {
	params    Parameters
	gadParams DigitGadgetParameters

	gadVec []*Element

	digitEmbedder []*crt.Embedder
	modEmbedder   *crt.Embedder

	pool *sync.Pool
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
			hasAux: params.HasAuxModulus(),
		}
		g.Lsh(g, uint(logDigitBase))
	}

	digitEmbedder := make([]*crt.Embedder, len(params.baseMod))
	for i := range params.baseMod {
		digitEmbedder[i] = crt.NewEmbedder([]*num.Modulus{baseMod}, params.baseMod[:i+1])
	}
	modEmbedder := crt.NewEmbedder(params.fullMod, []*num.Modulus{baseMod})

	return &digitDecomposer{
		params:    params,
		gadParams: params.gadgetParams.(DigitGadgetParameters),

		gadVec: gadVec,

		digitEmbedder: digitEmbedder,
		modEmbedder:   modEmbedder,

		pool: &sync.Pool{
			New: func() any {
				return NewPoly(params, false, false)
			},
		},
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
func (d *digitDecomposer) DecomposeLen(modLen int) int {
	modBitLen := 0.0
	for i := 0; i < modLen; i++ {
		modBitLen += num.Log2(d.params.baseMod[i].Value())
	}
	return int(math.Ceil(modBitLen / float64(d.gadParams.logDigitBase)))
}

// Decompose outputs the decomposition of p using the digit gadget.
func (d *digitDecomposer) Decompose(p *Element) *Vector {
	pOut := NewVectorCustom(
		p.Rank(), p.ModLen()+len(d.params.auxMod), d.params.HasAuxModulus(), d.DecomposeLen(p.ModLen()), false,
	)
	d.DecomposeTo(pOut, p)
	return pOut
}

// DecomposeTo decomposes p into pOut using the digit gadget.
func (d *digitDecomposer) DecomposeTo(pOut *Vector, p *Element) {
	modLen, auxLen, dcmpLen := p.ModLen(), len(d.params.auxMod), d.DecomposeLen(p.ModLen())

	if p.IsNTT() {
		panic("input(s) must be in standard form")
	} else if p.HasAuxModulus() {
		panic("input(s) must not have auxiliary modulus")
	} else if pOut.Len() != dcmpLen || pOut.ModLen() != modLen+auxLen {
		panic("output not consistent")
	}

	pBuf := d.pool.Get().(*Element)
	defer d.pool.Put(pBuf)

	pBuf = pBuf.WithModIdx(vec.Range(0, p.ModLen())...)
	pBuf.CopyFrom(p)

	base := uint64(1 << d.gadParams.logDigitBase)
	for i := 0; i < dcmpLen; i++ {
		d.digitEmbedder[modLen-1].EmbedTo(pOut.Value[i].Value.WithModIdx(0), pBuf.Value)
		d.modEmbedder.EmbedTo(pOut.Value[i].Value, pOut.Value[i].Value.WithModIdx(0))

		for j := 0; j < modLen; j++ {
			vec.SubTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], pOut.Value[i].Value.Coeffs[auxLen+j], d.params.baseMod[j])
			vec.MulScalarTo(pBuf.Value.Coeffs[j], pBuf.Value.Coeffs[j], num.Inv(base, d.params.baseMod[j]), d.params.baseMod[j])
		}
	}
}
