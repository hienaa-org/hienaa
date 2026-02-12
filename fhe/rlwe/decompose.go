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
	GadgetVector() []*crt.Element
	// DecomposeLen returns the length of the decomposition for a given modulus length.
	DecomposeLen(modLen int) int
	// Decompose decomposes p.
	Decompose(p *crt.Element) *Tensor
	// DecomposeTo decomposes p into pOut.
	DecomposeTo(pOut *Tensor, p *crt.Element)
}

// NewDecomposer creates a new [Decomposer] for the given parameters.
func NewDecomposer(p Parameters) Decomposer {
	switch p.gadgetParams.(type) {
	case RNSGadgetParameters:
		return newRNSDecomposer(p)
	case DigitGadgetParameters:
		return newDigitDecomposer(p)
	}
	panic("unsupported parameters")
}

// rnsDecomposer computes RNS decomposition.
type rnsDecomposer struct {
	params    Parameters
	gadparams RNSGadgetParameters

	gadVec    []*crt.Element
	embedders [][]*crt.Embedder
}

// newRNSDecomposer creates a new [rnsDecomposer].
func newRNSDecomposer(p Parameters) Decomposer {
	chunkSize := p.gadgetParams.(RNSGadgetParameters).chunkSize
	gadLen := p.GadgetLen()

	gadVec := make([]*crt.Element, gadLen)
	embedders := make([][]*crt.Embedder, gadLen)

	auxModBig := big.NewInt(1)
	for _, q := range p.auxModulus {
		auxModBig.Mul(auxModBig, new(big.Int).SetUint64(q.Value()))
	}

	for i := 0; i < gadLen; i++ {
		start := i * chunkSize
		end := min((i+1)*chunkSize, len(p.modulus))

		embedders[i] = make([]*crt.Embedder, end-start)
		for j := range embedders[i] {
			embedders[i][j] = crt.NewEmbedder(p.fullModulus, p.modulus[start:start+j+1])
		}

		g := big.NewInt(1)
		chunk := big.NewInt(1)
		for j, q := range p.modulus {
			if start <= j && j < end {
				chunk.Mul(chunk, new(big.Int).SetUint64(q.Value()))
			} else {
				g.Mul(g, new(big.Int).SetUint64(q.Value()))
			}
		}
		g.Mul(g, new(big.Int).ModInverse(g, chunk))
		g.Mul(g, auxModBig)

		gadVec[i] = crt.NewScalar(g, p.fullModulus)
	}

	return &rnsDecomposer{
		params:    p,
		gadparams: p.gadgetParams.(RNSGadgetParameters),

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
func (d *rnsDecomposer) GadgetVector() []*crt.Element {
	return d.gadVec
}

// DecomposeLen returns the length of the decomposition.
func (d *rnsDecomposer) DecomposeLen(modLen int) int {
	return int(math.Ceil(float64(modLen) / float64(d.gadparams.chunkSize)))
}

// Decompose decomposes a polynomial into a tensor using the RNS gadget.
func (d *rnsDecomposer) Decompose(p *crt.Element) *Tensor {
	modLen := p.ModLen()
	auxLen := len(d.params.auxModulus)
	resLen := d.DecomposeLen(modLen)

	pOut := NewTensorCustom(p.Rank(), modLen, auxLen, resLen, false)
	d.DecomposeTo(pOut, p)

	return pOut
}

// DecomposeTo decomposes p into pOut using the RNS gadget.
func (d *rnsDecomposer) DecomposeTo(pOut *Tensor, p *crt.Element) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	}

	modLen := p.ModLen()
	auxLen := len(d.params.auxModulus)
	resLen := d.DecomposeLen(modLen)

	if (pOut.Degree() != resLen) || (pOut.ModLen() != (modLen + auxLen)) {
		panic("inconsistent input(s)")
	}

	chunkSize := d.gadparams.chunkSize
	for i := 0; i < resLen; i++ {
		start := i * chunkSize
		end := min((i+1)*chunkSize, modLen)

		d.embedders[i][end-start-1].EmbedVecTo(pOut.Value[i].Coeffs, p.Coeffs[start:end])
	}
}

// digitDecomposer computes digit decomposition.
type digitDecomposer struct {
	params    Parameters
	gadParams DigitGadgetParameters

	gadVec []*crt.Element

	digitEmbedder []*crt.Embedder
	modEmbedder   *crt.Embedder

	pool *sync.Pool
}

// newDigitDecomposer creates a new [digitDecomposer].
func newDigitDecomposer(params Parameters) Decomposer {
	logDigitBase := params.gadgetParams.(DigitGadgetParameters).logDigitBase
	baseMod := num.NewModulus(1 << logDigitBase)

	gadLen := params.GadgetLen()
	gadVec := make([]*crt.Element, gadLen)

	g := big.NewInt(1)
	for _, q := range params.auxModulus {
		g.Mul(g, new(big.Int).SetUint64(q.Value()))
	}

	for i := 0; i < gadLen; i++ {
		gadVec[i] = crt.NewScalar(g, params.fullModulus)
		g.Lsh(g, uint(logDigitBase))
	}

	digitEmbedder := make([]*crt.Embedder, len(params.modulus))
	for i := range params.modulus {
		digitEmbedder[i] = crt.NewEmbedder([]*num.Modulus{baseMod}, params.modulus[:i+1])
	}
	modEmbedder := crt.NewEmbedder(params.fullModulus, []*num.Modulus{baseMod})

	return &digitDecomposer{
		params:    params,
		gadParams: params.gadgetParams.(DigitGadgetParameters),

		gadVec: gadVec,

		digitEmbedder: digitEmbedder,
		modEmbedder:   modEmbedder,

		pool: &sync.Pool{
			New: func() any {
				return crt.NewPoly(params.ringParams.Rank(), len(params.modulus))
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
func (d *digitDecomposer) GadgetVector() []*crt.Element {
	return d.gadVec
}

// DecomposeLen returns the length of the decomposition.
func (d *digitDecomposer) DecomposeLen(modLen int) int {
	modBitLen := 0.0
	for i := 0; i < modLen; i++ {
		modBitLen += num.Log2(d.params.modulus[i].Value())
	}
	return int(math.Ceil(modBitLen / float64(d.gadParams.logDigitBase)))
}

// Decompose outputs the decomposition of p using the digit gadget.
func (d *digitDecomposer) Decompose(p *crt.Element) *Tensor {
	modLen := p.ModLen()
	auxLen := len(d.params.auxModulus)
	decmpLen := d.DecomposeLen(modLen)

	pOut := NewTensorCustom(p.Rank(), modLen, auxLen, decmpLen, false)
	d.DecomposeTo(pOut, p)

	return pOut
}

// DecomposeTo decomposes p into pOut using the digit gadget.
func (d *digitDecomposer) DecomposeTo(pOut *Tensor, p *crt.Element) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	}

	modLen := p.ModLen()
	auxLen := len(d.params.auxModulus)
	decmpLen := d.DecomposeLen(modLen)

	if (pOut.Degree() != decmpLen) || (pOut.ModLen() != (modLen + auxLen)) {
		panic("inconsistent input(s)")
	}

	pBuf := d.pool.Get().(*crt.Element)
	defer d.pool.Put(pBuf)

	pBuf = pBuf.WithModIdx(vec.Range(0, modLen)...)
	pBuf.CopyFrom(p)

	base := uint64(1 << d.gadParams.logDigitBase)
	for i := 0; i < decmpLen; i++ {
		d.digitEmbedder[modLen-1].EmbedVecTo(pOut.Value[i].Coeffs[:1], pBuf.Coeffs)
		d.modEmbedder.EmbedVecTo(pOut.Value[i].Coeffs, pOut.Value[i].Coeffs[:1])

		for j := 0; j < modLen; j++ {
			vec.SubTo(pBuf.Coeffs[j], pBuf.Coeffs[j], pOut.Value[i].Coeffs[auxLen+j], d.params.modulus[j])
			vec.MulScalarTo(pBuf.Coeffs[j], pBuf.Coeffs[j], num.Inv(base, d.params.modulus[j]), d.params.modulus[j])
		}
	}
}
