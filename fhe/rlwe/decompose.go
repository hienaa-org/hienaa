package rlwe

import (
	"math"
	"math/big"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// GadgetType is the type of the gadget.
type GadgetType uint64

const (
	// RNS is the RNS gadget type.
	RNS GadgetType = iota
	// Digit is the digit gadget type.
	Digit
)

// GadgetParameters is the parameters for the gadget.
type GadgetParameters interface {
	// gadgetLen returns the length of the gadget.
	gadgetLen(p Parameters) int
}

// RNSGadgetParameters is the parameters for the RNS gadget.
type RNSGadgetParameters struct {
	// ChunkSize is the size of the chunk.
	ChunkSize int
	// Type is the type of the gadget.
	Type GadgetType
}

// NewRNSGadgetParameters creates a new [GadgetParameters] for the RNS gadget.
func NewRNSGadgetParameters(chunkSize int) RNSGadgetParameters {
	if chunkSize <= 0 {
		panic("NewRNSGadgetParameters: chunk size must be positive")
	}

	return RNSGadgetParameters{
		ChunkSize: chunkSize,
		Type:      RNS,
	}
}

// gadgetLen returns the length of the gadget.
func (g RNSGadgetParameters) gadgetLen(p Parameters) int {
	return int(math.Ceil(float64(len(p.modulus)) / float64(g.ChunkSize)))
}

type DigitGadgetParameters struct {
	// LogDigitBase is the logarithm of the base of the digit decomposition.
	LogDigitBase int
	// Type is the type of the gadget.
	Type GadgetType
}

// NewDigitGadgetParameters creates a new [GadgetParameters] for the digit gadget.
func NewDigitGadgetParameters(logDigitBase int) DigitGadgetParameters {
	if logDigitBase <= 0 || logDigitBase > num.MaxModulusBits {
		panic("NewDigitGadgetParameters: log digit base must be between 1 and max modulus bits")
	}

	return DigitGadgetParameters{
		LogDigitBase: logDigitBase,
		Type:         Digit,
	}
}

// gadgetLen returns the length of the gadget.
func (g DigitGadgetParameters) gadgetLen(p Parameters) int {
	modBitLen := 0.0
	for _, m := range p.modulus {
		modBitLen += math.Log2(float64(m.Value()))
	}
	return int(math.Ceil(modBitLen / float64(g.LogDigitBase)))
}

// maxGadgetLen returns the maximum length of the gadget.
func maxGadgetLen(p Parameters) int {
	switch p.gadgetParams.(type) {
	case RNSGadgetParameters:
		return int(math.Ceil(float64(len(p.modulus)) / float64(p.gadgetParams.(RNSGadgetParameters).ChunkSize)))
	case DigitGadgetParameters:
		modBitLen := 0.0
		for _, m := range p.modulus {
			modBitLen += math.Log2(float64(m.Value()))
		}
		return int(math.Ceil(modBitLen / float64(p.gadgetParams.(DigitGadgetParameters).LogDigitBase)))
	default:
		panic("maxGadgetLen: invalid gadget type")
	}
}

// Decomposer is a struct that decomposes a polynomial into a tensor.
type Decomposer interface {
	GadgetParams() GadgetParameters
	GadgetVector() []crt.Scalar
	DecomposeLen(modLen int) int
	Decompose(p *crt.Poly) *Tensor
	DecomposeTo(pOut *Tensor, p *crt.Poly)
	SafeCopy() Decomposer
}

// NewDecomposer creates a new [Decomposer] for the given parameters.
func NewDecomposer(p Parameters) Decomposer {
	switch p.gadgetParams.(type) {
	case RNSGadgetParameters:
		return newRNSDecomposer(p)
	case DigitGadgetParameters:
		return newDigitDecomposer(p)
	default:
		panic("NewDecomposer: invalid gadget type")
	}
}

// rnsDecomposer is a struct that decomposes a polynomial into a tensor using the RNS gadget.
type rnsDecomposer struct {
	schemeParams Parameters
	gadgetParams RNSGadgetParameters

	gadgetVector []crt.Scalar
	embedders    [][]*crt.Embedder
}

// newRNSDecomposer creates a new [rnsDecomposer] for the given parameters.
func newRNSDecomposer(p Parameters) Decomposer {
	chunkSize := p.gadgetParams.(RNSGadgetParameters).ChunkSize
	gadgetLen := p.gadgetParams.gadgetLen(p)
	gadgetMod := append(p.auxModulus, p.modulus...)

	gadgetVector := make([]crt.Scalar, gadgetLen)
	embedders := make([][]*crt.Embedder, gadgetLen)

	auxModBig := big.NewInt(1)
	for _, modi := range p.auxModulus {
		auxModBig.Mul(auxModBig, big.NewInt(int64(modi.Value())))
	}

	var gadgeti, chunki *big.Int
	for i := 0; i < gadgetLen; i++ {
		lo := i * chunkSize
		hi := int(math.Min(float64((i+1)*chunkSize), float64(len(p.modulus))))

		embedders[i] = make([]*crt.Embedder, hi-lo)
		for j := 0; j < hi-lo; j++ {
			embedders[i][j] = crt.NewEmbedder(gadgetMod, p.modulus[lo:lo+j+1])
		}

		gadgeti = big.NewInt(1)
		chunki = big.NewInt(1)
		for j, modj := range p.modulus {
			if lo <= j && j < hi {
				chunki.Mul(chunki, big.NewInt(int64(modj.Value())))
			} else {
				gadgeti.Mul(gadgeti, big.NewInt(int64(modj.Value())))
			}
		}
		gadgeti.Mul(gadgeti, big.NewInt(0).ModInverse(gadgeti, chunki))
		gadgeti.Mul(gadgeti, auxModBig)

		gadgetVector[i] = crt.NewScalar(gadgeti, gadgetMod)
	}

	return &rnsDecomposer{
		schemeParams: p,
		gadgetParams: p.gadgetParams.(RNSGadgetParameters),

		gadgetVector: gadgetVector,
		embedders:    embedders,
	}
}

// Params returns the parameters of the decomposer.
func (d *rnsDecomposer) GadgetParams() GadgetParameters {
	return d.gadgetParams
}

// GadgetVector returns the gadget vector.
func (d *rnsDecomposer) GadgetVector() []crt.Scalar {
	return d.gadgetVector
}

// DecomposeLen returns the length of the decomposition.
func (d *rnsDecomposer) DecomposeLen(modLen int) int {
	return int(math.Ceil(float64(modLen) / float64(d.gadgetParams.ChunkSize)))
}

// Decompose decomposes a polynomial into a tensor using the RNS gadget.
func (d *rnsDecomposer) Decompose(p *crt.Poly) *Tensor {
	modLen := p.ModLen()
	auxLen := len(d.schemeParams.auxModulus)
	resLen := d.DecomposeLen(modLen)

	pOut := NewTensorCustom(p.Rank(), modLen, auxLen, resLen, false)
	d.DecomposeTo(pOut, p)

	return pOut
}

// DecomposeTo decomposes p into pOut using the RNS gadget.
func (d *rnsDecomposer) DecomposeTo(pOut *Tensor, p *crt.Poly) {
	if p.IsNTT {
		panic("rnsDecomposer: polynomial must be in the coefficient form")
	}

	modLen := p.ModLen()
	auxLen := len(d.schemeParams.auxModulus)
	resLen := d.DecomposeLen(modLen)

	if pOut.Degree() != resLen {
		panic("rnsDecomposer: inconsistent tensor degree")
	}

	if pOut.ModLen() != (modLen + auxLen) {
		panic("rnsDecomposer: inconsistent tensor modulus length")
	}

	chunkSize := d.gadgetParams.ChunkSize
	for i := 0; i < resLen; i++ {
		lo := i * chunkSize
		hi := int(math.Min(float64((i+1)*chunkSize), float64(modLen)))

		d.embedders[i][hi-lo-1].EmbedVecTo(pOut.Value[i].Coeffs, p.Coeffs[lo:hi])
	}
}

// SafeCopy returns a thread-safe copy.
func (d *rnsDecomposer) SafeCopy() Decomposer {
	embedders := make([][]*crt.Embedder, len(d.embedders))
	for i := range d.embedders {
		embedders[i] = make([]*crt.Embedder, len(d.embedders[i]))
		for j := range d.embedders[i] {
			embedders[i][j] = d.embedders[i][j].SafeCopy()
		}
	}

	return &rnsDecomposer{
		schemeParams: d.schemeParams,
		gadgetParams: d.gadgetParams,

		gadgetVector: d.gadgetVector,
		embedders:    embedders,
	}
}

// digitDecomposer is a struct that decomposes a polynomial into a tensor using the digit gadget.
type digitDecomposer struct {
	schemeParams Parameters
	gadgetParams DigitGadgetParameters

	bufPoly      *crt.Poly
	gadgetVector []crt.Scalar

	digitEmbedder []*crt.Embedder
	modEmbedder   *crt.Embedder
}

// GadgetParams returns the gadget parameters.
func (d *digitDecomposer) GadgetParams() GadgetParameters {
	return d.gadgetParams
}

// newDigitDecomposer creates a new [digitDecomposer] for the given parameters.
func newDigitDecomposer(p Parameters) Decomposer {
	logDigitBase := p.gadgetParams.(DigitGadgetParameters).LogDigitBase
	baseMod := num.NewModulus(1 << logDigitBase)

	gadgetLen := p.gadgetParams.gadgetLen(p)
	gadgetVector := make([]crt.Scalar, gadgetLen)

	gadgetMod := append(p.auxModulus, p.modulus...)
	gadgeti := big.NewInt(1)
	for _, modi := range p.auxModulus {
		gadgeti.Mul(gadgeti, big.NewInt(int64(modi.Value())))
	}

	for i := 0; i < gadgetLen; i++ {
		gadgetVector[i] = crt.NewScalar(gadgeti, gadgetMod)
		gadgeti.Lsh(gadgeti, uint(logDigitBase))
	}

	digitEmbedder := make([]*crt.Embedder, len(p.modulus))
	for i := range p.modulus {
		digitEmbedder[i] = crt.NewEmbedder([]*num.Modulus{baseMod}, p.modulus[:i+1])
	}
	modEmbedder := crt.NewEmbedder(gadgetMod, []*num.Modulus{baseMod})

	return &digitDecomposer{
		schemeParams: p,
		gadgetParams: p.gadgetParams.(DigitGadgetParameters),

		bufPoly:      crt.NewPoly(p.ringParams.Rank(), len(p.modulus)),
		gadgetVector: gadgetVector,

		digitEmbedder: digitEmbedder,
		modEmbedder:   modEmbedder,
	}
}

// Params returns the parameters of the decomposer.
func (d *digitDecomposer) Params() GadgetParameters {
	return d.gadgetParams
}

// GadgetVector returns the gadget vector.
func (d *digitDecomposer) GadgetVector() []crt.Scalar {
	return d.gadgetVector
}

// DecomposeLen returns the length of the decomposition.
func (d *digitDecomposer) DecomposeLen(modLen int) int {
	modBitLen := 0.0
	for i := 0; i < modLen; i++ {
		modBitLen += math.Log2(float64(d.schemeParams.modulus[i].Value()))
	}
	return int(math.Ceil(modBitLen / float64(d.gadgetParams.LogDigitBase)))
}

// Decompose outputs the decomposition of p using the digit gadget.
func (d *digitDecomposer) Decompose(p *crt.Poly) *Tensor {
	modLen := p.ModLen()
	auxLen := len(d.schemeParams.auxModulus)
	decmpLen := d.DecomposeLen(modLen)

	pOut := NewTensorCustom(p.Rank(), modLen, auxLen, decmpLen, false)
	d.DecomposeTo(pOut, p)

	return pOut
}

// DecomposeTo decomposes p into pOut using the digit gadget.
func (d *digitDecomposer) DecomposeTo(pOut *Tensor, p *crt.Poly) {
	if p.IsNTT {
		panic("digitDecomposer: polynomial must be in the coefficient form")
	}

	modLen := p.ModLen()
	auxLen := len(d.schemeParams.auxModulus)
	decmpLen := d.DecomposeLen(modLen)

	if pOut.Degree() != decmpLen {
		panic("digitDecomposer: inconsistent tensor degree")
	}

	if pOut.ModLen() != (modLen + auxLen) {
		panic("digitDecomposer: inconsistent tensor modulus length")
	}

	buf := &crt.Poly{
		Coeffs: d.bufPoly.Coeffs[:modLen],
	}
	buf.CopyFrom(p)

	gadgetBase := uint64(1 << d.gadgetParams.LogDigitBase)
	for i := 0; i < decmpLen; i++ {
		d.digitEmbedder[modLen-1].EmbedVecTo(pOut.Value[i].Coeffs[:1], buf.Coeffs)
		d.modEmbedder.EmbedVecTo(pOut.Value[i].Coeffs, pOut.Value[i].Coeffs[:1])

		for j := 0; j < modLen; j++ {
			vec.SubTo(buf.Coeffs[j], buf.Coeffs[j], pOut.Value[i].Coeffs[auxLen+j], d.schemeParams.modulus[j])
			vec.ScalarMulTo(buf.Coeffs[j], buf.Coeffs[j], num.Inv(gadgetBase, d.schemeParams.modulus[j]), d.schemeParams.modulus[j])
		}
	}
}

// SafeCopy returns a thread-safe copy.
func (d *digitDecomposer) SafeCopy() Decomposer {
	digitEmbedder := make([]*crt.Embedder, len(d.digitEmbedder))
	for i := range d.digitEmbedder {
		digitEmbedder[i] = d.digitEmbedder[i].SafeCopy()
	}

	return &digitDecomposer{
		schemeParams: d.schemeParams,
		gadgetParams: d.gadgetParams,

		bufPoly:      crt.NewPoly(d.schemeParams.ringParams.Rank(), len(d.schemeParams.modulus)),
		gadgetVector: d.gadgetVector,

		digitEmbedder: digitEmbedder,
		modEmbedder:   d.modEmbedder.SafeCopy(),
	}
}
