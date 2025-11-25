package rlwe

import (
	"math"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// GadgetType is the type of the gadget.
type GadgetType uint64

const (
	// RNS is the RNS gadget type.
	RNS GadgetType = iota
	// Digit is the digit gadget type.
	Digit
)

// GadgetParameters is the parameters for gadget encryption.
type GadgetParameters struct {
	// Length is the length of the gadget vector.
	Length int
	// Type is the type of the gadget.
	Type GadgetType
}

// NewRNSGadgetParameters creates a new [GadgetParameters] for the RNS gadget.
func NewRNSGadgetParameters(length int) GadgetParameters {
	return GadgetParameters{
		Length: length,
		Type:   RNS,
	}
}

// NewDigitGadgetParameters creates a new [GadgetParameters] for the digit gadget.
func NewDigitGadgetParameters(length int) GadgetParameters {
	return GadgetParameters{
		Length: length,
		Type:   Digit,
	}
}

// Decomposer is a struct that decomposes a polynomial into a tensor.
type Decomposer interface {
	Params() GadgetParameters
	Decompose(p *crt.Poly) *Tensor
}

// NewDecomposer creates a new [Decomposer] for the given parameters.
func NewDecomposer(g GadgetParameters, p Parameters) Decomposer {
	switch g.Type {
	case RNS:
		return newRNSDecomposer(g, p)
	case Digit:
		return newDigitDecomposer(g, p)
	default:
		panic("NewDecomposer: invalid gadget type")
	}
}

// rnsDecomposer is a struct that decomposes a polynomial into a tensor using the RNS gadget.
type rnsDecomposer struct {
	params    GadgetParameters
	embedders []*crt.Embedder
}

// newRNSDecomposer creates a new [rnsDecomposer] for the given parameters.
func newRNSDecomposer(g GadgetParameters, p Parameters) Decomposer {
	dlen := int(math.Ceil(float64(len(p.modulus)) / float64(g.Length)))
	embedders := make([]*crt.Embedder, dlen)
	for i := 0; i < dlen-1; i++ {
		embedders[i] = crt.NewEmbedder(p.modulus[i*g.Length:(i+1)*g.Length], p.modulus)
	}
	embedders[dlen-1] = crt.NewEmbedder(p.modulus[(dlen-1)*g.Length:], p.modulus)

	return &rnsDecomposer{
		params:    g,
		embedders: embedders,
	}
}

// Params returns the parameters of the decomposer.
func (d *rnsDecomposer) Params() GadgetParameters {
	return d.params
}

// Decompose decomposes a polynomial into a tensor using the RNS gadget.
func (d *rnsDecomposer) Decompose(p *crt.Poly) *Tensor {
	return nil
}

func (d *rnsDecomposer) DecomposeTo(pOut *Tensor, p *crt.Poly) {
	if p.IsNTT {
		panic("rnsDecomposer: polynomial must be in the coefficient form")
	}
}

// digitDecomposer is a struct that decomposes a polynomial into a tensor using the digit gadget.
type digitDecomposer struct {
	params   GadgetParameters
	embedder *crt.Embedder
}

// newDigitDecomposer creates a new [digitDecomposer] for the given parameters.
func newDigitDecomposer(g GadgetParameters, p Parameters) Decomposer {
	bitLen := 0.0
	for _, m := range p.modulus {
		bitLen += math.Log2(float64(m.Value()))
	}

	baseLen := int(math.Ceil(bitLen / float64(g.Length)))

	if baseLen >= num.MaxModulusBits {
		panic("newDigitDecomposer: base length exceeds max modulus bits")
	}

	embedder := crt.NewEmbedder([]*num.Modulus{num.NewModulus(1 << baseLen)}, p.modulus)

	return &digitDecomposer{
		params:   g,
		embedder: embedder,
	}
}

// Params returns the parameters of the decomposer.
func (d *digitDecomposer) Params() GadgetParameters {
	return d.params
}

// Decompose decomposes a polynomial into a tensor using the digit gadget.
func (d *digitDecomposer) Decompose(p *crt.Poly) *Tensor {
	return nil
}

func (d *digitDecomposer) DecomposeTo(pOut *Tensor, p *crt.Poly) {
	if p.IsNTT {
		panic("digitDecomposer: polynomial must be in the coefficient form")
	}
}
