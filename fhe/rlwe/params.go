package rlwe

import (
	"math"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
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

// GadgetParametersLiteral is a basic config for [GadgetParameters].
// Must be compiled to [GadgetParameters] before usage.
//
// Currently, two types of gadgets are supported:
// [RNSGadgetParametersLiteral] and [DigitGadgetParametersLiteral].
type GadgetParametersLiteral interface {
	// Compile compiles [GadgetParametersLiteral] to [GadgetParameters].
	// Panics when non-sensible values exist.
	Compile() GadgetParameters
	isGagdetParametersLiteral()
}

// RNSGadgetParametersLiteral is the [GadgetParametersLiteral] for the RNS decomposition.
type RNSGadgetParametersLiteral struct {
	// ChunkSize is the size of the chunk.
	ChunkSize int
}

func (p RNSGadgetParametersLiteral) Compile() GadgetParameters {
	if p.ChunkSize <= 0 {
		panic("chunkSize must be positive")
	}

	return RNSGadgetParameters{
		chunkSize: p.ChunkSize,
	}
}

func (RNSGadgetParametersLiteral) isGagdetParametersLiteral() {}

// DigitGadgetParametersLiteral is the [GadgetParametersLiteral] for the digit decomposition.
type DigitGadgetParametersLiteral struct {
	// LogDigitBase is the logarithm of the base of the digit decomposition.
	LogDigitBase int
}

func (p DigitGadgetParametersLiteral) Compile() GadgetParameters {
	if p.LogDigitBase <= 0 || p.LogDigitBase > num.MaxModulusBits {
		panic("logDigitBase must be between 1 and max modulus bits")
	}

	return DigitGadgetParameters{
		logDigitBase: p.LogDigitBase,
	}
}

func (DigitGadgetParametersLiteral) isGagdetParametersLiteral() {}

// GadgetParameters is read-only parameters for gadget decomposition.
//
// Currenty, two types of gadgets are supported:
// [RNSGadgetParameters] and [DigitGadgetParameters].
type GadgetParameters interface {
	isGadgetParameters()
}

// RNSGadgetParameters is the [GadgetParameters] for the RNS decomposition.
type RNSGadgetParameters struct {
	// ChunkSize is the size of the chunk.
	chunkSize int
}

// ChunkSize is the size of the chunk.
func (p RNSGadgetParameters) ChunkSize() int {
	return p.chunkSize
}

func (RNSGadgetParameters) isGadgetParameters() {}

// DigitGadgetParameters is the [GadgetParameters] for the digit decomposition.
type DigitGadgetParameters struct {
	// LogDigitBase is the logarithm of the base of the digit decomposition.
	logDigitBase int
}

// LogDigitBase is the logarithm of the base of the digit decomposition.
func (p DigitGadgetParameters) LogDigitBase() int {
	return p.logDigitBase
}

func (DigitGadgetParameters) isGadgetParameters() {}

// ParametersLiteral is a basic config for [Parameters].
// Must be compiled to [Parameters] before usage.
//
// # Warning
//
// Unless you are a cryptographic expert, DO NOT set these by yourself;
// always use the default parameters provided.
type ParametersLiteral struct {
	// RingParams is the parameters for underlying ring.
	RingParams dft.RingParameters
	// BaseModulus is the modulus for encryption.
	BaseModulus []*num.Modulus
	// AuxModulus is the auxiliary or "special" modulus.
	AuxModulus []*num.Modulus

	// GadgetParams is the parameters for the gadget.
	// If nil, RNS gadget with chunk size 1 is used.
	GadgetParams GadgetParametersLiteral
	// SecretKeyParams is the parameters for secret key sampler.
	SecretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	NoiseParams crt.SamplerParameters
}

// Compile compiles [ParametersLiteral] to [Parameters].
// Panics when non-sensible values exist.
func (p ParametersLiteral) Compile() Parameters {
	gadgetParams := p.GadgetParams
	if gadgetParams == nil {
		gadgetParams = RNSGadgetParametersLiteral{ChunkSize: 1}
	}

	fullMod := make([]*num.Modulus, 0, len(p.AuxModulus)+len(p.BaseModulus))
	fullMod = append(fullMod, p.AuxModulus...)
	fullMod = append(fullMod, p.BaseModulus...)

	_ = p.SecretKeyParams.Sampler()
	_ = p.NoiseParams.Sampler()

	return Parameters{
		crtOp:   crt.NewOperator(p.RingParams, fullMod),
		fullMod: fullMod,
		baseMod: fullMod[len(p.AuxModulus):],
		auxMod:  fullMod[:len(p.AuxModulus)],

		gadgetParams:    gadgetParams.Compile(),
		secretKeyParams: p.SecretKeyParams,
		noiseParams:     p.NoiseParams,
	}
}

// Parameters is read-only parameters for the BGV scheme.
type Parameters struct {
	// crtOp is an underlying [crt.Operator].
	crtOp *crt.Operator
	// fullMod is the full modulus chain.
	fullMod []*num.Modulus
	// BaseModulus is the Modulus for encryption.
	baseMod []*num.Modulus
	// AuxModulus is the auxiliary or "special" modulus.
	auxMod []*num.Modulus

	// GadgetParams is the parameters for the gadget.
	gadgetParams GadgetParameters
	// SecretKeyParams is the parameters for secret key sampler.
	secretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	noiseParams crt.SamplerParameters
}

// Operator is the underlying [crt.Operator].
func (p Parameters) Operator() *crt.Operator {
	return p.crtOp
}

// RingParams is the parameters for underlying ring.
func (p Parameters) RingParams() dft.RingParameters {
	return p.crtOp.Params()
}

// Rank returns the rank of the ring.
func (p Parameters) Rank() int {
	return p.crtOp.Params().Rank()
}

// BaseModulus is the modulus for encryption.
func (p Parameters) BaseModulus() []*num.Modulus {
	return p.baseMod
}

// HasAuxModulus returns whether an auxiliary modulus exists.
func (p Parameters) HasAuxModulus() bool {
	return len(p.auxMod) > 0
}

// AuxModulus is the auxiliary or "special" modulus.
func (p Parameters) AuxModulus() []*num.Modulus {
	return p.auxMod
}

// FullModulus is the full modulus chain.
// Equals to AuxModulus || BaseModulus.
func (p Parameters) FullModulus() []*num.Modulus {
	return p.fullMod
}

// GadgetParams is the parameters for the gadget.
func (p Parameters) GadgetParams() GadgetParameters {
	return p.gadgetParams
}

// GadgetLen returns the length of the gadget decomposition.
func (p Parameters) GadgetLen() int {
	switch gadgetParams := p.gadgetParams.(type) {
	case RNSGadgetParameters:
		return int(math.Ceil(float64(len(p.baseMod)) / float64(gadgetParams.chunkSize)))
	case DigitGadgetParameters:
		modBitLen := 0.0
		for _, q := range p.baseMod {
			modBitLen += num.Log2(q.Value())
		}
		return int(math.Ceil(modBitLen / float64(gadgetParams.logDigitBase)))
	}
	return 0
}

// SecretKeyParams is the parameters for secret key sampler.
func (p Parameters) SecretKeyParams() crt.SamplerParameters {
	return p.secretKeyParams
}

// NoiseParams is the parameters for noise sampler.
func (p Parameters) NoiseParams() crt.SamplerParameters {
	return p.noiseParams
}
