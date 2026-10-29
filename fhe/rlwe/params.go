package rlwe

import (
	"math/big"
	"math/bits"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
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
	isGadgetParamsLiteral()
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

func (RNSGadgetParametersLiteral) isGadgetParamsLiteral() {}

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

func (DigitGadgetParametersLiteral) isGadgetParamsLiteral() {}

// GadgetParameters is read-only parameters for gadget decomposition.
//
// Currently, two types of gadgets are supported:
// [RNSGadgetParameters] and [DigitGadgetParameters].
type GadgetParameters interface {
	isGadgetParams()
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

func (RNSGadgetParameters) isGadgetParams() {}

// DigitGadgetParameters is the [GadgetParameters] for the digit decomposition.
type DigitGadgetParameters struct {
	// LogDigitBase is the logarithm of the base of the digit decomposition.
	logDigitBase int
}

// LogDigitBase is the logarithm of the base of the digit decomposition.
func (p DigitGadgetParameters) LogDigitBase() int {
	return p.logDigitBase
}

func (DigitGadgetParameters) isGadgetParams() {}

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
	if len(p.BaseModulus) == 0 {
		panic("BaseModulus should be non-empty")
	}

	gadgetParamsLiteral := p.GadgetParams
	if gadgetParamsLiteral == nil {
		gadgetParamsLiteral = RNSGadgetParametersLiteral{ChunkSize: 1}
	}
	gadgetParams := gadgetParamsLiteral.Compile()

	fullMod := make([]*num.Modulus, 0, len(p.AuxModulus)+len(p.BaseModulus))
	fullMod = append(fullMod, p.AuxModulus...)
	fullMod = append(fullMod, p.BaseModulus...)

	tmpBig := new(big.Int)
	baseModBig := big.NewInt(1)
	for i := range p.BaseModulus {
		tmpBig.SetUint64(p.BaseModulus[i].Value())
		baseModBig.Mul(baseModBig, tmpBig)
	}

	auxModBig := big.NewInt(1)
	for i := range p.AuxModulus {
		tmpBig.SetUint64(p.AuxModulus[i].Value())
		auxModBig.Mul(auxModBig, tmpBig)
	}

	fullModBig := new(big.Int).Mul(baseModBig, auxModBig)

	_ = p.SecretKeyParams.Sampler()
	_ = p.NoiseParams.Sampler()

	var gadgetVec []*Element
	switch gadgetParams := gadgetParams.(type) {
	case RNSGadgetParameters:
		gadgetVec = make([]*Element, num.DivCeil(len(p.BaseModulus), gadgetParams.chunkSize))
		g, chunk := new(big.Int), new(big.Int)
		for i := range gadgetVec {
			start := i * gadgetParams.chunkSize
			end := min((i+1)*gadgetParams.chunkSize, len(p.BaseModulus))

			g.SetUint64(1)
			chunk.SetUint64(1)

			for j := range p.BaseModulus {
				tmpBig.SetUint64(p.BaseModulus[j].Value())
				if start <= j && j < end {
					chunk.Mul(chunk, tmpBig)
				} else {
					g.Mul(g, tmpBig)
				}
			}

			tmpBig.ModInverse(g, chunk)
			g.Mul(g, tmpBig)
			g.Mod(g, baseModBig)
			g.Mul(g, auxModBig)

			gadgetVec[i] = &Element{
				Value:  crt.NewScalarFrom(g, fullMod),
				auxLen: len(p.AuxModulus),
			}
		}

	case DigitGadgetParameters:
		for i := range p.BaseModulus {
			if p.BaseModulus[i].Value()%2 == 0 {
				panic("base modulus should be odd in digit decomposition")
			}
		}

		modBitLen := 0
		for _, q := range p.BaseModulus {
			modBitLen += bits.Len64(q.Value() - 1)
		}
		gadgetVec = make([]*Element, num.DivCeil(modBitLen, gadgetParams.logDigitBase))

		g := new(big.Int).Set(auxModBig)
		if len(p.AuxModulus) == 0 {
			g.SetUint64(1)
		}

		for i := range gadgetVec {
			gadgetVec[i] = &Element{
				Value:  crt.NewScalarFrom(g, fullMod),
				auxLen: len(p.AuxModulus),
			}
			g.Lsh(g, uint(gadgetParams.logDigitBase))
		}
	}

	if len(p.AuxModulus) == 0 {
		auxModBig.SetUint64(0)
	}

	return Parameters{
		crtOp:   crt.NewOperator(p.RingParams, fullMod),
		fullMod: fullMod,
		baseMod: fullMod[len(p.AuxModulus):],
		auxMod:  fullMod[:len(p.AuxModulus)],

		fullModBig: fullModBig,
		baseModBig: baseModBig,
		auxModBig:  auxModBig,

		gadgetParams:    gadgetParams,
		gadgetVec:       gadgetVec,
		secretKeyParams: p.SecretKeyParams,
		noiseParams:     p.NoiseParams,
	}
}

// Parameters is read-only parameters for RLWE.
type Parameters struct {
	// crtOp is an underlying [crt.Operator].
	crtOp *crt.Operator
	// fullMod is the full modulus chain.
	fullMod []*num.Modulus
	// BaseModulus is the Modulus for encryption.
	baseMod []*num.Modulus
	// AuxModulus is the auxiliary or "special" modulus.
	auxMod []*num.Modulus

	// fullModBig is a product of fullMod.
	fullModBig *big.Int
	// baseModBig is a product of baseMod.
	baseModBig *big.Int
	// auxModBig is a product of auxMod.
	auxModBig *big.Int

	// GadgetParams is the parameters for the gadget.
	gadgetParams GadgetParameters
	// gadgetVec is the gadget vector.
	gadgetVec []*Element
	// SecretKeyParams is the parameters for secret key sampler.
	secretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	noiseParams crt.SamplerParameters
}

// Op is the underlying [crt.Operator].
func (p Parameters) Op() *crt.Operator {
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

// BaseMod is the modulus for encryption.
func (p Parameters) BaseMod() []*num.Modulus {
	return p.baseMod
}

// AuxMod is the auxiliary or "special" modulus.
func (p Parameters) AuxMod() []*num.Modulus {
	return p.auxMod
}

// FullMod is the full modulus chain.
// Equals to AuxMod || BaseMod.
func (p Parameters) FullMod() []*num.Modulus {
	return p.fullMod
}

// BaseModLen is the length of base modulus.
func (p Parameters) BaseModLen() int {
	return len(p.baseMod)
}

// AuxModLen is the length of auxiliary modulus.
func (p Parameters) AuxModLen() int {
	return len(p.auxMod)
}

// FullModLen is the length of full modulus chain.
func (p Parameters) FullModLen() int {
	return len(p.fullMod)
}

// BaseModBig is the product of BaseMod.
func (p Parameters) BaseModBig() *big.Int {
	return p.baseModBig
}

// AuxModBig is the product of AuxMod.
func (p Parameters) AuxModBig() *big.Int {
	return p.auxModBig
}

// FullModBig is the product of FullMod.
func (p Parameters) FullModBig() *big.Int {
	return p.fullModBig
}

// GadgetParams is the parameters for the gadget.
func (p Parameters) GadgetParams() GadgetParameters {
	return p.gadgetParams
}

// GadgetLen returns the length of the gadget decomposition.
func (p Parameters) GadgetLen() int {
	return len(p.gadgetVec)
}

// GadgetVec returns the gadget vector.
func (p Parameters) GadgetVec() []*Element {
	return p.gadgetVec
}

// SecretKeyParams is the parameters for secret key sampler.
func (p Parameters) SecretKeyParams() crt.SamplerParameters {
	return p.secretKeyParams
}

// NoiseParams is the parameters for noise sampler.
func (p Parameters) NoiseParams() crt.SamplerParameters {
	return p.noiseParams
}
