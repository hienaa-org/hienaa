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
	// LogBaseMod is the bit length of base modulus.
	LogBaseMod float64
	// LogAuxMod is the bit length of auxililary modulus.
	LogAuxMod float64

	// baseMod is the base modulus. This overrides LogBaseModulus.
	baseMod []*num.Modulus
	// auxMod is the auxiliary modulus. This overrides LogAuxModulus.
	auxMod []*num.Modulus

	// GadgetParams is the parameters for the gadget.
	// If nil, RNS gadget with chunk size same as max(1, auxiliary modulus length) is used.
	GadgetParams GadgetParametersLiteral
	// SecretKeyParams is the parameters for secret key sampler.
	SecretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	NoiseParams crt.SamplerParameters
}

// WithBaseMod overrides LogBaseMod and sets the base modulus directly.
// Must be called before compile.
func (p ParametersLiteral) WithBaseMod[Q uint64 | *num.Modulus](baseMod []Q) ParametersLiteral {
	if baseMod == nil {
		p.baseMod = nil
		return p
	}

	switch baseMod := any(baseMod).(type) {
	case []uint64:
		p.baseMod = make([]*num.Modulus, len(baseMod))
		for i := range baseMod {
			p.baseMod[i] = num.NewModulus(baseMod[i])
		}
	case []*num.Modulus:
		p.baseMod = baseMod
	}
	return p
}

// WithAuxMod overrides LogAuxMod and sets the auxiliary modulus directly.
// Must be called before compile.
func (p ParametersLiteral) WithAuxMod[Q uint64 | *num.Modulus](auxMod []Q) ParametersLiteral {
	if auxMod == nil {
		p.auxMod = nil
		return p
	}

	switch auxMod := any(auxMod).(type) {
	case []uint64:
		p.auxMod = make([]*num.Modulus, len(auxMod))
		for i := range auxMod {
			p.auxMod[i] = num.NewModulus(auxMod[i])
		}
	case []*num.Modulus:
		p.auxMod = auxMod
	}
	return p
}

// Compile compiles [ParametersLiteral] to [Parameters].
// Panics when non-sensible values exist.
func (p ParametersLiteral) Compile() Parameters {
	var baseMod, auxMod []*num.Modulus
	switch {
	case p.baseMod == nil && p.auxMod == nil:
		baseMod, auxMod = FindNTTPrimes(p.RingParams, float64(p.LogBaseMod), float64(p.LogAuxMod))

	case p.baseMod == nil && p.auxMod != nil:
		auxMod = p.auxMod
		auxModMap := make(map[uint64]struct{}, len(auxMod))
		for i := range auxMod {
			auxModMap[auxMod[i].Value()] = struct{}{}
		}
		baseMod = prevPrimes(p.RingParams, p.LogBaseMod, auxModMap)

	case p.baseMod != nil && p.auxMod == nil:
		baseMod = p.baseMod
		baseModMap := make(map[uint64]struct{}, len(baseMod))
		for i := range baseMod {
			baseModMap[baseMod[i].Value()] = struct{}{}
		}
		auxMod = prevPrimes(p.RingParams, p.LogAuxMod, baseModMap)

	case p.baseMod != nil && p.auxMod != nil:
		baseMod, auxMod = p.baseMod, p.auxMod
	}

	if len(baseMod) == 0 {
		panic("base modulus must be nonempty")
	}

	fullMod := make([]*num.Modulus, 0, len(auxMod)+len(baseMod))
	fullMod = append(fullMod, auxMod...)
	fullMod = append(fullMod, baseMod...)

	tmpBig := new(big.Int)
	baseModBig := big.NewInt(1)
	for i := range baseMod {
		tmpBig.SetUint64(baseMod[i].Value())
		baseModBig.Mul(baseModBig, tmpBig)
	}

	auxModBig := big.NewInt(1)
	for i := range auxMod {
		tmpBig.SetUint64(auxMod[i].Value())
		auxModBig.Mul(auxModBig, tmpBig)
	}

	fullModBig := new(big.Int).Mul(baseModBig, auxModBig)

	_ = p.SecretKeyParams.Sampler()
	_ = p.NoiseParams.Sampler()

	gadgetParamsLiteral := p.GadgetParams
	if gadgetParamsLiteral == nil {
		gadgetParamsLiteral = RNSGadgetParametersLiteral{ChunkSize: max(1, len(auxMod))}
	}
	gadgetParams := gadgetParamsLiteral.Compile()

	var gadgetVec []*Element
	switch gadgetParams := gadgetParams.(type) {
	case RNSGadgetParameters:
		gadgetVec = make([]*Element, num.DivCeil(len(baseMod), gadgetParams.chunkSize))
		g, chunk := new(big.Int), new(big.Int)
		for i := range gadgetVec {
			start := i * gadgetParams.chunkSize
			end := min((i+1)*gadgetParams.chunkSize, len(baseMod))

			g.SetUint64(1)
			chunk.SetUint64(1)

			for j := range baseMod {
				tmpBig.SetUint64(baseMod[j].Value())
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
				auxLen: len(auxMod),
			}
		}

	case DigitGadgetParameters:
		for i := range baseMod {
			if baseMod[i].Value()&1 == 0 {
				panic("base modulus should be odd in digit decomposition")
			}
		}

		modBitLen := 0
		for _, q := range baseMod {
			modBitLen += bits.Len64(q.Value() - 1)
		}
		gadgetVec = make([]*Element, num.DivCeil(modBitLen, gadgetParams.logDigitBase))

		g := new(big.Int).Set(auxModBig)
		if len(auxMod) == 0 {
			g.SetUint64(1)
		}

		for i := range gadgetVec {
			gadgetVec[i] = &Element{
				Value:  crt.NewScalarFrom(g, fullMod),
				auxLen: len(auxMod),
			}
			g.Lsh(g, uint(gadgetParams.logDigitBase))
		}
	}

	if len(auxMod) == 0 {
		auxModBig.SetUint64(0)
	}

	return Parameters{
		crtOp:   crt.NewOperator(p.RingParams, fullMod),
		fullMod: fullMod,
		baseMod: fullMod[len(auxMod):],
		auxMod:  fullMod[:len(auxMod)],

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
