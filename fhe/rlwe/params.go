package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

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
	// Modulus is the modulus for encryption.
	Modulus []*num.Modulus
	// AuxModulus is the auxiliary or "special" modulus.
	AuxModulus []*num.Modulus

	// GadgetParams is the parameters for the gadget.
	GadgetParams GadgetParameters
	// SecretKeyParams is the parameters for secret key sampler.
	SecretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	NoiseParams crt.SamplerParameters
}

// Compile compiles [ParametersLiteral] to [Parameters].
// Panics when non-sensible values exist.
func (p ParametersLiteral) Compile() Parameters {
	// If GadgetParams is not set, use the default RNS gadget parameters.
	var gadgetParams GadgetParameters
	if p.GadgetParams != nil {
		gadgetParams = p.GadgetParams
	} else {
		gadgetParams = NewRNSGadgetParameters(1)
	}

	rP := p.RingParams
	_ = p.SecretKeyParams.Sampler()
	_ = p.NoiseParams.Sampler()

	return Parameters{
		ringParams: rP,
		modulus:    p.Modulus,
		auxModulus: p.AuxModulus,

		gadgetParams:    gadgetParams,
		secretKeyParams: p.SecretKeyParams,
		noiseParams:     p.NoiseParams,
	}
}

// Parameters is read-only parameters for the BGV scheme.
type Parameters struct {
	// RingParams is the parameters for underlying ring.
	ringParams dft.RingParameters
	// Modulus is the modulus for encryption.
	modulus []*num.Modulus
	// AuxModulus is the auxiliary or "special" modulus.
	auxModulus []*num.Modulus

	// GadgetParams is the parameters for the gadget.
	gadgetParams GadgetParameters
	// SecretKeyParams is the parameters for secret key sampler.
	secretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	noiseParams crt.SamplerParameters
}

func (p Parameters) RingParams() dft.RingParameters {
	return p.ringParams
}

func (p Parameters) Modulus() []*num.Modulus {
	return p.modulus
}

func (p Parameters) AuxModulus() []*num.Modulus {
	return p.auxModulus
}

func (p Parameters) GadgetParams() GadgetParameters {
	return p.gadgetParams
}

func (p Parameters) SecretKeyParams() crt.SamplerParameters {
	return p.secretKeyParams
}

func (p Parameters) NoiseParams() crt.SamplerParameters {
	return p.noiseParams
}
