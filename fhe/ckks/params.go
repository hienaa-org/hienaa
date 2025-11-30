package ckks

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

	// SecretKeyParams is the parameters for secret key sampler.
	SecretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	NoiseParams crt.SamplerParameters
}

// MustCompile compiles [ParametersLiteral] to [Parameters].
// Panics when non-sensible values exist.
func (p ParametersLiteral) Compile() Parameters {
	_ = p.SecretKeyParams.Sampler()
	_ = p.NoiseParams.Sampler()

	return Parameters{
		ringParams: p.RingParams,
		modulus:    p.Modulus,
		auxModulus: p.AuxModulus,

		secretKeyParams: p.SecretKeyParams,
		noiseParams:     p.NoiseParams,
	}
}

// Parameters is read-only parameters for the CKKS scheme.
type Parameters struct {
	// RingParams is the parameters for underlying ring.
	ringParams dft.RingParameters
	// Modulus is the modulus for encryption.
	modulus []*num.Modulus
	// AuxModulus is the auxiliary or "special" modulus.
	auxModulus []*num.Modulus

	// SecretKeyParams is the parameters for secret key sampler.
	secretKeyParams crt.SamplerParameters
	// NoiseParams is the parameters for noise sampler.
	noiseParams crt.SamplerParameters
}
