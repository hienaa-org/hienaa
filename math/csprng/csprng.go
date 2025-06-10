// Package csprng implements various samplers used throughout HE schemes.
// All samplers in this package are CSPRNG (Cryptographically Secure Pseudo-Random Generator),
// hence the name.
package csprng

// VecSamplerBuilder is an interface for parameters of [VecSampler].
type VecSamplerBuilder interface {
	// Build returns a new [VecSampler] instance.
	Build() VecSampler
}

// VecSampler is an interface for sampling uint64 vectors.
type VecSampler interface {
	// SampleVecTo samples and fills the provided uint64 vector.
	SampleVecTo(vec []uint64)
	// Copy returns a thread-safe copy of the sampler.
	Copy() VecSampler
}
