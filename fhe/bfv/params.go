package bfv

import (
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

// EstimType is the type of the noise estimation.
type EstimType uint64

const (
	// TypeVariance is the variance of the noise.
	TypeVariance EstimType = iota
	// TypeWorstCase is the inf-norm bound of the noise.
	TypeWorstCase
)

// ParametersLiteral is a basic config for [Parameters].
// Must be compiled to [Parameters] before usage.
//
// # Warning
//
// Unless you are a cryptographic expert, DO NOT set these by yourself;
// always use the default parameters provided.
type ParametersLiteral struct {
	// RLWEParams is the parameters for the RLWE scheme.
	RLWEParams rlwe.ParametersLiteral

	// MessageModulus is the modulus for the messages and plaintexts.
	MessageModulus *num.Modulus

	// IsLevelled is whether the scheme is levelled.
	IsLevelled bool

	// EstimType is the type of the noise estimation.
	EstimType EstimType
}

// MustCompile compiles [ParametersLiteral] to [Parameters].
// Panics when non-sensible values exist.
func (p ParametersLiteral) Compile() Parameters {
	return Parameters{
		rlweParams:     p.RLWEParams.Compile(),
		messageModulus: p.MessageModulus,
		isLevelled:     p.IsLevelled,
		estimType:      p.EstimType,
	}
}

// Parameters is read-only parameters for the BFV scheme.
type Parameters struct {
	// rlweParams is the parameters for the RLWE scheme.
	rlweParams rlwe.Parameters

	// messageModulus is the modulus for the messages and plaintexts.
	messageModulus *num.Modulus

	// isLevelled is whether the scheme is levelled.
	isLevelled bool

	// estimType is the type of the noise estimation.
	estimType EstimType
}

// RLWEParams returns the parameters for the RLWE scheme.
func (p Parameters) RLWEParams() rlwe.Parameters {
	return p.rlweParams
}

// MessageModulus returns the modulus for the messages and plaintexts.
func (p Parameters) MessageModulus() *num.Modulus {
	return p.messageModulus
}

// IsLevelled returns whether the scheme is levelled.
func (p Parameters) IsLevelled() bool {
	return p.isLevelled
}

// EstimType returns the type of the noise estimation.
func (p Parameters) EstimType() EstimType {
	return p.estimType
}
