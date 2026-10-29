package rlwe

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// Encryptor encrypts and decrypts RLWE ciphertexts.
type Encryptor struct {
	params Parameters

	sk *SecretKey
	pk *PublicKey

	plainOp *PlainOperator
	dcmp    *Decomposer

	maskSampler   *crt.UniformSampler
	auxKeySampler crt.Sampler
	noiseSampler  crt.Sampler

	poolScalar *pool.Pool[*Element]
	poolPoly   *pool.Pool[*Element]
	poolRLWE   *pool.Pool[*Ciphertext]
}

// NewEncryptor creates a new [Encryptor].
// Secret key is automatically sampled.
func NewEncryptor(params Parameters) *Encryptor {
	sk := &SecretKey{
		Value:  params.secretKeyParams.Sampler().Sample(params.Rank(), params.fullMod),
		auxLen: params.AuxModLen(),
	}
	params.crtOp.FwdNTTTo(sk.Value, sk.Value)

	return newEncryptorWithKey(params, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] with supplied key.
// The supplied key is copied.
// The keys are assumed to be in full modulus.
func NewEncryptorWithKey[K *SecretKey | *PublicKey](params Parameters, key K) *Encryptor {
	switch key := any(key).(type) {
	case *SecretKey:
		checkElementShape(params, (*Element)(key))
		if key.FullModLen() != params.FullModLen() || key.Value.Type != crt.TypePoly {
			panic("inconsistent inputs")
		}

		sk := key.Copy()
		if sk.Form() == dft.FormCoeff {
			params.crtOp.FwdNTTTo(sk.Value, sk.Value)
		}
		return newEncryptorWithKey(params, sk)

	case *PublicKey:
		checkCiphertextShape(params, (*Ciphertext)(key))
		if key.FullModLen() != params.FullModLen() {
			panic("inconsistent inputs")
		}

		pk := key.Copy()
		if pk.Form() == dft.FormCoeff {
			params.crtOp.FwdNTTTo(pk.Body.Value, pk.Body.Value)
			params.crtOp.FwdNTTTo(pk.Mask.Value, pk.Mask.Value)
		}
		return newEncryptorWithKey(params, pk)
	}

	return nil
}

func newEncryptorWithKey[K *SecretKey | *PublicKey](params Parameters, key K) *Encryptor {
	enc := &Encryptor{
		params: params,

		plainOp: NewPlainOperator(params),
		dcmp:    NewDecomposer(params),

		maskSampler:  crt.UniformSamplerParameters{}.Sampler().(*crt.UniformSampler),
		noiseSampler: params.noiseParams.Sampler(),

		poolScalar: pool.NewPool(func() *Element {
			return NewScalar(params, params.AuxModLen() > 0)
		}),
		poolPoly: pool.NewPool(func() *Element {
			return NewPoly(params, params.AuxModLen() > 0, dft.FormCoeff)
		}),
		poolRLWE: pool.NewPool(func() *Ciphertext {
			return NewCiphertext(params, params.AuxModLen() > 0, dft.FormCoeff)
		}),
	}

	switch key := any(key).(type) {
	case *SecretKey:
		enc.sk = key
	case *PublicKey:
		enc.pk = key
		enc.auxKeySampler = params.secretKeyParams.Sampler()
	}

	return enc
}

// Params returns the [Parameters] of the encryptor.
func (enc *Encryptor) Params() Parameters {
	return enc.params
}

// PublicKey returns the [PublicKey] of the encryptor.
// It returns nil if enc is not a public-key encryptor.
func (enc *Encryptor) PublicKey() *PublicKey {
	return enc.pk
}

// PlainOp returns the [PlainOperator] of the encryptor.
func (enc *Encryptor) PlainOp() *PlainOperator {
	return enc.plainOp
}

// Decomposer returns the [Decomposer] of the encryptor.
func (enc *Encryptor) Decomposer() *Decomposer {
	return enc.dcmp
}

// HasPublicKey returns if enc is a public key encryptor.
func (enc *Encryptor) HasPublicKey() bool {
	return enc.pk != nil
}
