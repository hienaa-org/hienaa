package rlwe

import "github.com/hienaa-org/hienaa/math/crt"

type Encryptor interface {
	SampleRlwe() *Ciphertext
	Encrypt(p *crt.Poly) *Ciphertext
}

type SKEncryptor struct {
}

type PKEncryptor struct {
}
