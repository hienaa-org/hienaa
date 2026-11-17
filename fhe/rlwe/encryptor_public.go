package rlwe

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
)

// PublicEncryptor encrypts RLWE ciphertexts with a public key.
type PublicEncryptor struct {
	params Parameters

	pk *Ciphertext

	plainOp *PlainOperator
	dcmp    *Decomposer

	auxKeySampler crt.Sampler
	noiseSampler  crt.Sampler

	poolScalar *pool.Pool[*Element]
	poolPoly   *pool.Pool[*Element]
	poolRLWE   *pool.Pool[*Ciphertext]
}

// NewPublicEncryptor creates a new [PublicEncryptor] with the supplied public key.
// The supplied key is copied.
// The key is assumed to be in full modulus.
func NewPublicEncryptor(params Parameters, key *PublicKey) *PublicEncryptor {
	checkCiphertextShape(params, (*Ciphertext)(key))
	if key.FullModLen() != params.FullModLen() {
		panic("inconsistent inputs")
	}

	pk := key.Copy()
	if pk.Form() == dft.FormCoeff {
		params.crtOp.FwdNTTTo(pk.Body.Value, pk.Body.Value)
		params.crtOp.FwdNTTTo(pk.Mask.Value, pk.Mask.Value)
	}

	return &PublicEncryptor{
		params: params,

		pk: (*Ciphertext)(pk),

		plainOp: NewPlainOperator(params),
		dcmp:    NewDecomposer(params),

		auxKeySampler: params.secretKeyParams.Sampler(),
		noiseSampler:  params.noiseParams.Sampler(),

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
}

// Params returns the [Parameters] of the encryptor.
func (enc *PublicEncryptor) Params() Parameters {
	return enc.params
}

// PlainOp returns the [PlainOperator] of the encryptor.
func (enc *PublicEncryptor) PlainOp() *PlainOperator {
	return enc.plainOp
}

// Decomposer returns the [Decomposer] of the encryptor.
func (enc *PublicEncryptor) Decomposer() *Decomposer {
	return enc.dcmp
}

// PublicKey returns the [PublicKey] of the encryptor.
func (enc *PublicEncryptor) PublicKey() *PublicKey {
	return (*PublicKey)(enc.pk)
}

// Encrypt returns an encryption of pt.
// If pt is nil, it outputs a zero encryption.
func (enc *PublicEncryptor) Encrypt[C *Ciphertext | *GadgetEncryption | *RGSWEncryption](pt *Element, outForm dft.Form) C {
	baseLen, auxLen := enc.params.BaseModLen(), enc.params.AuxModLen()
	if pt != nil {
		baseLen, auxLen = pt.BaseModLen(), pt.AuxModLen()
	}

	var ctOut C
	switch any(ctOut).(type) {
	case *Ciphertext:
		ctOut = any(NewCiphertextCustom(enc.params.Rank(), baseLen, auxLen, dft.FormCoeff)).(C)
	case *GadgetEncryption:
		ctOut = any(NewGadgetEncryptionCustom(enc.params.Rank(), baseLen, auxLen, enc.params.GadgetLen(), dft.FormCoeff)).(C)
	case *RGSWEncryption:
		ctOut = any(NewRGSWEncryptionCustom(enc.params.Rank(), baseLen, auxLen, enc.params.GadgetLen(), dft.FormCoeff)).(C)
	}
	enc.EncryptTo(ctOut, pt, outForm)
	return ctOut
}

// EncryptTo encrypts pt to ctOut.
// If pt is nil, ctOut is a zero encryption.
func (enc *PublicEncryptor) EncryptTo[C *Ciphertext | *GadgetEncryption | *RGSWEncryption](ctOut C, pt *Element, outForm dft.Form) {
	ctOutRLWE, ctOutIsRLWE := any(ctOut).(*Ciphertext)
	ctOutGadEnc, ctOutIsGadEnc := any(ctOut).(*GadgetEncryption)
	ctOutRGSW, ctOutIsRGSW := any(ctOut).(*RGSWEncryption)

	switch {
	case ctOutIsRLWE:
		enc.encryptPublicKeyTo(ctOutRLWE, pt, outForm)

	case ctOutIsGadEnc:
		checkCiphertextShape(enc.params, ctOutGadEnc)
		if pt == nil {
			for i := 0; i < enc.params.GadgetLen(); i++ {
				enc.EncryptTo(ctOutGadEnc.Value[i], pt, outForm)
			}
			return
		}

		var ptBufFull *Element
		var ptCopyFull *Element
		switch {
		case pt.Type() == crt.TypeScalar:
			ptBufFull = enc.poolScalar.Get()
			ptCopyFull = enc.poolScalar.Get()
			defer enc.poolScalar.Put(ptBufFull)
			defer enc.poolScalar.Put(ptCopyFull)
		case pt.Type() == crt.TypePoly:
			ptBufFull = enc.poolPoly.Get()
			ptCopyFull = enc.poolPoly.Get()
			defer enc.poolPoly.Put(ptBufFull)
			defer enc.poolPoly.Put(ptCopyFull)
		}
		ptBuf := ptBufFull.WithModLen(ctOutGadEnc.BaseModLen(), ctOutGadEnc.AuxModLen())
		ptCopy := ptCopyFull.WithModLen(ctOutGadEnc.BaseModLen(), ctOutGadEnc.AuxModLen())
		ptCopy.CopyFrom(pt)

		for i := 0; i < enc.params.GadgetLen(); i++ {
			enc.plainOp.MulTo(ptBuf, enc.params.gadgetVec[i].WithModLen(ctOutGadEnc.BaseModLen(), ctOutGadEnc.AuxModLen()), ptCopy)
			enc.EncryptTo(ctOutGadEnc.Value[i], ptBuf, outForm)
		}

	case ctOutIsRGSW:
		checkCiphertextShape(enc.params, ctOutRGSW)
		if pt == nil {
			for i := 0; i < enc.params.GadgetLen(); i++ {
				enc.EncryptTo(ctOutRGSW.BodyStack.Value[i], pt, outForm)
				enc.EncryptTo(ctOutRGSW.MaskStack.Value[i], pt, outForm)
			}
			return
		}

		var ptBufFull *Element
		var ptCopyFull *Element
		switch {
		case pt.Type() == crt.TypeScalar:
			ptBufFull = enc.poolScalar.Get()
			ptCopyFull = enc.poolScalar.Get()
			defer enc.poolScalar.Put(ptBufFull)
			defer enc.poolScalar.Put(ptCopyFull)
		case pt.Type() == crt.TypePoly:
			ptBufFull = enc.poolPoly.Get()
			ptCopyFull = enc.poolPoly.Get()
			defer enc.poolPoly.Put(ptBufFull)
			defer enc.poolPoly.Put(ptCopyFull)
		}
		ptBuf := ptBufFull.WithModLen(ctOutRGSW.BaseModLen(), ctOutRGSW.AuxModLen())
		ptCopy := ptCopyFull.WithModLen(ctOutRGSW.BaseModLen(), ctOutRGSW.AuxModLen())

		ptCopy.CopyFrom(pt)
		if pt.Type() == crt.TypePoly && pt.Form() != outForm {
			switch {
			case outForm == dft.FormCoeff:
				enc.plainOp.InvNTTTo(ptCopy, ptCopy)
			case outForm == dft.FormNTT:
				enc.plainOp.FwdNTTTo(ptCopy, ptCopy)
			}
		}

		for i := 0; i < enc.params.GadgetLen(); i++ {
			enc.plainOp.MulTo(ptBuf, enc.params.gadgetVec[i].WithModLen(ctOutRGSW.BaseModLen(), ctOutRGSW.AuxModLen()), ptCopy)
			enc.EncryptTo(ctOutRGSW.BodyStack.Value[i], ptBuf, outForm)

			enc.EncryptTo(ctOutRGSW.MaskStack.Value[i], nil, outForm)
			enc.plainOp.AddTo(ctOutRGSW.MaskStack.Value[i].Mask, ctOutRGSW.MaskStack.Value[i].Mask, ptBuf)
		}
	}
}

func (enc *PublicEncryptor) encryptPublicKeyTo(ctOut *Ciphertext, pt *Element, outForm dft.Form) {
	checkCiphertextShape(enc.params, ctOut)
	if pt != nil {
		checkElementShape(enc.params, pt)
	}

	ptBuf := pt
	if pt != nil {
		var ptBufFull *Element
		switch {
		case pt.Type() == crt.TypeScalar:
			ptBufFull = enc.poolScalar.Get()
			defer enc.poolScalar.Put(ptBufFull)
		case pt.Type() == crt.TypePoly:
			ptBufFull = enc.poolPoly.Get()
			defer enc.poolPoly.Put(ptBufFull)
		}
		ptBuf = ptBufFull.WithModLen(pt.BaseModLen(), pt.AuxModLen())
		ptBuf.CopyFrom(pt)
	}

	ctMod := enc.params.fullMod[enc.params.AuxModLen()-ctOut.AuxModLen() : enc.params.AuxModLen()+ctOut.BaseModLen()]

	auxKeyFull := enc.poolPoly.Get()
	defer enc.poolPoly.Put(auxKeyFull)
	auxKey := auxKeyFull.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen())

	enc.auxKeySampler.SampleTo(auxKey.Value, ctMod)
	enc.plainOp.FwdNTTTo(auxKey, auxKey)

	ctOut.Clear()
	ctOut.SetForm(dft.FormNTT)
	enc.plainOp.MulSubTo(ctOut.Mask, auxKey, enc.pk.Mask.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen()))
	enc.plainOp.MulSubTo(ctOut.Body, auxKey, enc.pk.Body.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen()))

	noiseFull := enc.poolRLWE.Get()
	defer enc.poolRLWE.Put(noiseFull)
	noise := noiseFull.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen())

	enc.noiseSampler.SampleTo(noise.Body.Value, ctMod)
	enc.noiseSampler.SampleTo(noise.Mask.Value, ctMod)

	if ptBuf != nil {
		if ptBuf.Type() == crt.TypePoly && ptBuf.Form() == dft.FormNTT {
			enc.plainOp.AddTo(ctOut.Body, ctOut.Body, ptBuf)
		} else {
			enc.plainOp.AddTo(noise.Body, noise.Body, ptBuf)
		}
	}

	switch {
	case outForm == dft.FormCoeff:
		enc.plainOp.InvNTTTo(ctOut.Mask, ctOut.Mask)
		enc.plainOp.InvNTTTo(ctOut.Body, ctOut.Body)
	case outForm == dft.FormNTT:
		enc.plainOp.FwdNTTTo(noise.Body, noise.Body)
		enc.plainOp.FwdNTTTo(noise.Mask, noise.Mask)
	}

	enc.plainOp.AddTo(ctOut.Body, ctOut.Body, noise.Body)
	enc.plainOp.AddTo(ctOut.Mask, ctOut.Mask, noise.Mask)
}
