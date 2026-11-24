package rlwe

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encryptor encrypts and decrypts RLWE ciphertexts with a secret key.
type Encryptor struct {
	params Parameters

	sk *Element

	plainOp *PlainOperator
	dcmp    *Decomposer

	maskSampler  *crt.UniformSampler
	noiseSampler crt.Sampler

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

// NewEncryptorWithKey creates a new [Encryptor] with the supplied secret key.
// The supplied key is copied.
// The key is assumed to be in full modulus.
func NewEncryptorWithKey(params Parameters, key *SecretKey) *Encryptor {
	checkElementShape(params, (*Element)(key))
	if key.FullModLen() != params.FullModLen() || key.Value.Type != crt.TypePoly {
		panic("inconsistent inputs")
	}

	sk := key.Copy()
	if sk.Form() == dft.FormCoeff {
		params.crtOp.FwdNTTTo(sk.Value, sk.Value)
	}
	return newEncryptorWithKey(params, sk)
}

func newEncryptorWithKey(params Parameters, key *SecretKey) *Encryptor {
	return &Encryptor{
		params: params,

		sk: (*Element)(key),

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
}

// Params returns the [Parameters] of the encryptor.
func (enc *Encryptor) Params() Parameters {
	return enc.params
}

// PlainOp returns the [PlainOperator] of the encryptor.
func (enc *Encryptor) PlainOp() *PlainOperator {
	return enc.plainOp
}

// Decomposer returns the [Decomposer] of the encryptor.
func (enc *Encryptor) Decomposer() *Decomposer {
	return enc.dcmp
}

// Key returns the [SecretKey] of the encryptor.
func (enc *Encryptor) Key() *SecretKey {
	return (*SecretKey)(enc.sk)
}

// Encrypt returns an encryption of pt.
// If pt is nil, it outputs a zero encryption.
func (enc *Encryptor) Encrypt[C *Ciphertext | *GadgetEncryption | *RGSWEncryption](pt *Element, outForm dft.Form) C {
	baseLen, auxLen := enc.params.BaseModLen(), enc.params.AuxModLen()
	if pt != nil {
		baseLen, auxLen = pt.BaseModLen(), pt.AuxModLen()
	}

	var ctOut C
	switch any(ctOut).(type) {
	case *Ciphertext:
		ctOut = any(NewCiphertextCustom(enc.params.Rank(), baseLen, auxLen, dft.FormCoeff)).(C)
	case *GadgetEncryption:
		ctOut = any(NewGadgetEncryptionCustom(enc.params.Rank(), baseLen, enc.params.AuxModLen(), enc.params.GadgetLen(), dft.FormCoeff)).(C)
	case *RGSWEncryption:
		ctOut = any(NewRGSWEncryptionCustom(enc.params.Rank(), baseLen, enc.params.AuxModLen(), enc.params.GadgetLen(), dft.FormCoeff)).(C)
	}
	enc.EncryptTo(ctOut, pt, outForm)
	return ctOut
}

// EncryptTo encrypts pt to ctOut.
// If pt is nil, ctOut is a zero encryption.
func (enc *Encryptor) EncryptTo[C *Ciphertext | *GadgetEncryption | *RGSWEncryption](ctOut C, pt *Element, outForm dft.Form) {
	ctOutRLWE, ctOutIsRLWE := any(ctOut).(*Ciphertext)
	ctOutGadEnc, ctOutIsGadEnc := any(ctOut).(*GadgetEncryption)
	ctOutRGSW, ctOutIsRGSW := any(ctOut).(*RGSWEncryption)

	switch {
	case ctOutIsRLWE:
		enc.encryptSecretKeyTo(ctOutRLWE, pt, outForm)

	case ctOutIsGadEnc:
		checkCiphertextShape(enc.params, ctOutGadEnc)

		if pt == nil {
			for i := 0; i < enc.params.GadgetLen(); i++ {
				enc.EncryptTo(ctOutGadEnc.Value[i], pt, outForm)
			}
			return
		}

		checkElementShape(enc.params, pt)

		var ptBufFull, ptCopyFull *Element
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

		if pt.AuxModLen() == 0 {
			ptCopy.WithModLen(ptCopy.BaseModLen(), 0).CopyFrom(pt)
			for i := 0; i < ptCopy.AuxModLen(); i++ {
				clear(ptCopy.Value.Coeffs[i])
			}
			ptCopy.SetForm(pt.Form())
		} else {
			ptCopy.CopyFrom(pt)
		}

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

		checkElementShape(enc.params, pt)

		var ptBufFull, ptCopyFull *Element
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

		if pt.AuxModLen() == 0 {
			ptCopy.WithModLen(ptCopy.BaseModLen(), 0).CopyFrom(pt)
			for i := 0; i < ptCopy.AuxModLen(); i++ {
				clear(ptCopy.Value.Coeffs[i])
			}
			ptCopy.SetForm(pt.Form())
		} else {
			ptCopy.CopyFrom(pt)
		}

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

func (enc *Encryptor) encryptSecretKeyTo(ctOut *Ciphertext, pt *Element, outForm dft.Form) {
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

	enc.maskSampler.SampleTo(ctOut.Mask.Value, ctMod)
	ctOut.Mask.SetForm(dft.FormNTT)

	ctOut.Body.Clear()
	ctOut.Body.SetForm(dft.FormNTT)
	enc.plainOp.MulSubTo(ctOut.Body, ctOut.Mask, enc.sk.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen()))

	noiseFull := enc.poolPoly.Get()
	defer enc.poolPoly.Put(noiseFull)
	noise := noiseFull.WithModLen(ctOut.BaseModLen(), ctOut.AuxModLen())

	enc.noiseSampler.SampleTo(noise.Value, ctMod)

	if ptBuf != nil {
		if ptBuf.Type() == crt.TypePoly && ptBuf.Form() == dft.FormNTT {
			enc.plainOp.AddTo(ctOut.Body, ctOut.Body, ptBuf)
		} else {
			enc.plainOp.AddTo(noise, noise, ptBuf)
		}
	}

	switch {
	case outForm == dft.FormCoeff:
		enc.plainOp.InvNTTTo(ctOut.Mask, ctOut.Mask)
		enc.plainOp.InvNTTTo(ctOut.Body, ctOut.Body)
	case outForm == dft.FormNTT:
		enc.plainOp.FwdNTTTo(noise, noise)
	}

	enc.plainOp.AddTo(ctOut.Body, ctOut.Body, noise)
}

// Phase returns ct.Body + ct.Mask * sk.
// Output is always in Coeff form.
func (enc *Encryptor) Phase(ct *Ciphertext) *Element {
	pOut := NewPolyCustom(enc.params.Rank(), ct.BaseModLen(), ct.AuxModLen(), dft.FormCoeff)
	enc.PhaseTo(pOut, ct)
	return pOut
}

// PhaseTo computes pOut = ct.Body + ct.Mask * sk.
// Output is always in Coeff form.
func (enc *Encryptor) PhaseTo(pOut *Element, ct *Ciphertext) {
	checkElementShape(enc.params, pOut)
	checkCiphertextShape(enc.params, ct)
	if pOut.Type() != crt.TypePoly {
		panic("inconsistent types")
	}

	ctNTTFull := enc.poolRLWE.Get()
	defer enc.poolRLWE.Put(ctNTTFull)
	ctNTT := ctNTTFull.WithModLen(ct.BaseModLen(), ct.AuxModLen())

	switch {
	case ct.Form() == dft.FormCoeff:
		enc.plainOp.FwdNTTTo(ctNTT.Body, ct.Body)
		enc.plainOp.FwdNTTTo(ctNTT.Mask, ct.Mask)
	case ct.Form() == dft.FormNTT:
		ctNTT.CopyFrom(ct)
	}

	pOut.CopyFrom(ctNTT.Body)
	enc.plainOp.MulAddTo(pOut, ctNTT.Mask, enc.sk.WithModLen(ctNTT.BaseModLen(), ctNTT.AuxModLen()))
	enc.plainOp.InvNTTTo(pOut, pOut)
}

// GenPublicKey generates a new [PublicKey].
// The output is always in NTT form.
func (enc *Encryptor) GenPublicKey() *PublicKey {
	pk := NewPublicKey(enc.params, enc.params.AuxModLen() > 0, dft.FormCoeff)
	enc.EncryptTo((*Ciphertext)(pk), nil, dft.FormNTT)
	return pk
}

// GenRelinKey generates a new [RelinKey].
// The output is always in NTT form.
func (enc *Encryptor) GenRelinKey() *RelinKey {
	rlk := NewRelinKey(enc.params, dft.FormNTT)
	enc.EncryptTo((*GadgetEncryption)(rlk), nil, dft.FormNTT)
	for i := 0; i < enc.params.GadgetLen(); i++ {
		enc.plainOp.MulAddTo(rlk.Value[i].Mask, enc.params.gadgetVec[i], enc.sk)
	}
	return rlk
}

// GenAutKey generates a new [AutKey] corresponding to idx.
// The output is always in NTT form.
func (enc *Encryptor) GenAutKey(idx int) *AutKey {
	if !enc.plainOp.CanAut(idx) {
		panic("invalid automorphism index")
	}

	idxInv := 1
	if enc.params.RingParams().CycloIdx() != 1 {
		idxMod := enc.params.RingParams().CycloIdxMod()
		idxInv = int(num.Inv(num.Reduce(idx, idxMod), idxMod))
	}

	skAut := enc.poolPoly.Get()
	defer enc.poolPoly.Put(skAut)
	enc.plainOp.AutTo((*Element)(skAut), (*Element)(enc.sk), idxInv)

	atk := NewAutKey(enc.params, idx, dft.FormNTT)
	for i := 0; i < enc.params.GadgetLen(); i++ {
		enc.maskSampler.SampleTo(atk.Value.Value[i].Mask.Value, enc.params.fullMod)
		atk.Value.Value[i].Mask.SetForm(dft.FormNTT)

		enc.noiseSampler.SampleTo(atk.Value.Value[i].Body.Value, enc.params.fullMod)

		enc.plainOp.FwdNTTTo(atk.Value.Value[i].Body, atk.Value.Value[i].Body)
		enc.plainOp.MulSubTo(atk.Value.Value[i].Body, atk.Value.Value[i].Mask, skAut)
		enc.plainOp.MulAddTo(atk.Value.Value[i].Body, enc.params.gadgetVec[i], enc.sk)
	}

	return atk
}

// GenKeySwitchKey generates a new [KeySwitchKey] corresponding to switching skIn -> sk.
// The output is always in NTT form.
func (enc *Encryptor) GenKeySwitchKey(skIn *SecretKey) *KeySwitchKey {
	ksk := NewKeySwitchKey(enc.params, dft.FormNTT)
	enc.EncryptTo((*GadgetEncryption)(ksk), (*Element)(skIn), dft.FormNTT)
	return ksk
}
