package rlwe

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encryptor encrypts/decrypts RLWE ciphertexts.
type Encryptor struct {
	params Parameters

	sk *SecretKey

	plainOp *PlainOperator
	crtOp   crt.Operator
	dcmp    Decomposer

	uSampler   crt.Sampler
	auxSampler crt.Sampler
	eSampler   crt.Sampler

	sPool  *sync.Pool
	ptPool *sync.Pool
	ctPool *sync.Pool
}

// NewEncryptor creates a new [Encryptor].
// Secret key is automatically sampled.
func NewEncryptor(params Parameters) *Encryptor {
	skSampler := params.secretKeyParams.Sampler()
	sk := &SecretKey{
		Value:  skSampler.Sample(params.Rank(), params.fullMod),
		auxLen: len(params.auxMod),
	}
	params.crtOp.FwdNTTTo(sk.Value, sk.Value)

	return newEncryptorWithKey(params, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] from the given secret key.
func NewEncryptorWithKey(params Parameters, sk *SecretKey) *Encryptor {
	checkShape(len(params.baseMod), len(params.auxMod), (*Element)(sk))

	skCopy := sk.Copy()
	if !skCopy.Value.IsNTT {
		params.crtOp.FwdNTTTo(skCopy.Value, skCopy.Value)
	}

	return newEncryptorWithKey(params, skCopy)
}

// newEncryptorWithKey creates a new [Encryptor] from the given secret key in NTT form, without copying.
func newEncryptorWithKey(params Parameters, skNTT *SecretKey) *Encryptor {
	return &Encryptor{
		params: params,

		sk: skNTT,

		uSampler: crt.UniformSamplerParameters{}.Sampler(),
		eSampler: params.noiseParams.Sampler(),

		plainOp: NewPlainOperator(params),
		crtOp:   params.crtOp,
		dcmp:    NewDecomposer(params),

		sPool: &sync.Pool{
			New: func() any {
				return NewScalar(params, params.HasAuxModulus())
			},
		},
		ptPool: &sync.Pool{
			New: func() any {
				return NewPoly(params, params.HasAuxModulus(), false)
			},
		},
		ctPool: &sync.Pool{
			New: func() any {
				return NewCiphertext(params, params.HasAuxModulus(), false)
			},
		},
	}
}

// SecretKey returns the secret key.
func (e *Encryptor) SecretKey() *SecretKey {
	return e.sk
}

// SampleRLWE returns a new RLWE encryption of zero.
func (e *Encryptor) SampleRLWE(hasAux, isNTT bool) *Ciphertext {
	ctOut := NewCiphertext(e.params, hasAux, isNTT)
	e.SampleRLWETo(ctOut, isNTT)
	return ctOut
}

// SampleRLWECustom samples a new RLWE encryption of zero with the given parameters.
func (e *Encryptor) SampleRLWECustom(baseLen, auxLen int, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(e.params.Rank(), baseLen, auxLen, isNTT)
	e.SampleRLWETo(ctOut, isNTT)
	return ctOut
}

// SampleRLWETo samples a new RLWE encryption of zero to ctOut.
func (e *Encryptor) SampleRLWETo(ctOut *Ciphertext, isNTT bool) {
	var modStart, modEnd int
	baseLen := ctOut.BaseModLen()
	auxLen := ctOut.AuxModLen()

	modStart = len(e.params.auxMod) - auxLen
	modEnd = len(e.params.auxMod) + baseLen

	e.uSampler.SampleTo(ctOut.Mask.Value, e.params.fullMod[modStart:modEnd])
	ctOut.Mask.Value.IsNTT = true

	e.eSampler.SampleTo(ctOut.Body.Value, e.params.fullMod[modStart:modEnd])
	e.plainOp.FwdNTTTo(ctOut.Body, ctOut.Body)

	skValue := (*Element)(e.sk.WithModLen(baseLen, auxLen))
	e.plainOp.MulSubTo(ctOut.Body, skValue, ctOut.Mask)

	if !isNTT {
		e.plainOp.InvNTTTo(ctOut.Body, ctOut.Body)
		e.plainOp.InvNTTTo(ctOut.Mask, ctOut.Mask)
	}
}

// Phase performs Phase(c).
func (e *Encryptor) Phase(c *Ciphertext) *Element {
	baseLen := c.BaseModLen()
	auxLen := c.AuxModLen()

	eOut := NewElement(e.params.Rank(), baseLen, auxLen, false)
	e.PhaseTo(eOut, c)
	return eOut
}

// PhaseTo performs Phase(c) and stores the result in pOut.
func (e *Encryptor) PhaseTo(eOut *Element, c *Ciphertext) {
	if e.sk == nil {
		panic("secret key is not set")
	}

	baseLen := c.BaseModLen()
	auxLen := c.AuxModLen()

	cPhase := e.ctPool.Get().(*Ciphertext)
	defer e.ctPool.Put(cPhase)
	cPhase = cPhase.WithModLen(baseLen, auxLen)

	keyMod := (*Element)(e.sk.WithModLen(baseLen, auxLen))

	if c.IsNTT() {
		cPhase.Body.CopyFrom(c.Body)
		cPhase.Mask.CopyFrom(c.Mask)
	} else {
		e.plainOp.FwdNTTTo(cPhase.Body, c.Body)
		e.plainOp.FwdNTTTo(cPhase.Mask, c.Mask)
	}

	e.plainOp.MulAddTo(cPhase.Body, keyMod, cPhase.Mask)
	e.plainOp.InvNTTTo(eOut, cPhase.Body)
}

// Encrypt encrypts pt.
func (e *Encryptor) Encrypt(pt *Element, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(e.params.Rank(), pt.BaseModLen(), pt.AuxModLen(), isNTT)
	e.EncryptTo(ctOut, pt, isNTT)
	return ctOut
}

// EncryptTo encrypts pt to ctOut.
func (e *Encryptor) EncryptTo(ctOut *Ciphertext, pt *Element, isNTT bool) {
	e.SampleRLWETo(ctOut, isNTT)

	if pt.Value.Type() == crt.TypeScalar || pt.IsNTT() == isNTT {
		e.plainOp.AddTo(ctOut.Body, ctOut.Body, pt)
	} else {
		ptBuf := e.ptPool.Get().(*Element)
		defer e.ptPool.Put(ptBuf)

		ptBuf = ptBuf.WithModLen(pt.BaseModLen(), pt.AuxModLen())

		if isNTT {
			e.plainOp.FwdNTTTo(ptBuf, pt)
		} else {
			e.plainOp.InvNTTTo(ptBuf, pt)
		}
		e.plainOp.AddTo(ctOut.Body, ctOut.Body, ptBuf)
	}
}

// GadgetEncrypt encrypts pt.
func (e *Encryptor) GadgetEncrypt(pt *Element, isNTT bool) *GadgetEncryption {
	gOut := NewGadgetEncryptionCustom(e.params.Rank(), pt.BaseModLen(), pt.AuxModLen(), e.params.GadgetLen(), true)
	e.GadgetEncryptTo(gOut, pt, isNTT)
	return gOut
}

// GadgetEncryptTo encrypts pt to ctOut.
func (e *Encryptor) GadgetEncryptTo(ctOut *GadgetEncryption, pt *Element, isNTT bool) {
	if pt.BaseModLen() != len(e.params.baseMod) || pt.AuxModLen() != len(e.params.auxMod) {
		panic("inconsistent input plaintext")
	} else if ctOut.GadgetLen() != e.params.GadgetLen() {
		panic("inconsistent output gadget encryption")
	}

	var ptBuf *Element
	var ptMul *Element

	if pt.Value.Type() == crt.TypeScalar {
		ptBuf = e.sPool.Get().(*Element)
		ptMul = e.sPool.Get().(*Element)
		defer e.sPool.Put(ptBuf)
		defer e.sPool.Put(ptMul)

		ptBuf.CopyFrom(pt)
	} else {
		ptBuf = e.ptPool.Get().(*Element)
		ptMul = e.ptPool.Get().(*Element)
		defer e.ptPool.Put(ptBuf)
		defer e.ptPool.Put(ptMul)

		if isNTT == pt.IsNTT() {
			ptBuf.CopyFrom(pt)
		} else if isNTT {
			e.plainOp.FwdNTTTo(ptBuf, pt)
		} else {
			e.plainOp.InvNTTTo(ptBuf, pt)
		}
	}

	gadVec := e.dcmp.GadgetVector()
	for i := 0; i < e.params.GadgetLen(); i++ {
		e.SampleRLWETo(ctOut.Value[i], isNTT)
		e.plainOp.MulTo(ptMul, ptBuf, gadVec[i])
		e.plainOp.AddTo(ctOut.Value[i].Body, ctOut.Value[i].Body, ptMul)
	}
}

// RGSWEncrypt encrypts pt.
func (e *Encryptor) RGSWEncrypt(pt *Element, isNTT bool) *RGSW {
	ctOut := NewRGSWCustom(e.params.Rank(), pt.BaseModLen(), pt.AuxModLen(), e.params.GadgetLen(), true)
	e.RGSWEncryptTo(ctOut, pt, isNTT)
	return ctOut
}

// RGSWencryptTo encrypts pt to ctOut.
func (e *Encryptor) RGSWEncryptTo(ctOut *RGSW, pt *Element, isNTT bool) {
	if pt.BaseModLen() != len(e.params.baseMod) || pt.AuxModLen() != len(e.params.auxMod) {
		panic("inconsistent input plaintext")
	} else if ctOut.GadgetLen() != e.params.GadgetLen() {
		panic("inconsistent output RGSW ciphertext")
	}

	var ptBuf *Element
	var ptMul *Element

	if pt.Value.Type() == crt.TypeScalar {
		ptBuf = e.sPool.Get().(*Element)
		ptMul = e.sPool.Get().(*Element)
		defer e.sPool.Put(ptBuf)
		defer e.sPool.Put(ptMul)

		ptBuf.CopyFrom(pt)
	} else {
		ptBuf = e.ptPool.Get().(*Element)
		ptMul = e.ptPool.Get().(*Element)
		defer e.ptPool.Put(ptBuf)
		defer e.ptPool.Put(ptMul)

		if isNTT == pt.IsNTT() {
			ptBuf.CopyFrom(pt)
		} else if isNTT {
			e.plainOp.FwdNTTTo(ptBuf, pt)
		} else {
			e.plainOp.InvNTTTo(ptBuf, pt)
		}
	}

	gadVec := e.dcmp.GadgetVector()
	for i := 0; i < e.params.GadgetLen(); i++ {
		e.plainOp.MulTo(ptMul, ptBuf, gadVec[i])

		e.SampleRLWETo(ctOut.Body.Value[i], isNTT)
		e.plainOp.AddTo(ctOut.Body.Value[i].Body, ctOut.Body.Value[i].Body, ptMul)

		e.SampleRLWETo(ctOut.Mask.Value[i], isNTT)
		e.plainOp.AddTo(ctOut.Mask.Value[i].Mask, ctOut.Mask.Value[i].Mask, ptMul)
	}
}

// NewRelinKey creates a new relinearisation key.
func (e *Encryptor) NewRelinKey() *RelinKey {
	rlk := NewGadgetEncryption(e.params, true)

	skValue := (*Element)(e.sk)
	for i := 0; i < e.params.GadgetLen(); i++ {
		e.SampleRLWETo(rlk.Value[i], true)
		e.plainOp.MulAddTo(rlk.Value[i].Mask, skValue, e.dcmp.GadgetVector()[i])
	}

	return (*RelinKey)(rlk)
}

// NewKeySwitchKey creates a new key switch key.
func (e *Encryptor) NewKeySwitchKey(skNew *SecretKey) *KeySwitchKey {
	return (*KeySwitchKey)(e.GadgetEncrypt((*Element)(skNew), true))
}

// NewAutomorphismKey creates a new automorphism key.
func (e *Encryptor) NewAutomorphismKey(idx int) *AutomorphismKey {
	if e.sk == nil {
		panic("secret key is not set")
	}

	gLen := e.dcmp.Params().GadgetLen()
	baseLen := len(e.params.baseMod)
	auxLen := len(e.params.auxMod)

	pOp := e.plainOp

	atkVal := make([]*Ciphertext, gLen)
	for i := 0; i < gLen; i++ {
		atkVal[i] = NewCiphertextCustom(e.params.Rank(), baseLen, auxLen, true)
	}

	skAut := e.ptPool.Get().(*Element)
	defer e.ptPool.Put(skAut)

	idxInv := int(num.Inv(uint64(idx), num.NewModulus(e.params.RingParams().CycloOrder())))
	pOp.AutTo(skAut, (*Element)(e.sk), idxInv)

	gadVec := e.dcmp.GadgetVector()
	for i := 0; i < gLen; i++ {
		e.uSampler.SampleTo(atkVal[i].Mask.Value, e.params.fullMod)
		atkVal[i].Mask.Value.IsNTT = true

		e.eSampler.SampleTo(atkVal[i].Body.Value, e.params.fullMod)
		atkVal[i].Body.Value.IsNTT = false

		pOp.FwdNTTTo(atkVal[i].Body, atkVal[i].Body)
		pOp.MulSubTo(atkVal[i].Body, skAut, atkVal[i].Mask)
		pOp.MulAddTo(atkVal[i].Body, (*Element)(e.sk), gadVec[i])
	}

	return &AutomorphismKey{
		Value: atkVal,
		Idx:   idx,
	}
}
