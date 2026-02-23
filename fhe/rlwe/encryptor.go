package rlwe

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/vec"
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

	ptPool *sync.Pool
	ctPool *sync.Pool
}

// NewEncryptor creates a new [Encryptor].
// Secret key is automatically sampled.
func NewEncryptor(params Parameters) *Encryptor {
	skSampler := params.secretKeyParams.Sampler()
	sk := &SecretKey{
		Value:  skSampler.Sample(params.RingParams().Rank(), params.fullMod),
		hasAux: params.HasAuxModulus(),
	}
	params.crtOp.FwdNTTTo(sk.Value, sk.Value)

	return newEncryptorWithKey(params, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] from the given secret key.
func NewEncryptorWithKey(params Parameters, sk *SecretKey) *Encryptor {
	checkShape(len(params.fullMod), params.HasAuxModulus(), (*Element)(sk))

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
	ctOut := NewCiphertext(e.params, hasAux, false)
	e.SampleRLWETo(ctOut, isNTT)
	return ctOut
}

// SampleRLWETo samples a new RLWE encryption of zero to ctOut.
func (e *Encryptor) SampleRLWETo(ctOut *Ciphertext, isNTT bool) {
	var modStart, modEnd int
	if ctOut.HasAuxModulus() {
		modStart, modEnd = 0, ctOut.ModLen()
	} else {
		modStart, modEnd = len(e.params.auxMod), len(e.params.auxMod)+ctOut.ModLen()
	}

	e.uSampler.SampleTo(ctOut.Mask.Value, e.params.fullMod[modStart:modEnd])
	ctOut.Mask.Value.IsNTT = true

	e.eSampler.SampleTo(ctOut.Body.Value, e.params.fullMod[modStart:modEnd])
	e.plainOp.FwdNTTTo(ctOut.Body, ctOut.Body)

	skValue := (*Element)(e.sk.WithModIdx(vec.Range(modStart, modEnd)...))
	e.plainOp.MulSubTo(ctOut.Body, skValue, ctOut.Mask)

	if !isNTT {
		e.plainOp.InvNTTTo(ctOut.Body, ctOut.Body)
		e.plainOp.InvNTTTo(ctOut.Mask, ctOut.Mask)
	}
}

// Encrypt encrypts pt.
func (e *Encryptor) Encrypt(pt *Element, isNTT bool) *Ciphertext {
	modLen := pt.ModLen()
	if e.params.HasAuxModulus() {
		modLen -= len(e.params.auxMod)
	}

	ctOut := NewCiphertextCustom(e.params.RingParams().Rank(), modLen, pt.HasAuxModulus(), false)
	e.EncryptTo(ctOut, pt, isNTT)
	return ctOut
}

// EncryptTo encrypts pt to ctOut.
func (e *Encryptor) EncryptTo(ctOut *Ciphertext, pt *Element, isNTT bool) {
	e.SampleRLWETo(ctOut, isNTT)

	if pt.IsNTT() == isNTT {
		e.plainOp.AddTo(ctOut.Body, ctOut.Body, pt)
	} else {
		ptBuf := e.ptPool.Get().(*Element)
		defer e.ptPool.Put(ptBuf)

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
	gadLen := e.params.GadgetLen()
	modLen := pt.ModLen()
	if e.params.HasAuxModulus() {
		modLen -= len(e.params.auxMod)
	}

	gOut := NewGadgetEncryptionCustom(e.params.RingParams().Rank(), modLen, e.params.HasAuxModulus(), gadLen, true)
	e.GadgetEncryptTo(gOut, pt, isNTT)
	return gOut
}

// GadgetEncryptTo encrypts pt to ctOut.
func (e *Encryptor) GadgetEncryptTo(ctOut *GadgetEncryption, pt *Element, isNTT bool) {
	if ctOut.GadgetLen() != e.params.GadgetLen() || ctOut.ModLen() != pt.ModLen() {
		panic("inconsistent output")
	}

	ptBuf := e.ptPool.Get().(*Element)
	defer e.ptPool.Put(ptBuf)

	if pt.IsNTT() == isNTT {
		ptBuf.CopyFrom(pt)
	} else {
		if isNTT {
			e.plainOp.FwdNTTTo(ptBuf, pt)
		} else {
			e.plainOp.InvNTTTo(ptBuf, pt)
		}
	}

	gadVec := e.dcmp.GadgetVector()
	for i := 0; i < e.params.GadgetLen(); i++ {
		e.SampleRLWETo(ctOut.Value[i], isNTT)
		e.plainOp.MulAddTo(ctOut.Value[i].Body, ptBuf, gadVec[i])
	}
}

// RGSWencryptTo encrypts pt to ctOut.
func (e *Encryptor) RGSWEncryptTo(ctOut *RGSW, pt *Element, isNTT bool) {
	if ctOut.GadgetLen() != e.params.GadgetLen() || ctOut.ModLen() != pt.ModLen() {
		panic("inconsistent output")
	}

	ptBuf := e.ptPool.Get().(*Element)
	defer e.ptPool.Put(ptBuf)

	if pt.IsNTT() == isNTT {
		ptBuf.CopyFrom(pt)
	} else {
		if isNTT {
			e.plainOp.FwdNTTTo(ptBuf, pt)
		} else {
			e.plainOp.InvNTTTo(ptBuf, pt)
		}
	}

	gadVec := e.dcmp.GadgetVector()
	for i := 0; i < e.params.GadgetLen(); i++ {
		e.SampleRLWETo(ctOut.Body.Value[i], isNTT)
		e.plainOp.MulAddTo(ctOut.Body.Value[i].Body, ptBuf, gadVec[i])

		e.SampleRLWETo(ctOut.Mask.Value[i], isNTT)
		e.plainOp.MulAddTo(ctOut.Mask.Value[i].Mask, ptBuf, gadVec[i])
	}
}

// // NewRelinKey creates a new relinearisation key.
// func (e *Encryptor) NewRelinKey() *RelinKey {
// 	rlk := NewGadgetEncryption(e.params, true)

// 	gParams := e.params.gadgetParams
// 	gLen := gParams.gadgetLen(e.params)
// 	modLen := len(e.params.modulus)
// 	if e.params.auxModulus != nil {
// 		modLen += len(e.params.auxModulus)
// 	}

// 	eval := e.op.PlainOp.SubEvaluatorAt(rlk.Value[0].HasAux, modLen)
// 	gadVec := e.op.Decmp.GadgetVector()

// 	for i := 0; i < gLen; i++ {
// 		e.SampleRlweTo(rlk.Value[i], true)
// 		eval.ScalarMulAddTo(rlk.Value[i].Mask, e.sk.Value, gadVec[i])
// 	}

// 	return &RelinKey{
// 		Value: rlk.Value,
// 	}
// }

// // NewKeySwitchKey creates a new key switch key.
// func (e *Encryptor) NewKeySwitchKey(skNew *SecretKey) *KeySwitchKey {
// 	pt := &PlainPoly{
// 		Value:  skNew.Value,
// 		HasAux: skNew.HasAux,
// 	}
// 	ksk := e.GadgetEncrypt(pt, true)

// 	return &KeySwitchKey{
// 		Value: ksk.Value,
// 	}
// }

// // NewAutomorphismKey creates a new automorphism key.
// func (e *Encryptor) NewAutomorphismKey(idx int) *AutomorphismKey {
// 	if e.sk == nil {
// 		panic("secret key is not set")
// 	}

// 	gParams := e.params.gadgetParams
// 	gLen := gParams.gadgetLen(e.params)
// 	modLen := len(e.params.modulus)
// 	var auxLen int
// 	if e.params.auxModulus != nil {
// 		auxLen = len(e.params.auxModulus)
// 	}

// 	atkVal := make([]*Ciphertext, gLen)
// 	for i := 0; i < gLen; i++ {
// 		atkVal[i] = NewCiphertextCustom(e.params.ringParams.Rank(), modLen, auxLen, true)
// 	}

// 	eval := e.op.PlainOp.SubEvaluatorAt(atkVal[0].HasAux, modLen+auxLen)
// 	idxInv := int(num.Inv(uint64(idx), num.NewModulus(e.params.ringParams.CycloOrder())))
// 	skAut := e.buf.pEnc
// 	eval.AutTo(skAut, e.sk.Value, idxInv)

// 	gadVec := e.op.Decmp.GadgetVector()
// 	for i := 0; i < gLen; i++ {
// 		e.uSampler.SampleTo(atkVal[i].Mask, eval.Modulus())
// 		atkVal[i].Mask.IsNTT = true

// 		e.eSampler.SampleTo(atkVal[i].Body, eval.Modulus())
// 		atkVal[i].Body.IsNTT = false

// 		eval.FwdNTTTo(atkVal[i].Body, atkVal[i].Body)
// 		eval.MulSubTo(atkVal[i].Body, skAut, atkVal[i].Mask)
// 		eval.ScalarMulAddTo(atkVal[i].Body, e.sk.Value, gadVec[i])
// 	}

// 	return &AutomorphismKey{
// 		Value: atkVal,
// 		Idx:   idx,
// 	}
// }
