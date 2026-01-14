package rlwe

import (
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// encryptorBuffer is a buffer for the [Encryptor].
type encryptorBuffer struct {
	pSample *crt.Poly
	pEnc    *crt.Poly

	cPhase *Ciphertext
}

func newEncryptorBuffer(params Parameters) encryptorBuffer {
	return encryptorBuffer{
		pSample: crt.NewPoly(params.ringParams.Rank(), len(params.auxModulus)+len(params.modulus)),
		pEnc:    crt.NewPoly(params.ringParams.Rank(), len(params.auxModulus)+len(params.modulus)),

		cPhase: NewCiphertext(params, true, false),
	}
}

// Encryptor is a RLWE encryptor for secret key encryption.
type Encryptor struct {
	// params is the parameters.
	params Parameters
	// sk is the secret key.
	sk *SecretKey
	// pk is the public key.
	pk *PublicKey

	// uSampler is the uniform sampler.
	uSampler crt.Sampler
	// rSampler is the sampler for the temporary key.
	rSampler crt.Sampler
	// eSampler is the sampler for the error.
	eSampler crt.Sampler

	// op is the operator.
	op *Operator

	// buf is the buffer for the encryptor.
	buf encryptorBuffer
}

// NewEncryptor creates a new [Encryptor] from the given parameters.
func NewEncryptor(params Parameters) *Encryptor {
	var mod []*num.Modulus
	if params.auxModulus != nil {
		mod = append(params.auxModulus, params.modulus...)
	} else {
		mod = params.modulus
	}

	ss := params.secretKeyParams.Sampler()
	sk := &SecretKey{
		Value:  ss.Sample(params.ringParams.Rank(), mod),
		HasAux: params.auxModulus != nil,
	}

	us := crt.UniformSamplerParameters{}.Sampler()
	es := params.noiseParams.Sampler()

	op := NewOperator(params)
	if !sk.Value.IsNTT {
		op.PlainOp.Eval.FwdNTTTo(sk.Value, sk.Value)
	}

	return &Encryptor{
		params: params,
		sk:     sk,
		pk:     nil,

		uSampler: us,
		eSampler: es,

		op: op,

		buf: newEncryptorBuffer(params),
	}
}

// NewEncryptorWithSK creates a new [Encryptor] from the given secret key.
func NewEncryptorWithSK(params Parameters, sk *SecretKey) *Encryptor {
	rP := params.ringParams
	mod := params.modulus
	aux := params.auxModulus
	cmod := append(aux, mod...)

	if sk.Value.ModLen() != len(cmod) || sk.Value.Rank() != rP.Rank() {
		panic("inconsistent secret key")
	}
	skCopy := &SecretKey{
		Value:  sk.Value.Copy(),
		HasAux: sk.HasAux,
	}

	us := crt.UniformSamplerParameters{}.Sampler()
	es := params.noiseParams.Sampler()

	op := NewOperator(params)
	if !skCopy.Value.IsNTT {
		op.PlainOp.Eval.FwdNTTTo(skCopy.Value, skCopy.Value)
	}

	return &Encryptor{
		params: params,
		sk:     skCopy,
		pk:     nil,

		uSampler: us,
		eSampler: es,

		op: op,

		buf: newEncryptorBuffer(params),
	}
}

// NewEncryptorWithPK creates a new [Encryptor] from the given public key.
func NewEncryptorWithPK(params Parameters, pk *PublicKey) *Encryptor {
	rP := params.ringParams
	mod := params.modulus
	aux := params.auxModulus
	cmod := append(aux, mod...)

	if pk.Body.ModLen() != len(cmod) || pk.Body.Rank() != rP.Rank() {
		panic("inconsistent public key")
	}
	pkCopy := &PublicKey{
		Body:   pk.Body.Copy(),
		Mask:   pk.Mask.Copy(),
		HasAux: pk.HasAux,
	}

	rs := params.secretKeyParams.Sampler()
	es := params.noiseParams.Sampler()

	op := NewOperator(params)
	if !pkCopy.Body.IsNTT {
		op.Eval.FwdNTTTo(pkCopy.Body, pkCopy.Body)
	}
	if !pkCopy.Mask.IsNTT {
		op.Eval.FwdNTTTo(pkCopy.Mask, pkCopy.Mask)
	}

	return &Encryptor{
		params: params,
		sk:     nil,
		pk:     pkCopy,

		rSampler: rs,
		eSampler: es,

		op: op,

		buf: newEncryptorBuffer(params),
	}
}

// NewPublicKey creates a new public key.
func (e *Encryptor) NewPublicKey() *PublicKey {
	var pk *Ciphertext
	if e.params.auxModulus == nil {
		modLen := len(e.params.modulus)
		pk = e.SampleRlweCustom(modLen, false, true)
	} else {
		modLen := len(e.params.auxModulus) + len(e.params.modulus)
		pk = e.SampleRlweCustom(modLen, true, true)
	}

	return &PublicKey{
		Body:   pk.Body,
		Mask:   pk.Mask,
		HasAux: pk.HasAux,
	}
}

// SecretKey returns the secret key.
func (e *Encryptor) SecretKey() *SecretKey {
	return e.sk
}

// PublicKey returns the public key.
func (e *Encryptor) PublicKey() *PublicKey {
	return e.pk
}

func (e *Encryptor) SampleRlwe(params Parameters, hasAux bool, isNTT bool) *Ciphertext {
	c := NewCiphertext(params, hasAux, isNTT)
	e.SampleRlweTo(c, isNTT)
	return c
}

// SampleRlweCustom samples a new RLWE ciphertext with the given modulus length, auxiliary flag, and NTT flag.
func (e *Encryptor) SampleRlweCustom(modLen int, hasAux bool, isNTT bool) *Ciphertext {
	var auxLen int
	if hasAux {
		if e.params.auxModulus == nil {
			panic("Auxiliary modulus is not set.")
		}
		auxLen = len(e.params.auxModulus)
		modLen -= auxLen
	}
	c := NewCiphertextCustom(e.params.ringParams.Rank(), modLen, auxLen, isNTT)
	e.SampleRlweTo(c, isNTT)
	return c
}

// SampleRlweTo samples a new RLWE ciphertext and stores the result in cOut.
func (e *Encryptor) SampleRlweTo(cOut *Ciphertext, isNTT bool) {
	rank := cOut.Body.Rank()
	modLen := cOut.Body.ModLen()
	if e.params.auxModulus != nil {
		modLen -= len(e.params.auxModulus)
	}

	if rank != e.params.ringParams.Rank() || modLen > len(e.params.modulus) {
		panic("inconsistent input ciphertext")
	} else if e.params.auxModulus == nil && cOut.HasAux {
		panic("auxiliary modulus is not set")
	}

	if e.sk != nil {
		e.sampleRlweToWithSK(cOut, isNTT)
	} else {
		e.sampleRlweToWithPK(cOut, isNTT)
	}
}

// sampleRlweToWithSK samples a new RLWE ciphertext with the secret key and stores the result in cOut.
func (e *Encryptor) sampleRlweToWithSK(cOut *Ciphertext, isNTT bool) {
	modLen := cOut.Body.ModLen()
	eval := e.op.PlainOp.SubEvaluatorAt(cOut.HasAux, modLen)
	mod := eval.Modulus()

	var begin, end int
	if cOut.HasAux {
		begin = 0
		end = modLen
	} else {
		begin = len(e.params.auxModulus)
		end = begin + modLen
	}
	skMod := e.sk.Value.WithModIdx(vec.Range(begin, end)...)

	e.uSampler.SampleTo(cOut.Mask, mod)
	cOut.Mask.IsNTT = true

	e.eSampler.SampleTo(cOut.Body, mod)
	cOut.Body.IsNTT = false

	eval.FwdNTTTo(cOut.Body, cOut.Body)
	eval.MulSubTo(cOut.Body, skMod, cOut.Mask)

	if !isNTT {
		eval.InvNTTTo(cOut.Body, cOut.Body)
		eval.InvNTTTo(cOut.Mask, cOut.Mask)
	}
}

// sampleRlweToWithPK samples a new RLWE ciphertext with the public key and stores the result in cOut.
func (e *Encryptor) sampleRlweToWithPK(cOut *Ciphertext, isNTT bool) {
	modLen := cOut.Body.ModLen()
	eval := e.op.PlainOp.SubEvaluatorAt(cOut.HasAux, modLen)
	mod := eval.Modulus()

	var begin, end int
	if cOut.HasAux {
		begin = 0
		end = modLen
	} else {
		begin = len(e.params.auxModulus)
		end = begin + modLen
	}
	pMod := e.buf.pSample.WithModIdx(vec.Range(begin, end)...)
	pkMod := &Ciphertext{
		Body: e.pk.Body.WithModIdx(vec.Range(begin, end)...),
		Mask: e.pk.Mask.WithModIdx(vec.Range(begin, end)...),
	}

	e.rSampler.SampleTo(pMod, mod)
	pMod.IsNTT = false
	eval.FwdNTTTo(pMod, pMod)

	e.eSampler.SampleTo(cOut.Mask, mod)
	cOut.Mask.IsNTT = false
	eval.FwdNTTTo(cOut.Mask, cOut.Mask)
	eval.MulAddTo(cOut.Mask, pMod, pkMod.Mask)

	e.eSampler.SampleTo(cOut.Body, mod)
	cOut.Body.IsNTT = false
	eval.FwdNTTTo(cOut.Body, cOut.Body)
	eval.MulAddTo(cOut.Body, pMod, pkMod.Body)

	if !isNTT {
		eval.InvNTTTo(cOut.Body, cOut.Body)
		eval.InvNTTTo(cOut.Mask, cOut.Mask)
	}
}

// Phase performs Phase(c).
func (e *Encryptor) Phase(c *Ciphertext) *PlainPoly {
	modLen := c.Body.ModLen()
	var auxLen int
	if c.HasAux {
		if e.params.auxModulus == nil {
			panic("auxiliary modulus is not set")
		}

		auxLen = len(e.params.auxModulus)
		modLen -= auxLen
	}

	pOut := NewPlainPolyCustom(e.params.ringParams.Rank(), modLen, auxLen, false)
	e.PhaseTo(pOut, c)
	return pOut
}

// PhaseTo performs Phase(c) and stores the result in pOut.
func (e *Encryptor) PhaseTo(pOut *PlainPoly, c *Ciphertext) {
	if e.sk == nil {
		panic("secret key is not set")
	} else if !c.IsCompatible(pOut) {
		panic("inconsistent output plaintext")
	} else if !c.CheckSanity() {
		panic("inconsistent ciphertext")
	}

	modLen := c.Body.ModLen()
	eval := e.op.PlainOp.SubEvaluatorAt(c.HasAux, modLen)

	var begin, end int
	if c.HasAux {
		begin = 0
		end = modLen
	} else {
		begin = len(e.params.auxModulus)
		end = begin + modLen
	}

	cPhase := e.buf.cPhase.WithModIdx(vec.Range(begin, end)...)
	keyMod := e.sk.Value.WithModIdx(vec.Range(begin, end)...)

	if c.Body.IsNTT {
		cPhase.Body.CopyFrom(c.Body)
		cPhase.Mask.CopyFrom(c.Mask)
	} else {
		eval.FwdNTTTo(cPhase.Body, c.Body)
		eval.FwdNTTTo(cPhase.Mask, c.Mask)
	}

	eval.MulAddTo(cPhase.Body, keyMod, cPhase.Mask)
	eval.InvNTTTo(pOut.Value, cPhase.Body)
	pOut.HasAux = c.HasAux
}

// ScalarEncrypt performs ScalarEncrypt(p).
func (e *Encryptor) ScalarEncrypt(s *PlainScalar, isNTT bool) *Ciphertext {
	modLen := s.ModLen()
	var auxLen int
	if s.HasAux {
		if e.params.auxModulus == nil {
			panic("auxiliary modulus is not set")
		}
		auxLen = len(e.params.auxModulus)
		modLen -= auxLen
	}

	c := NewCiphertextCustom(e.params.ringParams.Rank(), modLen, auxLen, false)
	e.ScalarEncryptTo(c, s, isNTT)
	return c
}

// ScalarEncryptTo performs ScalarEncrypt(s) and stores the result in cOut.
func (e *Encryptor) ScalarEncryptTo(cOut *Ciphertext, s *PlainScalar, isNTT bool) {
	if !cOut.IsCompatibleScalar(s) {
		panic("inconsistent output ciphertext")
	}

	eval := e.op.SubEvaluatorAt(s.HasAux, s.ModLen())
	e.SampleRlweTo(cOut, isNTT)
	eval.ScalarAddTo(cOut.Body, cOut.Body, s.Value)
}

// Encrypt performs Encrypt(p).
func (e *Encryptor) Encrypt(p *PlainPoly, isNTT bool) *Ciphertext {
	modLen := p.ModLen()
	var auxLen int
	if p.HasAux {
		if e.params.auxModulus == nil {
			panic("auxiliary modulus is not set")
		}
		auxLen = len(e.params.auxModulus)
		modLen -= auxLen
	}
	c := NewCiphertextCustom(e.params.ringParams.Rank(), modLen, auxLen, false)
	e.EncryptTo(c, p, isNTT)
	return c
}

// EncryptTo performs Encrypt(p) and stores the result in cOut.
func (e *Encryptor) EncryptTo(cOut *Ciphertext, p *PlainPoly, isNTT bool) {
	if !cOut.IsCompatible(p) {
		panic("inconsistent output ciphertext")
	}

	modLen := p.ModLen()
	eval := e.op.PlainOp.SubEvaluatorAt(p.HasAux, modLen)
	e.SampleRlweTo(cOut, isNTT)

	buf := e.buf.pEnc.WithModIdx(vec.Range(0, modLen)...)
	if isNTT && !p.Value.IsNTT {
		eval.FwdNTTTo(buf, p.Value)
	} else if !isNTT && p.Value.IsNTT {
		eval.InvNTTTo(buf, p.Value)
	} else {
		buf.CopyFrom(p.Value)
	}

	eval.AddTo(cOut.Body, buf, cOut.Body)
}

// ScalarGadgetEncrypt performs ScalarGadgetEncrypt(p).
func (e *Encryptor) ScalarGadgetEncrypt(p *PlainScalar, isNTT bool) *GadgetEncryption {
	g := NewGadgetEncryption(e.params, false)
	e.ScalarGadgetEncryptTo(g, p, isNTT)
	return g
}

// ScalarGadgetEncryptTo performs ScalarGadgetEncrypt(s) and stores the result in gOut.
func (e *Encryptor) ScalarGadgetEncryptTo(gOut *GadgetEncryption, p *PlainScalar, isNTT bool) {
	gParams := e.params.gadgetParams
	gLen := gParams.gadgetLen(e.params)
	modLen := len(e.params.modulus)
	if e.params.auxModulus != nil {
		modLen += len(e.params.auxModulus)
	}

	if gOut.GadgetLen() != gLen || gOut.ModLen() != modLen {
		panic("inconsistent output gadget encryption")
	} else if !gOut.Value[0].IsCompatibleScalar(p) {
		panic("inconsistent input plaintext")
	}

	eval := e.op.SubEvaluatorAt(p.HasAux, modLen)
	gadVec := e.op.Decmp.GadgetVector()
	mod := eval.Modulus()

	mulScalar := make([]uint64, modLen)
	for i := 0; i < gLen; i++ {
		e.SampleRlweTo(gOut.Value[i], isNTT)
		crt.MulScalarTo(mulScalar, p.Value, gadVec[i], mod)
		eval.ScalarAddTo(gOut.Value[i].Body, gOut.Value[i].Body, mulScalar)
	}
}

// GadgetEncrypt performs GadgetEncrypt(p).
func (e *Encryptor) GadgetEncrypt(p *PlainPoly, isNTT bool) *GadgetEncryption {
	g := NewGadgetEncryption(e.params, false)
	e.GadgetEncryptTo(g, p, isNTT)
	return g
}

// GadgetEncryptTo performs GadgetEncrypt(p) and stores the result in gOut.
func (e *Encryptor) GadgetEncryptTo(gOut *GadgetEncryption, p *PlainPoly, isNTT bool) {
	gParams := e.params.gadgetParams
	gLen := gParams.gadgetLen(e.params)
	modLen := len(e.params.modulus)
	if e.params.auxModulus != nil {
		modLen += len(e.params.auxModulus)
	}

	if gOut.GadgetLen() != gLen || gOut.ModLen() != modLen {
		panic("inconsistent output gadget encryption")
	} else if !gOut.Value[0].IsCompatible(p) {
		panic("inconsistent input plaintext")
	}

	eval := e.op.PlainOp.SubEvaluatorAt(p.HasAux, modLen)
	gadVec := e.op.Decmp.GadgetVector()

	buf := e.buf.pEnc
	if isNTT && !p.Value.IsNTT {
		eval.FwdNTTTo(buf, p.Value)
	} else if !isNTT && p.Value.IsNTT {
		eval.InvNTTTo(buf, p.Value)
	} else {
		buf.CopyFrom(p.Value)
	}

	for i := 0; i < gLen; i++ {
		e.SampleRlweTo(gOut.Value[i], isNTT)
		eval.ScalarMulAddTo(gOut.Value[i].Body, buf, gadVec[i])
	}
}

// ScalarRGSWEncrypt performs ScalarRGSWEncrypt(p).
func (e *Encryptor) ScalarRGSWEncrypt(p *PlainScalar, isNTT bool) *RGSW {
	r := NewRGSW(e.params, false)
	e.ScalarRGSWEncryptTo(r, p, isNTT)
	return r
}

// ScalarRGSWEncryptTo performs ScalarRGSWEncrypt(p) and stores the result in rOut.
func (e *Encryptor) ScalarRGSWEncryptTo(rOut *RGSW, p *PlainScalar, isNTT bool) {
	gParams := e.params.gadgetParams
	gLen := gParams.gadgetLen(e.params)
	modLen := len(e.params.modulus)
	if e.params.auxModulus != nil {
		modLen += len(e.params.auxModulus)
	}

	if rOut.GadgetLen() != gLen || rOut.ModLen() != modLen {
		panic("inconsistent output gadget encryption")
	} else if !rOut.Body.Value[0].IsCompatibleScalar(p) {
		panic("inconsistent input plaintext")
	}

	eval := e.op.SubEvaluatorAt(p.HasAux, modLen)
	gadVec := e.op.Decmp.GadgetVector()
	mod := eval.Modulus()

	mulScalar := make([]uint64, modLen)
	for i := 0; i < gLen; i++ {
		crt.MulScalarTo(mulScalar, p.Value, gadVec[i], mod)

		e.SampleRlweTo(rOut.Body.Value[i], isNTT)
		eval.ScalarAddTo(rOut.Body.Value[i].Body, rOut.Body.Value[i].Body, mulScalar)

		e.SampleRlweTo(rOut.Mask.Value[i], isNTT)
		eval.ScalarAddTo(rOut.Mask.Value[i].Mask, rOut.Mask.Value[i].Mask, mulScalar)
	}
}

// RGSWEncrypt performs RGSWEncrypt(p).
func (e *Encryptor) RGSWEncrypt(p *PlainPoly, isNTT bool) *RGSW {
	r := NewRGSW(e.params, false)
	e.RGSWEncryptTo(r, p, isNTT)
	return r
}

// RGSWEncryptTo performs RGSWEncrypt(p) and stores the result in rOut.
func (e *Encryptor) RGSWEncryptTo(rOut *RGSW, p *PlainPoly, isNTT bool) {
	gParams := e.params.gadgetParams
	gLen := gParams.gadgetLen(e.params)
	modLen := len(e.params.modulus)
	if e.params.auxModulus != nil {
		modLen += len(e.params.auxModulus)
	}

	if rOut.GadgetLen() != gLen || rOut.ModLen() != modLen {
		panic("inconsistent output gadget encryption")
	} else if !rOut.Body.Value[0].IsCompatible(p) {
		panic("inconsistent input plaintext")
	}

	eval := e.op.PlainOp.SubEvaluatorAt(p.HasAux, modLen)
	gadVec := e.op.Decmp.GadgetVector()

	buf := e.buf.pEnc
	if isNTT && !p.Value.IsNTT {
		eval.FwdNTTTo(buf, p.Value)
	} else if !isNTT && p.Value.IsNTT {
		eval.InvNTTTo(buf, p.Value)
	} else {
		buf.CopyFrom(p.Value)
	}

	for i := 0; i < gLen; i++ {
		e.SampleRlweTo(rOut.Body.Value[i], isNTT)
		eval.ScalarMulAddTo(rOut.Body.Value[i].Body, buf, gadVec[i])

		e.SampleRlweTo(rOut.Mask.Value[i], isNTT)
		eval.ScalarMulAddTo(rOut.Mask.Value[i].Mask, buf, gadVec[i])
	}
}

// NewRelinKey creates a new relinearisation key.
func (e *Encryptor) NewRelinKey() *RelinKey {
	rlk := NewGadgetEncryption(e.params, true)

	gParams := e.params.gadgetParams
	gLen := gParams.gadgetLen(e.params)
	modLen := len(e.params.modulus)
	if e.params.auxModulus != nil {
		modLen += len(e.params.auxModulus)
	}

	eval := e.op.PlainOp.SubEvaluatorAt(rlk.Value[0].HasAux, modLen)
	gadVec := e.op.Decmp.GadgetVector()

	for i := 0; i < gLen; i++ {
		e.SampleRlweTo(rlk.Value[i], true)
		eval.ScalarMulAddTo(rlk.Value[i].Mask, e.sk.Value, gadVec[i])
	}

	return &RelinKey{
		Value: rlk.Value,
	}
}

// NewKeySwitchKey creates a new key switch key.
func (e *Encryptor) NewKeySwitchKey(skNew *SecretKey) *KeySwitchKey {
	pt := &PlainPoly{
		Value:  skNew.Value,
		HasAux: skNew.HasAux,
	}
	ksk := e.GadgetEncrypt(pt, true)

	return &KeySwitchKey{
		Value: ksk.Value,
	}
}

// NewAutomorphismKey creates a new automorphism key.
func (e *Encryptor) NewAutomorphismKey(idx int) *AutomorphismKey {
	if e.sk == nil {
		panic("secret key is not set")
	}

	gParams := e.params.gadgetParams
	gLen := gParams.gadgetLen(e.params)
	modLen := len(e.params.modulus)
	var auxLen int
	if e.params.auxModulus != nil {
		auxLen = len(e.params.auxModulus)
	}

	atkVal := make([]*Ciphertext, gLen)
	for i := 0; i < gLen; i++ {
		atkVal[i] = NewCiphertextCustom(e.params.ringParams.Rank(), modLen, auxLen, true)
	}

	eval := e.op.PlainOp.SubEvaluatorAt(atkVal[0].HasAux, modLen+auxLen)
	idxInv := int(num.Inv(uint64(idx), num.NewModulus(e.params.ringParams.CycloOrder())))
	skAut := e.buf.pEnc
	eval.AutTo(skAut, e.sk.Value, idxInv)

	gadVec := e.op.Decmp.GadgetVector()
	for i := 0; i < gLen; i++ {
		e.uSampler.SampleTo(atkVal[i].Mask, eval.Modulus())
		atkVal[i].Mask.IsNTT = true

		e.eSampler.SampleTo(atkVal[i].Body, eval.Modulus())
		atkVal[i].Body.IsNTT = false

		eval.FwdNTTTo(atkVal[i].Body, atkVal[i].Body)
		eval.MulSubTo(atkVal[i].Body, skAut, atkVal[i].Mask)
		eval.ScalarMulAddTo(atkVal[i].Body, e.sk.Value, gadVec[i])
	}

	return &AutomorphismKey{
		Value: atkVal,
		Idx:   idx,
	}
}
