package heint

import (
	"github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

// Encryptor encrypts/decrypts [*Ciphertext] and [*Plaintext].
type Encryptor struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	pOp  *rlwe.PlainOperator
	enc  *rlwe.Encryptor
	pack pack.IntPacker
	ecd  *Encoder

	scFacs []*rlwe.Element

	ePool *rlwe.ElementPool
}

// NewEncryptor creates a new [Encryptor].
func NewEncryptor(params rlwe.Parameters, msgMod *num.Modulus) *Encryptor {
	skSampler := params.SecretKeyParams().Sampler()
	skValue := skSampler.Sample(params.Rank(), params.FullModulus())
	sk := (*rlwe.SecretKey)(rlwe.NewElementFrom(skValue, len(params.AuxModulus())))
	params.Operator().FwdNTTTo(skValue, skValue)

	return NewEncryptorWithKey(params, msgMod, sk)
}

// NewEncryptorWithKey creates a new [Encryptor] with a secret key.
func NewEncryptorWithKey(params rlwe.Parameters, msgMod *num.Modulus, skNTT *rlwe.SecretKey) *Encryptor {
	baseMod := params.BaseModulus()
	scFacs := computeScalingFactor(baseMod, msgMod)

	return &Encryptor{
		params: params,
		msgMod: msgMod,

		pOp:  rlwe.NewPlainOperator(params),
		enc:  rlwe.NewEncryptorWithKey(params, skNTT),
		pack: pack.NewIntPacker(params.RingParams(), msgMod),
		ecd:  NewEncoder(params, msgMod),

		scFacs: scFacs,

		ePool: rlwe.NewElementPool(params, params.HasAuxModulus(), true),
	}
}

// Parameters returns the parameters.
func (e *Encryptor) Parameters() rlwe.Parameters {
	return e.params
}

// SecretKey returns the secret key.
func (e *Encryptor) SecretKey() *rlwe.SecretKey {
	return e.enc.SecretKey()
}

// NewRelinKey creates a new relinearisation key.
func (e *Encryptor) NewRelinKey() *rlwe.RelinKey {
	return e.enc.NewRelinKey()
}

// NewKeySwitchKey creates a new key switch key.
func (e *Encryptor) NewKeySwitchKey(skNew *rlwe.SecretKey) *rlwe.KeySwitchKey {
	return e.enc.NewKeySwitchKey(skNew)
}

// NewAutomorphismKey creates a new automorphism key.
func (e *Encryptor) NewAutomorphismKey(idx int) *rlwe.AutomorphismKey {
	return e.enc.NewAutomorphismKey(idx)
}

// NewRotationKey creates a new rotation key.
func (e *Encryptor) NewRotationKey(idx []int) *rlwe.AutomorphismKey {
	autIdx := e.pack.RotIdxToAutIdx(idx)
	return e.enc.NewAutomorphismKey(autIdx)
}

// Encrypt encrypts the message v.
func (e *Encryptor) Encrypt(v []uint64, isNTT bool) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertext(e.params, false, isNTT)
	e.EncryptTo(ctOut, v, isNTT)
	return ctOut
}

// EncryptCustom encrypts the message v with custom parameters.
func (e *Encryptor) EncryptCustom(v []uint64, baseLen int, isNTT bool) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertextCustom(e.params.Rank(), baseLen, 0, isNTT)
	e.EncryptTo(ctOut, v, isNTT)
	return ctOut
}

// EncryptTo encrypts the message v to ctOut.
func (e *Encryptor) EncryptTo(ctOut *rlwe.Ciphertext, v []uint64, isNTT bool) {
	if ctOut.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	}

	baseLen := ctOut.BaseModLen()
	var pt *rlwe.Element
	switch len(v) {
	case 1:
		pt = e.ePool.Get(crt.TypeScalar)
	case ctOut.Rank():
		pt = e.ePool.Get(crt.TypePoly)
	default:
		panic("invalid input length")
	}
	defer e.ePool.Put(pt)
	pt = pt.WithModLen(baseLen, 0)

	e.ecd.ScaleEncodeTo(pt, v, isNTT)
	e.enc.EncryptTo(ctOut, pt, isNTT)
}

// EncryptElement encrypts the element e.
func (e *Encryptor) EncryptElement(eIn *rlwe.Element, isNTT bool) *rlwe.Ciphertext {
	ctOut := rlwe.NewCiphertextCustom(e.params.Rank(), eIn.BaseModLen(), 0, isNTT)
	e.EncryptElementTo(ctOut, eIn, isNTT)
	return ctOut
}

// EncryptElementTo encrypts the element e to ctOut.
func (e *Encryptor) EncryptElementTo(ctOut *rlwe.Ciphertext, eIn *rlwe.Element, isNTT bool) {
	if eIn.AuxModLen() > 0 {
		panic("auxiliary modulus length should be zero")
	}

	baseLen := eIn.BaseModLen()

	pt := e.ePool.Get(eIn.Type())
	defer e.ePool.Put(pt)
	pt = pt.WithModLen(baseLen, 0)

	pt.CopyFrom(eIn)
	e.pOp.MulTo(pt, pt, e.scFacs[baseLen-1])
	e.enc.EncryptTo(ctOut, pt, isNTT)
}

// Decrypt decrypts the ciphertext ct.
func (e *Encryptor) Decrypt(ct *rlwe.Ciphertext) []uint64 {
	vOut := make([]uint64, ct.Rank())
	e.DecryptTo(vOut, ct)
	return vOut
}

// DecryptTo decrypts the ciphertext ct to vOut.
func (e *Encryptor) DecryptTo(vOut []uint64, ct *rlwe.Ciphertext) {
	if len(vOut) != 1 && len(vOut) != ct.Rank() {
		panic("inconsistent output")
	}

	baseLen := ct.BaseModLen()

	var pt *rlwe.Element
	switch len(vOut) {
	case 1:
		pt = e.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		pt = e.ePool.Get(crt.TypePoly)
	}
	defer e.ePool.Put(pt)
	pt = pt.WithModLen(baseLen, 0)
	ptScale := pt.WithModLen(1, 0)

	e.PhaseTo(pt, ct)
	e.ecd.ScaleDecodeTo(ptScale.Value.Coeffs[0], pt)
	copy(vOut, ptScale.Value.Coeffs[0][:len(vOut)])
}

// Phase performs Phase(ct).
func (e *Encryptor) Phase(ct *rlwe.Ciphertext) *rlwe.Element {
	return e.enc.Phase(ct)
}

// PhaseTo performs Phase(ct) and stores the result in pt.
func (e *Encryptor) PhaseTo(eOut *rlwe.Element, ct *rlwe.Ciphertext) {
	e.enc.PhaseTo((*rlwe.Element)(eOut), ct)
}

// Noise returns the noise of the ciphertext.
func (e *Encryptor) Noise(ct *rlwe.Ciphertext) *rlwe.Element {
	noise := rlwe.NewElement(ct.Rank(), ct.BaseModLen(), 0, false)
	e.NoiseTo(noise, ct)
	return noise
}

// NoiseTo stores the noise of the ciphertext in eOut.
func (e *Encryptor) NoiseTo(eOut *rlwe.Element, ct *rlwe.Ciphertext) {
	if !eOut.IsConsistent(ct.Body) {
		panic("inconsistent output")
	}

	baseLen := ct.BaseModLen()

	pt := e.ePool.Get(crt.TypePoly)
	vEcd := e.ePool.Get(crt.TypePoly)
	defer e.ePool.Put(pt)
	defer e.ePool.Put(vEcd)
	pt = pt.WithModLen(baseLen, 0)
	vEcd = vEcd.WithModLen(baseLen, 0)
	v := vEcd.Value.WithModIdx(0)

	e.PhaseTo(pt, ct)
	e.ecd.ScaleDecodeTo(v.Coeffs[0], pt)
	e.ecd.ScaleEncodeTo(vEcd, v.Coeffs[0], false)
	e.pOp.SubTo(eOut, pt, vEcd)
}
