package crt

import (
	"math/big"
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/float128"
	"github.com/hienaa-org/hienaa/math/num"
)

// Embedder embeds a polynomial into different modulus.
// In other words, it computes
//
//	[p]_modIn -> [p]_modOut
//
// It uses HPS-like algorithm, so the computation is exact.
type Embedder struct {
	// modIn is the input modulus.
	modIn []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus

	// compInv is the modular inverse of the compliment of the input modulus limb.
	compInv []uint64
	// compInvS is the Shoup multiplication form of the modular inverse of the compliment of the input modulus limb.
	compInvS []uint64

	// comp is the compliment of the input modulus limb.
	comp [][]uint64
	// compS is the Shoup form of comp.
	compS [][]uint64

	// negMod is the negative of the input modulus modulo the output modulus limb.
	negMod []uint64
	// negModS is the Shoup form of negMod.
	negModS []uint64

	// inv is the floating-point approximation of the inverse of the input modulus limb.
	inv []float128.Float128

	// idx holds the index of the input modulus limb if it overlaps with the output modulus limb.
	// For example, if modOut[i] = modIn[j], then idx[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	idx []int

	buf embedderBuffer
}

// embedderBuffer is a buffer for [Embedder].
type embedderBuffer struct {
	// f128 is a buffer for the floating-point number.
	// Always has length 8.
	f128 []float128.Float128
	// u64 is a buffer for the 64-bit integers.
	// Always has length 8.
	u64 []uint64
	// in is a buffer for input coefficient.
	// Always has length [len(modIn)][8].
	in [][]uint64
}

// NewEmbedder creates a new [Embedder].
func NewEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *Embedder {
	if !isCoprime(modIn) || !isCoprime(modOut) {
		panic("modulus must be coprime")
	}

	compInv := make([]uint64, len(modIn))
	compInvS := make([]uint64, len(modIn))

	inv := make([]float128.Float128, len(modIn))

	invRat := big.NewRat(0, 1)
	for i := 0; i < len(modIn); i++ {
		compInv[i] = 1
		for j := 0; j < len(modIn); j++ {
			if i != j {
				compInv[i] = num.Mul(compInv[i], num.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
		}
		compInvS[i] = num.SForm(compInv[i], modIn[i])

		invRat.Denom().SetUint64(modIn[i].Value())
		invRat.Num().SetInt64(1)
		inv[i] = float128.FromRat(invRat)
	}

	comp := make([][]uint64, len(modOut))
	compS := make([][]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		comp[i] = make([]uint64, len(modIn))
		compS[i] = make([]uint64, len(modIn))
	}

	negMod := make([]uint64, len(modOut))
	negModS := make([]uint64, len(modOut))
	idx := make([]int, len(modOut))

	for i := 0; i < len(modOut); i++ {
		negMod[i] = modOut[i].Value() - 1
		idx[i] = -1

		for j := 0; j < len(modIn); j++ {
			comp[i][j] = 1
			for k := 0; k < len(modIn); k++ {
				if j != k {
					comp[i][j] = num.Mul(comp[i][j], modIn[k].Value(), modOut[i])
				}
			}
			compS[i][j] = num.SForm(comp[i][j], modOut[i])

			negMod[i] = num.Mul(negMod[i], modIn[j].Value(), modOut[i])

			if modOut[i].Value() == modIn[j].Value() {
				idx[i] = j
			}
		}

		negModS[i] = num.SForm(negMod[i], modOut[i])
	}

	return &Embedder{
		modIn:  modIn,
		modOut: modOut,

		compInv:  compInv,
		compInvS: compInvS,

		comp:  comp,
		compS: compS,

		negMod:  negMod,
		negModS: negModS,

		inv: inv,

		idx: idx,

		buf: newEmbedderBuffer(modIn),
	}
}

// newEmbedderBuffer creates a new [embedderBuffer].
func newEmbedderBuffer(modIn []*num.Modulus) embedderBuffer {
	in := make([][]uint64, len(modIn))
	for i := range in {
		in[i] = make([]uint64, 8)
	}

	return embedderBuffer{
		f128: make([]float128.Float128, 8),
		u64:  make([]uint64, 8),
		in:   in,
	}
}

// Embed returns the embedding of e to the output modulus.
// If e.ModLen() < len(emb.modIn), it only embeds the first e.ModLen() elements.
func (emb *Embedder) Embed(e *Element) *Element {
	eOut := NewPoly(e.Rank(), e.ModLen())
	emb.EmbedTo(eOut, e)
	return eOut
}

// EmbedTo embeds e to eOut.
// If e.ModLen() < len(emb.modIn) or eOut.ModLen() < len(emb.modOut),
// it only embeds the first emb.ModLen() elements to eOut.ModLen() elements.
func (emb *Embedder) EmbedTo(eOut, e *Element) {
	checkBinaryOperable(e.Rank(), e.ModLen(), eOut, e)
	if e.Type() == TypePoly && e.IsNTT {
		panic("input(s) must be in standard form")
	}

	emb.EmbedVecTo(eOut.Coeffs, e.Coeffs)
	eOut.IsNTT = false
}

// EmbedVec returns the embedding of v to the output modulus.
// If len(v) < len(emb.modIn), it only embeds the first len(v) elements.
func (emb *Embedder) EmbedVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(emb.modOut))
	for i := 0; i < len(emb.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	emb.EmbedVecTo(vOut, v)
	return vOut
}

// EmbedVecTo embeds v to vOut.
// If len(vOut) < len(emb.modOut),
// it only embeds to len(vOut) elements.
func (emb *Embedder) EmbedVecTo(vOut, v [][]uint64) {
	M := (len(v[0]) >> 3) << 3

	inLen, outLen := len(v), min(len(vOut), len(emb.modOut))

	if inLen != len(emb.modIn) {
		panic("input(s) not consistent")
	}

	if inLen == 1 {
		qv := emb.modIn[0].Value()
		halfQv := qv >> 1

		bufIn := (*[8]uint64)(unsafe.Pointer(&emb.buf.in[0][0]))

		for i := 0; i < M; i += 8 {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[0][i]))
			copy(bufIn[:], wIn[:])

			for j := 0; j < outLen; j++ {
				wOut := (*[8]uint64)(unsafe.Pointer(&vOut[j][i]))

				if emb.idx[j] == 0 {
					copy(wOut[:], bufIn[:])
				} else {
					wOut[0] = reduceModInToModOutSigned(bufIn[0], emb.modOut[j], qv, halfQv)
					wOut[1] = reduceModInToModOutSigned(bufIn[1], emb.modOut[j], qv, halfQv)
					wOut[2] = reduceModInToModOutSigned(bufIn[2], emb.modOut[j], qv, halfQv)
					wOut[3] = reduceModInToModOutSigned(bufIn[3], emb.modOut[j], qv, halfQv)

					wOut[4] = reduceModInToModOutSigned(bufIn[4], emb.modOut[j], qv, halfQv)
					wOut[5] = reduceModInToModOutSigned(bufIn[5], emb.modOut[j], qv, halfQv)
					wOut[6] = reduceModInToModOutSigned(bufIn[6], emb.modOut[j], qv, halfQv)
					wOut[7] = reduceModInToModOutSigned(bufIn[7], emb.modOut[j], qv, halfQv)
				}
			}
		}

		for i := M; i < len(v[0]); i++ {
			vi := v[0][i]
			for j := 0; j < outLen; j++ {
				if emb.idx[j] == 0 {
					vOut[j][i] = vi
				} else {
					vOut[j][i] = reduceModInToModOutSigned(vi, emb.modOut[j], qv, halfQv)
				}
			}
		}
		return
	}

	f128 := (*[8]float128.Float128)(unsafe.Pointer(&emb.buf.f128[0]))
	i64 := (*[8]uint64)(unsafe.Pointer(&emb.buf.u64[0]))

	for k := 0; k < M; k += 8 {
		clear(f128[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&emb.buf.in[i][0]))

			compInv, compInvS := emb.compInv[i], emb.compInvS[i]
			inv := emb.inv[i]
			modIn := emb.modIn[i]

			bufIn[0] = num.SMul(wIn[0], compInv, compInvS, modIn)
			bufIn[1] = num.SMul(wIn[1], compInv, compInvS, modIn)
			bufIn[2] = num.SMul(wIn[2], compInv, compInvS, modIn)
			bufIn[3] = num.SMul(wIn[3], compInv, compInvS, modIn)

			bufIn[4] = num.SMul(wIn[4], compInv, compInvS, modIn)
			bufIn[5] = num.SMul(wIn[5], compInv, compInvS, modIn)
			bufIn[6] = num.SMul(wIn[6], compInv, compInvS, modIn)
			bufIn[7] = num.SMul(wIn[7], compInv, compInvS, modIn)

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(bufIn[0])), inv))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(bufIn[1])), inv))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(bufIn[2])), inv))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(bufIn[3])), inv))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(bufIn[4])), inv))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(bufIn[5])), inv))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(bufIn[6])), inv))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(bufIn[7])), inv))
		}

		i64[0] = float128.ToUint64(f128[0])
		i64[1] = float128.ToUint64(f128[1])
		i64[2] = float128.ToUint64(f128[2])
		i64[3] = float128.ToUint64(f128[3])

		i64[4] = float128.ToUint64(f128[4])
		i64[5] = float128.ToUint64(f128[5])
		i64[6] = float128.ToUint64(f128[6])
		i64[7] = float128.ToUint64(f128[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				copy(wOut[:], v[emb.idx[i]][k:k+8])
			} else {
				negMod, negModS := emb.negMod[i], emb.negModS[i]
				modOut := emb.modOut[i]

				wOut[0] = num.SMul(i64[0], negMod, negModS, modOut)
				wOut[1] = num.SMul(i64[1], negMod, negModS, modOut)
				wOut[2] = num.SMul(i64[2], negMod, negModS, modOut)
				wOut[3] = num.SMul(i64[3], negMod, negModS, modOut)

				wOut[4] = num.SMul(i64[4], negMod, negModS, modOut)
				wOut[5] = num.SMul(i64[5], negMod, negModS, modOut)
				wOut[6] = num.SMul(i64[6], negMod, negModS, modOut)
				wOut[7] = num.SMul(i64[7], negMod, negModS, modOut)

				comp, compS := emb.comp[i], emb.compS[i]

				for j := 0; j < inLen; j++ {
					bufIn := (*[8]uint64)(unsafe.Pointer(&emb.buf.in[j][0]))

					comp, compS := comp[j], compS[j]

					wOut[0] = num.Add(wOut[0], num.SMul(bufIn[0], comp, compS, modOut), modOut)
					wOut[1] = num.Add(wOut[1], num.SMul(bufIn[1], comp, compS, modOut), modOut)
					wOut[2] = num.Add(wOut[2], num.SMul(bufIn[2], comp, compS, modOut), modOut)
					wOut[3] = num.Add(wOut[3], num.SMul(bufIn[3], comp, compS, modOut), modOut)

					wOut[4] = num.Add(wOut[4], num.SMul(bufIn[4], comp, compS, modOut), modOut)
					wOut[5] = num.Add(wOut[5], num.SMul(bufIn[5], comp, compS, modOut), modOut)
					wOut[6] = num.Add(wOut[6], num.SMul(bufIn[6], comp, compS, modOut), modOut)
					wOut[7] = num.Add(wOut[7], num.SMul(bufIn[7], comp, compS, modOut), modOut)
				}
			}
		}
	}

	for k := M; k < len(v[0]); k++ {
		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			emb.buf.in[i][0] = num.SMul(v[i][k], emb.compInv[i], emb.compInvS[i], emb.modIn[i])
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(emb.buf.in[i][0])), emb.inv[i]))
		}
		i64[0] = float128.ToUint64(f128[0])

		for i := 0; i < outLen; i++ {
			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				vOut[i][k] = v[emb.idx[i]][k]
			} else {
				vOut[i][k] = num.SMul(i64[0], emb.negMod[i], emb.negModS[i], emb.modOut[i])
				for j := 0; j < inLen; j++ {
					vOut[i][k] = num.Add(vOut[i][k], num.SMul(emb.buf.in[j][0], emb.comp[i][j], emb.compS[i][j], emb.modOut[i]), emb.modOut[i])
				}
			}
		}
	}
}

// ModulusIn returns the input modulus.
func (emb *Embedder) ModulusIn() []*num.Modulus {
	return emb.modIn
}

// ModulusOut returns the output modulus.
func (emb *Embedder) ModulusOut() []*num.Modulus {
	return emb.modOut
}

// SafeCopy returns a thread-safe copy.
func (emb *Embedder) SafeCopy() *Embedder {
	return &Embedder{
		modIn:  emb.modIn,
		modOut: emb.modOut,

		compInv:  emb.compInv,
		compInvS: emb.compInvS,

		comp:  emb.comp,
		compS: emb.compS,

		negMod:  emb.negMod,
		negModS: emb.negModS,

		inv: emb.inv,

		idx: emb.idx,

		buf: newEmbedderBuffer(emb.modIn),
	}
}

// Scaler scales a polynomial to different modulus.
// In other words, it computes
//
//	[p]_modIn -> [(modOut / modIn) * p]_modOut
//
// It uses HPS-like algorithm, so the computation is exact.
type Scaler struct {
	// modIn is the input modulus.
	modIn []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus

	// compInv is the modular inverse of the compliment of the input modulus limb.
	compInv []uint64
	// compInvS is the Shoup form of compInv.
	compInvS []uint64

	// scInt is the integer part of the modulus scaling factor.
	scInt [][]uint64
	// scIntS is the Shoup form of modScInt.
	scIntS [][]uint64

	// scFrac is the fractional part of the modulus scaling factor.
	scFrac []float128.Float128

	buf scalerBuffer
}

// scalerBuffer is a buffer for [Scaler].
type scalerBuffer struct {
	// f128 is a buffer for the 128-bit floating-point number.
	// Always has length 8.
	f128     []float128.Float128
	fLo, fHi []float64
	// uHi is a buffer for the high bits of the 64-bit integers.
	// Always has length 8.
	uHi []uint64
	// uLo is a buffer for the low bits of the 64-bit integers.
	// Always has length 8.
	uLo []uint64
	// in is a buffer for input coefficient.
	// Always has length [len(modIn)][8].
	in [][]uint64
}

// NewScaler creates a new [Scaler].
func NewScaler(modOut []*num.Modulus, modIn []*num.Modulus) *Scaler {
	if !isCoprime(modIn) || !isCoprime(modOut) {
		panic("modulus must be coprime")
	}

	compInv := make([]uint64, len(modIn))
	compInvS := make([]uint64, len(modIn))

	for i := 0; i < len(modIn); i++ {
		compInv[i] = 1
		for j := 0; j < len(modIn); j++ {
			if i != j {
				compInv[i] = num.Mul(compInv[i], num.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
		}
		compInvS[i] = num.SForm(compInv[i], modIn[i])
	}

	scInt := make([][]uint64, len(modOut))
	scIntS := make([][]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		scInt[i] = make([]uint64, len(modIn))
		scIntS[i] = make([]uint64, len(modIn))
	}

	scFrac := make([]float128.Float128, len(modIn))

	modOutBig := big.NewInt(1)
	tmpInt := big.NewInt(0)
	for i := 0; i < len(modOut); i++ {
		modOutBig.Mul(modOutBig, tmpInt.SetUint64(modOut[i].Value()))
	}

	scRat := big.NewRat(0, 1)
	scIntBig := big.NewInt(0)
	for i := 0; i < len(modIn); i++ {
		scIntBig.Div(modOutBig, tmpInt.SetUint64(modIn[i].Value()))
		for j := 0; j < len(modOut); j++ {
			tmpInt.Mod(scIntBig, tmpInt.SetUint64(modOut[j].Value()))
			scInt[j][i] = tmpInt.Uint64()
			scIntS[j][i] = num.SForm(scInt[j][i], modOut[j])
		}

		tmpInt.Mul(scIntBig, tmpInt.SetUint64(modIn[i].Value()))
		tmpInt.Sub(modOutBig, tmpInt)

		scRat.Denom().SetUint64(modIn[i].Value())
		scRat.Num().Set(tmpInt)
		scFrac[i] = float128.FromRat(scRat)
	}

	return &Scaler{
		modIn:  modIn,
		modOut: modOut,

		compInv:  compInv,
		compInvS: compInvS,

		scInt:  scInt,
		scIntS: scIntS,

		scFrac: scFrac,

		buf: newScalerBuffer(modIn),
	}
}

// newScalerBuffer creates a new [scalerBuffer].
func newScalerBuffer(modIn []*num.Modulus) scalerBuffer {
	in := make([][]uint64, len(modIn))
	for i := range in {
		in[i] = make([]uint64, 8)
	}

	return scalerBuffer{
		f128: make([]float128.Float128, 8),
		fLo:  make([]float64, 8),
		fHi:  make([]float64, 8),
		uHi:  make([]uint64, 8),
		uLo:  make([]uint64, 8),
		in:   in,
	}
}

// Scale returns the scaled element of e.
func (sc *Scaler) Scale(e *Element) *Element {
	pOut := NewPoly(e.Rank(), e.ModLen())
	sc.ScaleTo(pOut, e)
	return pOut
}

// ScaleTo scales e to eOut.
func (sc *Scaler) ScaleTo(eOut, e *Element) {
	checkBinaryOperable(e.Rank(), e.ModLen(), eOut, e)
	if e.Type() == TypePoly && e.IsNTT {
		panic("input(s) must be in standard form")
	}

	sc.ScaleVecTo(eOut.Coeffs, e.Coeffs)
	eOut.IsNTT = false
}

// ScaleVec returns the scaled vector of v.
func (sc *Scaler) ScaleVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(sc.modOut))
	for i := 0; i < len(sc.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	sc.ScaleVecTo(vOut, v)
	return vOut
}

// ScaleVecTo scales v to vOut.
func (sc *Scaler) ScaleVecTo(vOut, v [][]uint64) {
	if len(v) != len(sc.modIn) || len(vOut) != len(sc.modOut) {
		panic("input(s) not consistent")
	}

	inLen, outLen := len(sc.modIn), len(sc.modOut)
	M := (len(v[0]) >> 3) << 3

	f128 := (*[8]float128.Float128)(unsafe.Pointer(&sc.buf.f128[0]))
	uHi := (*[8]uint64)(unsafe.Pointer(&sc.buf.uHi[0]))
	uLo := (*[8]uint64)(unsafe.Pointer(&sc.buf.uLo[0]))

	for k := 0; k < M; k += 8 {
		clear(sc.buf.f128)
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&sc.buf.in[i][0]))

			compInv, compInvS := sc.compInv[i], sc.compInvS[i]
			scFrac := sc.scFrac[i]
			modIn := sc.modIn[i]

			bufIn[0] = num.SMul(wIn[0], compInv, compInvS, modIn)
			bufIn[1] = num.SMul(wIn[1], compInv, compInvS, modIn)
			bufIn[2] = num.SMul(wIn[2], compInv, compInvS, modIn)
			bufIn[3] = num.SMul(wIn[3], compInv, compInvS, modIn)

			bufIn[4] = num.SMul(wIn[4], compInv, compInvS, modIn)
			bufIn[5] = num.SMul(wIn[5], compInv, compInvS, modIn)
			bufIn[6] = num.SMul(wIn[6], compInv, compInvS, modIn)
			bufIn[7] = num.SMul(wIn[7], compInv, compInvS, modIn)

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(bufIn[0])), scFrac))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(bufIn[1])), scFrac))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(bufIn[2])), scFrac))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(bufIn[3])), scFrac))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(bufIn[4])), scFrac))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(bufIn[5])), scFrac))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(bufIn[6])), scFrac))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(bufIn[7])), scFrac))
		}

		uHi[0], uLo[0] = float128.ToUint128(f128[0])
		uHi[1], uLo[1] = float128.ToUint128(f128[1])
		uHi[2], uLo[2] = float128.ToUint128(f128[2])
		uHi[3], uLo[3] = float128.ToUint128(f128[3])

		uHi[4], uLo[4] = float128.ToUint128(f128[4])
		uHi[5], uLo[5] = float128.ToUint128(f128[5])
		uHi[6], uLo[6] = float128.ToUint128(f128[6])
		uHi[7], uLo[7] = float128.ToUint128(f128[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			modOut := sc.modOut[i]

			wOut[0] = num.Reduce128(uHi[0], uLo[0], modOut)
			wOut[1] = num.Reduce128(uHi[1], uLo[1], modOut)
			wOut[2] = num.Reduce128(uHi[2], uLo[2], modOut)
			wOut[3] = num.Reduce128(uHi[3], uLo[3], modOut)

			wOut[4] = num.Reduce128(uHi[4], uLo[4], modOut)
			wOut[5] = num.Reduce128(uHi[5], uLo[5], modOut)
			wOut[6] = num.Reduce128(uHi[6], uLo[6], modOut)
			wOut[7] = num.Reduce128(uHi[7], uLo[7], modOut)

			modScInt, modScIntS := sc.scInt[i], sc.scIntS[i]

			for j := 0; j < inLen; j++ {
				bufIn := (*[8]uint64)(unsafe.Pointer(&sc.buf.in[j][0]))

				modScInt, modScIntS := modScInt[j], modScIntS[j]

				wOut[0] = num.Add(wOut[0], num.SMul(bufIn[0], modScInt, modScIntS, modOut), modOut)
				wOut[1] = num.Add(wOut[1], num.SMul(bufIn[1], modScInt, modScIntS, modOut), modOut)
				wOut[2] = num.Add(wOut[2], num.SMul(bufIn[2], modScInt, modScIntS, modOut), modOut)
				wOut[3] = num.Add(wOut[3], num.SMul(bufIn[3], modScInt, modScIntS, modOut), modOut)

				wOut[4] = num.Add(wOut[4], num.SMul(bufIn[4], modScInt, modScIntS, modOut), modOut)
				wOut[5] = num.Add(wOut[5], num.SMul(bufIn[5], modScInt, modScIntS, modOut), modOut)
				wOut[6] = num.Add(wOut[6], num.SMul(bufIn[6], modScInt, modScIntS, modOut), modOut)
				wOut[7] = num.Add(wOut[7], num.SMul(bufIn[7], modScInt, modScIntS, modOut), modOut)
			}
		}
	}

	for k := M; k < len(v[0]); k++ {
		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			sc.buf.in[i][0] = num.SMul(v[i][k], sc.compInv[i], sc.compInvS[i], sc.modIn[i])
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(sc.buf.in[i][0])), sc.scFrac[i]))
		}
		uHi[0], uLo[0] = float128.ToUint128(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.Reduce128(uHi[0], uLo[0], sc.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(sc.buf.in[j][0], sc.scInt[i][j], sc.scIntS[i][j], sc.modOut[i]), sc.modOut[i])
			}
		}
	}
}

// ModulusIn returns the input modulus.
func (sc *Scaler) ModulusIn() []*num.Modulus {
	return sc.modIn
}

// ModulusOut returns the output modulus.
func (sc *Scaler) ModulusOut() []*num.Modulus {
	return sc.modOut
}

// SafeCopy returns a thread-safe copy.
func (sc *Scaler) SafeCopy() *Scaler {
	return &Scaler{
		modIn:  sc.modIn,
		modOut: sc.modOut,

		compInv:  sc.compInv,
		compInvS: sc.compInvS,

		scInt:  sc.scInt,
		scIntS: sc.scIntS,

		scFrac: sc.scFrac,

		buf: newScalerBuffer(sc.modIn),
	}
}

// ScaleEmbedder scales a polynomial and embeds it into different modulus.
// In other words, it computes
//
//	[p]_modIn -> [scale * p]_modOut
//
// for some rational scale.
// It uses HPS-like algorithm, so the computation is exact.
type ScaleEmbedder struct {
	// modIn is the input modulus of the ScaleEmbedder.
	modIn []*num.Modulus
	// modOut is the output modulus of the ScaleEmbedder.
	modOut []*num.Modulus
	// scale is the scaling factor.
	scale *big.Rat

	// compInv is the modular inverse of the compliment of the input modulus limb.
	compInv []uint64
	// compInvS is the Shoup form of compInv.
	compInvS []uint64

	// inv is the floating-point approximation of the inverse of the input modulus limb.
	inv []float128.Float128

	// scInt is the integer part of the modulus scaling factor.
	scInt [][]uint64
	// scIntS is the Shoup form of modScInt.
	scIntS [][]uint64

	// scFrac is the floating-point approximation of the modulus scaling factor.
	scFrac []float128.Float128

	// ovfInt is the integer part for the constant to compute the overflow multiplied by the scaling factor.
	ovfInt []uint64
	// ovfIntS is the Shoup form of ovfInt.
	ovfIntS []uint64

	// ovfFrac is the floating-point approximation of the fractional part
	// for the constant to compute the overflow multiplied by the scaling factor.
	ovfFrac float128.Float128

	buf scalerBuffer
}

// NewScaleEmbedder creates a new [ScaleEmbedder].
func NewScaleEmbedder(modOut []*num.Modulus, modIn []*num.Modulus, scale *big.Rat) *ScaleEmbedder {
	if !isCoprime(modIn) || !isCoprime(modOut) {
		panic("modulus must be coprime")
	}

	inLen, outLen := len(modIn), len(modOut)

	compInv := make([]uint64, inLen)
	compInvS := make([]uint64, inLen)

	for i := 0; i < inLen; i++ {
		compInv[i] = 1
		for j := 0; j < inLen; j++ {
			if i != j {
				compInv[i] = num.Mul(compInv[i], num.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
		}
		compInvS[i] = num.SForm(compInv[i], modIn[i])
	}

	inv := make([]float128.Float128, inLen)

	modInBig := big.NewInt(1)
	tmpInt := big.NewInt(0)
	tmpRat := big.NewRat(0, 1)
	modScIntBig := big.NewInt(0)

	for i := 0; i < inLen; i++ {
		modInBig.Mul(modInBig, tmpInt.SetUint64(modIn[i].Value()))

		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Num().SetInt64(1)
		inv[i] = float128.FromRat(tmpRat)
	}

	scInt := make([][]uint64, outLen)
	scIntS := make([][]uint64, outLen)
	for i := 0; i < outLen; i++ {
		scInt[i] = make([]uint64, inLen)
		scIntS[i] = make([]uint64, inLen)
	}

	scFrac := make([]float128.Float128, inLen)

	for i := 0; i < inLen; i++ {
		tmpRat.Num().Set(modInBig)
		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Mul(tmpRat, scale)

		modScIntBig.Div(tmpRat.Num(), tmpRat.Denom())
		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(modScIntBig, tmpRat.Denom()))

		for j := 0; j < outLen; j++ {
			tmpInt.Mod(modScIntBig, tmpInt.SetUint64(modOut[j].Value()))
			scInt[j][i] = tmpInt.Uint64()
			scIntS[j][i] = num.SForm(scInt[j][i], modOut[j])
		}
		scFrac[i] = float128.FromRat(tmpRat)
	}

	ovfIntBig := big.NewInt(0)

	ovfInt := make([]uint64, outLen)
	ovfIntS := make([]uint64, outLen)

	tmpRat.Set(scale)
	tmpRat.Num().Mul(tmpRat.Num(), modInBig)
	ovfIntBig.Div(tmpRat.Num(), tmpRat.Denom())
	for i := 0; i < outLen; i++ {
		tmpInt.Mod(ovfIntBig, tmpInt.SetUint64(modOut[i].Value()))
		ovfInt[i] = tmpInt.Uint64()
		ovfIntS[i] = num.SForm(ovfInt[i], modOut[i])
	}

	tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(ovfIntBig, tmpRat.Denom()))
	ovfFrac := float128.FromRat(tmpRat)

	return &ScaleEmbedder{
		modIn:  modIn,
		modOut: modOut,
		scale:  scale,

		compInv:  compInv,
		compInvS: compInvS,

		inv: inv,

		scInt:  scInt,
		scIntS: scIntS,

		scFrac: scFrac,

		ovfInt:  ovfInt,
		ovfIntS: ovfIntS,

		ovfFrac: ovfFrac,

		buf: newScalerBuffer(modIn),
	}
}

// ScaleEmbed scales and embeds e to the output modulus.
func (sc *ScaleEmbedder) ScaleEmbed(e *Element) *Element {
	pOut := NewPoly(e.Rank(), e.ModLen())
	sc.ScaleEmbedTo(pOut, e)
	return pOut
}

// ScaleEmbedTo scales and embeds e to eOut.
func (sc *ScaleEmbedder) ScaleEmbedTo(eOut, e *Element) {
	checkBinaryOperable(e.Rank(), e.ModLen(), eOut, e)
	if e.Type() != TypePoly && e.IsNTT {
		panic("input(s) must be in standard form")
	}

	sc.ScaleEmbedVecTo(eOut.Coeffs, e.Coeffs)
	eOut.IsNTT = false
}

// ScaleEmbedVec scales and embeds v to the output modulus.
func (sc *ScaleEmbedder) ScaleEmbedVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(sc.modOut))
	for i := 0; i < len(sc.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	sc.ScaleEmbedVecTo(vOut, v)
	return vOut
}

// ScaleEmbedVecTo scales and embeds v to vOut.
func (sc *ScaleEmbedder) ScaleEmbedVecTo(vOut, v [][]uint64) {
	if len(v) != len(sc.modIn) || len(vOut) != len(sc.modOut) {
		panic("input(s) not consistent")
	}

	inLen, outLen := len(sc.modIn), len(sc.modOut)
	M := (len(v[0]) >> 3) << 3

	f128 := (*[8]float128.Float128)(unsafe.Pointer(&sc.buf.f128[0]))
	uHi := (*[8]uint64)(unsafe.Pointer(&sc.buf.uHi[0]))
	uLo := (*[8]uint64)(unsafe.Pointer(&sc.buf.uLo[0]))

	for k := 0; k < M; k += 8 {
		clear(f128[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&sc.buf.in[i][0]))

			compInv, compInvS, modIn := sc.compInv[i], sc.compInvS[i], sc.modIn[i]
			inv := sc.inv[i]

			bufIn[0] = num.SMul(wIn[0], compInv, compInvS, modIn)
			bufIn[1] = num.SMul(wIn[1], compInv, compInvS, modIn)
			bufIn[2] = num.SMul(wIn[2], compInv, compInvS, modIn)
			bufIn[3] = num.SMul(wIn[3], compInv, compInvS, modIn)

			bufIn[4] = num.SMul(wIn[4], compInv, compInvS, modIn)
			bufIn[5] = num.SMul(wIn[5], compInv, compInvS, modIn)
			bufIn[6] = num.SMul(wIn[6], compInv, compInvS, modIn)
			bufIn[7] = num.SMul(wIn[7], compInv, compInvS, modIn)

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(bufIn[0])), inv))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(bufIn[1])), inv))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(bufIn[2])), inv))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(bufIn[3])), inv))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(bufIn[4])), inv))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(bufIn[5])), inv))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(bufIn[6])), inv))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(bufIn[7])), inv))
		}

		uLo[0] = float128.ToUint64(f128[0])
		uLo[1] = float128.ToUint64(f128[1])
		uLo[2] = float128.ToUint64(f128[2])
		uLo[3] = float128.ToUint64(f128[3])

		uLo[4] = float128.ToUint64(f128[4])
		uLo[5] = float128.ToUint64(f128[5])
		uLo[6] = float128.ToUint64(f128[6])
		uLo[7] = float128.ToUint64(f128[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			intOv, intOvS := sc.ovfInt[i], sc.ovfIntS[i]
			modOut, modOutv := sc.modOut[i], sc.modOut[i].Value()

			wOut[0] = num.SMul(modOutv-uLo[0], intOv, intOvS, modOut)
			wOut[1] = num.SMul(modOutv-uLo[1], intOv, intOvS, modOut)
			wOut[2] = num.SMul(modOutv-uLo[2], intOv, intOvS, modOut)
			wOut[3] = num.SMul(modOutv-uLo[3], intOv, intOvS, modOut)

			wOut[4] = num.SMul(modOutv-uLo[4], intOv, intOvS, modOut)
			wOut[5] = num.SMul(modOutv-uLo[5], intOv, intOvS, modOut)
			wOut[6] = num.SMul(modOutv-uLo[6], intOv, intOvS, modOut)
			wOut[7] = num.SMul(modOutv-uLo[7], intOv, intOvS, modOut)

			modScInt, modScIntS := sc.scInt[i], sc.scIntS[i]

			for j := 0; j < inLen; j++ {
				bufIn := (*[8]uint64)(unsafe.Pointer(&sc.buf.in[j][0]))

				modScInt, modScIntS := modScInt[j], modScIntS[j]

				wOut[0] = num.Add(wOut[0], num.SMul(bufIn[0], modScInt, modScIntS, modOut), modOut)
				wOut[1] = num.Add(wOut[1], num.SMul(bufIn[1], modScInt, modScIntS, modOut), modOut)
				wOut[2] = num.Add(wOut[2], num.SMul(bufIn[2], modScInt, modScIntS, modOut), modOut)
				wOut[3] = num.Add(wOut[3], num.SMul(bufIn[3], modScInt, modScIntS, modOut), modOut)

				wOut[4] = num.Add(wOut[4], num.SMul(bufIn[4], modScInt, modScIntS, modOut), modOut)
				wOut[5] = num.Add(wOut[5], num.SMul(bufIn[5], modScInt, modScIntS, modOut), modOut)
				wOut[6] = num.Add(wOut[6], num.SMul(bufIn[6], modScInt, modScIntS, modOut), modOut)
				wOut[7] = num.Add(wOut[7], num.SMul(bufIn[7], modScInt, modScIntS, modOut), modOut)
			}
		}

		clear(f128[:])
		for i := 0; i < inLen; i++ {
			bufIn := (*[8]uint64)(unsafe.Pointer(&sc.buf.in[i][0]))

			modScFrac := sc.scFrac[i]

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(bufIn[0])), modScFrac))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(bufIn[1])), modScFrac))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(bufIn[2])), modScFrac))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(bufIn[3])), modScFrac))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(bufIn[4])), modScFrac))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(bufIn[5])), modScFrac))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(bufIn[6])), modScFrac))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(bufIn[7])), modScFrac))
		}

		f128[0] = float128.Sub(f128[0], float128.Mul(float128.FromInt64(int64(uLo[0])), sc.ovfFrac))
		f128[1] = float128.Sub(f128[1], float128.Mul(float128.FromInt64(int64(uLo[1])), sc.ovfFrac))
		f128[2] = float128.Sub(f128[2], float128.Mul(float128.FromInt64(int64(uLo[2])), sc.ovfFrac))
		f128[3] = float128.Sub(f128[3], float128.Mul(float128.FromInt64(int64(uLo[3])), sc.ovfFrac))

		f128[4] = float128.Sub(f128[4], float128.Mul(float128.FromInt64(int64(uLo[4])), sc.ovfFrac))
		f128[5] = float128.Sub(f128[5], float128.Mul(float128.FromInt64(int64(uLo[5])), sc.ovfFrac))
		f128[6] = float128.Sub(f128[6], float128.Mul(float128.FromInt64(int64(uLo[6])), sc.ovfFrac))
		f128[7] = float128.Sub(f128[7], float128.Mul(float128.FromInt64(int64(uLo[7])), sc.ovfFrac))

		uHi[0], uLo[0] = float128.ToUint128(f128[0])
		uHi[1], uLo[1] = float128.ToUint128(f128[1])
		uHi[2], uLo[2] = float128.ToUint128(f128[2])
		uHi[3], uLo[3] = float128.ToUint128(f128[3])

		uHi[4], uLo[4] = float128.ToUint128(f128[4])
		uHi[5], uLo[5] = float128.ToUint128(f128[5])
		uHi[6], uLo[6] = float128.ToUint128(f128[6])
		uHi[7], uLo[7] = float128.ToUint128(f128[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			modOut := sc.modOut[i]

			wOut[0] = add64To128Signed(wOut[0], uHi[0], uLo[0], modOut)
			wOut[1] = add64To128Signed(wOut[1], uHi[1], uLo[1], modOut)
			wOut[2] = add64To128Signed(wOut[2], uHi[2], uLo[2], modOut)
			wOut[3] = add64To128Signed(wOut[3], uHi[3], uLo[3], modOut)

			wOut[4] = add64To128Signed(wOut[4], uHi[4], uLo[4], modOut)
			wOut[5] = add64To128Signed(wOut[5], uHi[5], uLo[5], modOut)
			wOut[6] = add64To128Signed(wOut[6], uHi[6], uLo[6], modOut)
			wOut[7] = add64To128Signed(wOut[7], uHi[7], uLo[7], modOut)
		}
	}

	for k := M; k < len(v[0]); k++ {
		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			sc.buf.in[i][0] = num.SMul(v[i][k], sc.compInv[i], sc.compInvS[i], sc.modIn[i])

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(sc.buf.in[i][0])), sc.inv[i]))
		}

		uLo[0] = float128.ToUint64(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.SMul(sc.modOut[i].Value()-uLo[0], sc.ovfInt[i], sc.ovfIntS[i], sc.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(sc.buf.in[j][0], sc.scInt[i][j], sc.scIntS[i][j], sc.modOut[i]), sc.modOut[i])
			}
		}

		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(sc.buf.in[i][0])), sc.scFrac[i]))
		}

		f128[0] = float128.Sub(f128[0], float128.Mul(float128.FromInt64(int64(uLo[0])), sc.ovfFrac))

		uHi[0], uLo[0] = float128.ToUint128(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = add64To128Signed(vOut[i][k], uHi[0], uLo[0], sc.modOut[i])
		}
	}
}

// ModulusIn returns the input modulus.
func (sc *ScaleEmbedder) ModulusIn() []*num.Modulus {
	return sc.modIn
}

// ModulusOut returns the output modulus.
func (sc *ScaleEmbedder) ModulusOut() []*num.Modulus {
	return sc.modOut
}

// ScalingFactor returns the scaling factor.
func (sc *ScaleEmbedder) ScalingFactor() *big.Rat {
	return sc.scale
}

// SafeCopy returns a thread-safe copy.
func (sc *ScaleEmbedder) SafeCopy() *ScaleEmbedder {
	return &ScaleEmbedder{
		modIn:  sc.modIn,
		modOut: sc.modOut,

		compInv:  sc.compInv,
		compInvS: sc.compInvS,

		inv: sc.inv,

		scInt:  sc.scInt,
		scIntS: sc.scIntS,

		scFrac: sc.scFrac,

		ovfInt:  sc.ovfInt,
		ovfIntS: sc.ovfIntS,

		ovfFrac: sc.ovfFrac,

		buf: newScalerBuffer(sc.modIn),
	}
}
