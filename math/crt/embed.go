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

// Embed returns the embedding of p to the output modulus.
// If p.ModLen() < len(e.modIn), it only embeds the first p.ModLen() elements.
func (e *Embedder) Embed(p *Element) *Element {
	pOut := NewPoly(p.Rank(), p.ModLen())
	e.EmbedTo(pOut, p)
	return pOut
}

// EmbedTo embeds the p to pOut.
// If p.ModLen() < len(e.modIn) or pOut.ModLen() < len(e.modOut),
// it only embeds the first p.ModLen() elements to pOut.ModLen() elements.
func (e *Embedder) EmbedTo(pOut, p *Element) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	}

	e.EmbedVecTo(pOut.Coeffs, p.Coeffs)
	pOut.IsNTT = false
}

// EmbedVec returns the embedding of v to the output modulus.
// If len(v) < len(e.modIn), it only embeds the first len(v) elements.
func (e *Embedder) EmbedVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(e.modOut))
	for i := 0; i < len(e.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	e.EmbedVecTo(vOut, v)
	return vOut
}

// EmbedVecTo embeds v to vOut.
// If len(vOut) < len(e.modOut),
// it only embeds to len(vOut) elements.
func (e *Embedder) EmbedVecTo(vOut, v [][]uint64) {
	M := (len(v[0]) >> 3) << 3

	inLen, outLen := len(v), min(len(vOut), len(e.modOut))

	if inLen != len(e.modIn) {
		panic("input(s) not consistent")
	}

	if inLen == 1 {
		qv := e.modIn[0].Value()
		halfQv := qv >> 1

		bufIn := (*[8]uint64)(unsafe.Pointer(&e.buf.in[0][0]))

		for i := 0; i < M; i += 8 {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[0][i]))
			copy(bufIn[:], wIn[:])

			for j := 0; j < outLen; j++ {
				wOut := (*[8]uint64)(unsafe.Pointer(&vOut[j][i]))

				if e.idx[j] == 0 {
					copy(wOut[:], bufIn[:])
				} else {
					wOut[0] = reduceModInToModOutSigned(bufIn[0], e.modOut[j], qv, halfQv)
					wOut[1] = reduceModInToModOutSigned(bufIn[1], e.modOut[j], qv, halfQv)
					wOut[2] = reduceModInToModOutSigned(bufIn[2], e.modOut[j], qv, halfQv)
					wOut[3] = reduceModInToModOutSigned(bufIn[3], e.modOut[j], qv, halfQv)

					wOut[4] = reduceModInToModOutSigned(bufIn[4], e.modOut[j], qv, halfQv)
					wOut[5] = reduceModInToModOutSigned(bufIn[5], e.modOut[j], qv, halfQv)
					wOut[6] = reduceModInToModOutSigned(bufIn[6], e.modOut[j], qv, halfQv)
					wOut[7] = reduceModInToModOutSigned(bufIn[7], e.modOut[j], qv, halfQv)
				}
			}
		}

		for i := M; i < len(v[0]); i++ {
			vi := v[0][i]
			for j := 0; j < outLen; j++ {
				if e.idx[j] == 0 {
					vOut[j][i] = vi
				} else {
					vOut[j][i] = reduceModInToModOutSigned(vi, e.modOut[j], qv, halfQv)
				}
			}
		}
		return
	}

	f128 := (*[8]float128.Float128)(unsafe.Pointer(&e.buf.f128[0]))
	i64 := (*[8]uint64)(unsafe.Pointer(&e.buf.u64[0]))

	for k := 0; k < M; k += 8 {
		clear(f128[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&e.buf.in[i][0]))

			compInv, compInvS := e.compInv[i], e.compInvS[i]
			inv := e.inv[i]
			modIn := e.modIn[i]

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

			if 0 <= e.idx[i] && e.idx[i] < inLen {
				copy(wOut[:], v[e.idx[i]][k:k+8])
			} else {
				negMod, negModS := e.negMod[i], e.negModS[i]
				modOut := e.modOut[i]

				wOut[0] = num.SMul(i64[0], negMod, negModS, modOut)
				wOut[1] = num.SMul(i64[1], negMod, negModS, modOut)
				wOut[2] = num.SMul(i64[2], negMod, negModS, modOut)
				wOut[3] = num.SMul(i64[3], negMod, negModS, modOut)

				wOut[4] = num.SMul(i64[4], negMod, negModS, modOut)
				wOut[5] = num.SMul(i64[5], negMod, negModS, modOut)
				wOut[6] = num.SMul(i64[6], negMod, negModS, modOut)
				wOut[7] = num.SMul(i64[7], negMod, negModS, modOut)

				comp, compS := e.comp[i], e.compS[i]

				for j := 0; j < inLen; j++ {
					bufIn := (*[8]uint64)(unsafe.Pointer(&e.buf.in[j][0]))

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
			e.buf.in[i][0] = num.SMul(v[i][k], e.compInv[i], e.compInvS[i], e.modIn[i])
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(e.buf.in[i][0])), e.inv[i]))
		}
		i64[0] = float128.ToUint64(f128[0])

		for i := 0; i < outLen; i++ {
			if 0 <= e.idx[i] && e.idx[i] < inLen {
				vOut[i][k] = v[e.idx[i]][k]
			} else {
				vOut[i][k] = num.SMul(i64[0], e.negMod[i], e.negModS[i], e.modOut[i])
				for j := 0; j < inLen; j++ {
					vOut[i][k] = num.Add(vOut[i][k], num.SMul(e.buf.in[j][0], e.comp[i][j], e.compS[i][j], e.modOut[i]), e.modOut[i])
				}
			}
		}
	}
}

// ModulusIn returns the input modulus.
func (e *Embedder) ModulusIn() []*num.Modulus {
	return e.modIn
}

// ModulusOut returns the output modulus.
func (e *Embedder) ModulusOut() []*num.Modulus {
	return e.modOut
}

// SafeCopy returns a thread-safe copy.
func (e *Embedder) SafeCopy() *Embedder {
	return &Embedder{
		modIn:  e.modIn,
		modOut: e.modOut,

		compInv:  e.compInv,
		compInvS: e.compInvS,

		comp:  e.comp,
		compS: e.compS,

		negMod:  e.negMod,
		negModS: e.negModS,

		inv: e.inv,

		idx: e.idx,

		buf: newEmbedderBuffer(e.modIn),
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

// Scale returns the scaled polynomial of p.
func (s *Scaler) Scale(p *Element) *Element {
	pOut := NewPoly(p.Rank(), p.ModLen())
	s.ScaleTo(pOut, p)
	return pOut
}

// ScaleTo scales p to pOut.
func (s *Scaler) ScaleTo(pOut, p *Element) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	}

	s.ScaleVecTo(pOut.Coeffs, p.Coeffs)
	pOut.IsNTT = false
}

// ScaleVec returns the scaled vector of v.
func (s *Scaler) ScaleVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(s.modOut))
	for i := 0; i < len(s.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	s.ScaleVecTo(vOut, v)
	return vOut
}

// ScaleVecTo scales v to vOut.
func (s *Scaler) ScaleVecTo(vOut, v [][]uint64) {
	if len(v) != len(s.modIn) || len(vOut) != len(s.modOut) {
		panic("input(s) not consistent")
	}

	inLen, outLen := len(s.modIn), len(s.modOut)
	M := (len(v[0]) >> 3) << 3

	f128 := (*[8]float128.Float128)(unsafe.Pointer(&s.buf.f128[0]))
	uHi := (*[8]uint64)(unsafe.Pointer(&s.buf.uHi[0]))
	uLo := (*[8]uint64)(unsafe.Pointer(&s.buf.uLo[0]))

	for k := 0; k < M; k += 8 {
		clear(s.buf.f128)
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[i][0]))

			compInv, compInvS := s.compInv[i], s.compInvS[i]
			scFrac := s.scFrac[i]
			modIn := s.modIn[i]

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

			modOut := s.modOut[i]

			wOut[0] = num.Reduce128(uHi[0], uLo[0], modOut)
			wOut[1] = num.Reduce128(uHi[1], uLo[1], modOut)
			wOut[2] = num.Reduce128(uHi[2], uLo[2], modOut)
			wOut[3] = num.Reduce128(uHi[3], uLo[3], modOut)

			wOut[4] = num.Reduce128(uHi[4], uLo[4], modOut)
			wOut[5] = num.Reduce128(uHi[5], uLo[5], modOut)
			wOut[6] = num.Reduce128(uHi[6], uLo[6], modOut)
			wOut[7] = num.Reduce128(uHi[7], uLo[7], modOut)

			modScInt, modScIntS := s.scInt[i], s.scIntS[i]

			for j := 0; j < inLen; j++ {
				bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[j][0]))

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
			s.buf.in[i][0] = num.SMul(v[i][k], s.compInv[i], s.compInvS[i], s.modIn[i])
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(s.buf.in[i][0])), s.scFrac[i]))
		}
		uHi[0], uLo[0] = float128.ToUint128(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.Reduce128(uHi[0], uLo[0], s.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(s.buf.in[j][0], s.scInt[i][j], s.scIntS[i][j], s.modOut[i]), s.modOut[i])
			}
		}
	}
}

// ModulusIn returns the input modulus.
func (s *Scaler) ModulusIn() []*num.Modulus {
	return s.modIn
}

// ModulusOut returns the output modulus.
func (s *Scaler) ModulusOut() []*num.Modulus {
	return s.modOut
}

// SafeCopy returns a thread-safe copy.
func (s *Scaler) SafeCopy() *Scaler {
	return &Scaler{
		modIn:  s.modIn,
		modOut: s.modOut,

		compInv:  s.compInv,
		compInvS: s.compInvS,

		scInt:  s.scInt,
		scIntS: s.scIntS,

		scFrac: s.scFrac,

		buf: newScalerBuffer(s.modIn),
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

// ScaleEmbed scales and embeds p to the output modulus.
func (s *ScaleEmbedder) ScaleEmbed(p *Element) *Element {
	pOut := NewPoly(p.Rank(), p.ModLen())
	s.ScaleEmbedTo(pOut, p)
	return pOut
}

// ScaleEmbedTo scales and embeds p to pOut.
func (s *ScaleEmbedder) ScaleEmbedTo(pOut, p *Element) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	}

	s.ScaleEmbedVecTo(pOut.Coeffs, p.Coeffs)
	pOut.IsNTT = false
}

// ScaleEmbedVec scales and embeds v to the output modulus.
func (s *ScaleEmbedder) ScaleEmbedVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(s.modOut))
	for i := 0; i < len(s.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	s.ScaleEmbedVecTo(vOut, v)
	return vOut
}

// ScaleEmbedVecTo scales and embeds v to vOut.
func (s *ScaleEmbedder) ScaleEmbedVecTo(vOut, v [][]uint64) {
	if len(v) != len(s.modIn) || len(vOut) != len(s.modOut) {
		panic("input(s) not consistent")
	}

	inLen, outLen := len(s.modIn), len(s.modOut)
	M := (len(v[0]) >> 3) << 3

	f128 := (*[8]float128.Float128)(unsafe.Pointer(&s.buf.f128[0]))
	uHi := (*[8]uint64)(unsafe.Pointer(&s.buf.uHi[0]))
	uLo := (*[8]uint64)(unsafe.Pointer(&s.buf.uLo[0]))

	for k := 0; k < M; k += 8 {
		clear(f128[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[i][0]))

			compInv, compInvS, modIn := s.compInv[i], s.compInvS[i], s.modIn[i]
			inv := s.inv[i]

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

			intOv, intOvS := s.ovfInt[i], s.ovfIntS[i]
			modOut, modOutv := s.modOut[i], s.modOut[i].Value()

			wOut[0] = num.SMul(modOutv-uLo[0], intOv, intOvS, modOut)
			wOut[1] = num.SMul(modOutv-uLo[1], intOv, intOvS, modOut)
			wOut[2] = num.SMul(modOutv-uLo[2], intOv, intOvS, modOut)
			wOut[3] = num.SMul(modOutv-uLo[3], intOv, intOvS, modOut)

			wOut[4] = num.SMul(modOutv-uLo[4], intOv, intOvS, modOut)
			wOut[5] = num.SMul(modOutv-uLo[5], intOv, intOvS, modOut)
			wOut[6] = num.SMul(modOutv-uLo[6], intOv, intOvS, modOut)
			wOut[7] = num.SMul(modOutv-uLo[7], intOv, intOvS, modOut)

			modScInt, modScIntS := s.scInt[i], s.scIntS[i]

			for j := 0; j < inLen; j++ {
				bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[j][0]))

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
			bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[i][0]))

			modScFrac := s.scFrac[i]

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(bufIn[0])), modScFrac))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(bufIn[1])), modScFrac))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(bufIn[2])), modScFrac))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(bufIn[3])), modScFrac))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(bufIn[4])), modScFrac))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(bufIn[5])), modScFrac))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(bufIn[6])), modScFrac))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(bufIn[7])), modScFrac))
		}

		f128[0] = float128.Sub(f128[0], float128.Mul(float128.FromInt64(int64(uLo[0])), s.ovfFrac))
		f128[1] = float128.Sub(f128[1], float128.Mul(float128.FromInt64(int64(uLo[1])), s.ovfFrac))
		f128[2] = float128.Sub(f128[2], float128.Mul(float128.FromInt64(int64(uLo[2])), s.ovfFrac))
		f128[3] = float128.Sub(f128[3], float128.Mul(float128.FromInt64(int64(uLo[3])), s.ovfFrac))

		f128[4] = float128.Sub(f128[4], float128.Mul(float128.FromInt64(int64(uLo[4])), s.ovfFrac))
		f128[5] = float128.Sub(f128[5], float128.Mul(float128.FromInt64(int64(uLo[5])), s.ovfFrac))
		f128[6] = float128.Sub(f128[6], float128.Mul(float128.FromInt64(int64(uLo[6])), s.ovfFrac))
		f128[7] = float128.Sub(f128[7], float128.Mul(float128.FromInt64(int64(uLo[7])), s.ovfFrac))

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

			modOut := s.modOut[i]

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
			s.buf.in[i][0] = num.SMul(v[i][k], s.compInv[i], s.compInvS[i], s.modIn[i])

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(s.buf.in[i][0])), s.inv[i]))
		}

		uLo[0] = float128.ToUint64(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.SMul(s.modOut[i].Value()-uLo[0], s.ovfInt[i], s.ovfIntS[i], s.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(s.buf.in[j][0], s.scInt[i][j], s.scIntS[i][j], s.modOut[i]), s.modOut[i])
			}
		}

		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(s.buf.in[i][0])), s.scFrac[i]))
		}

		f128[0] = float128.Sub(f128[0], float128.Mul(float128.FromInt64(int64(uLo[0])), s.ovfFrac))

		uHi[0], uLo[0] = float128.ToUint128(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = add64To128Signed(vOut[i][k], uHi[0], uLo[0], s.modOut[i])
		}
	}
}

// ModulusIn returns the input modulus.
func (s *ScaleEmbedder) ModulusIn() []*num.Modulus {
	return s.modIn
}

// ModulusOut returns the output modulus.
func (s *ScaleEmbedder) ModulusOut() []*num.Modulus {
	return s.modOut
}

// ScalingFactor returns the scaling factor.
func (s *ScaleEmbedder) ScalingFactor() *big.Rat {
	return s.scale
}

// SafeCopy returns a thread-safe copy.
func (s *ScaleEmbedder) SafeCopy() *ScaleEmbedder {
	return &ScaleEmbedder{
		modIn:  s.modIn,
		modOut: s.modOut,

		compInv:  s.compInv,
		compInvS: s.compInvS,

		inv: s.inv,

		scInt:  s.scInt,
		scIntS: s.scIntS,

		scFrac: s.scFrac,

		ovfInt:  s.ovfInt,
		ovfIntS: s.ovfIntS,

		ovfFrac: s.ovfFrac,

		buf: newScalerBuffer(s.modIn),
	}
}
