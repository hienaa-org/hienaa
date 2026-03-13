package crt

import (
	"math/big"
	"sync"
	"unsafe"

	"github.com/hienaa-org/hienaa/math/internal/float128"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

const (
	// embedBatch is the batch size for embedding.
	embedBatch = 1 << logEmbedBatch
	// logEmbedBatch equals log2(embedBatch).
	logEmbedBatch = 8
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

	// idx holds the index of the input modulus limb if it overlaps with the output modulus limb.
	// For example, if modOut[i] = modIn[j], then idx[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	idx []int

	u64Pool *sync.Pool
}

// NewEmbedder creates a new [Embedder].
func NewEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *Embedder {
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

		idx: idx,

		u64Pool: &sync.Pool{
			New: func() any {
				return &[embedBatch]uint64{}
			},
		},
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
	if e.Type() == TypePoly && e.IsNTT {
		panic("input(s) must be in standard form")
	} else if eOut.Rank() != e.Rank() {
		panic("input(s) not consistent")
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
	M := (len(v[0]) >> logEmbedBatch) << logEmbedBatch
	L := unsafe.Sizeof(uint64(0))

	inLen, outLen := len(v), min(len(vOut), len(emb.modOut))

	if inLen != len(emb.modIn) {
		panic("input(s) not consistent")
	}

	if inLen == 1 {
		qv := emb.modIn[0].Value()
		halfQv := qv >> 1

		vBuf := emb.u64Pool.Get().(*[embedBatch]uint64)
		defer emb.u64Pool.Put(vBuf)

		r := unsafe.Pointer(unsafe.SliceData(v[0]))

		for i := 0; i < M; i += 8 {
			wIn := (*[8]uint64)(unsafe.Add(r, uintptr(i)*L))
			copy(vBuf[:], wIn[:])

			for j := 0; j < outLen; j++ {
				rOut := unsafe.Pointer(unsafe.SliceData(vOut[j]))
				wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))

				if emb.idx[j] == 0 {
					copy(wOut[:], vBuf[:])
				} else {
					wOut[0] = reduceModInToModOutSigned(vBuf[0], emb.modOut[j], qv, halfQv)
					wOut[1] = reduceModInToModOutSigned(vBuf[1], emb.modOut[j], qv, halfQv)
					wOut[2] = reduceModInToModOutSigned(vBuf[2], emb.modOut[j], qv, halfQv)
					wOut[3] = reduceModInToModOutSigned(vBuf[3], emb.modOut[j], qv, halfQv)

					wOut[4] = reduceModInToModOutSigned(vBuf[4], emb.modOut[j], qv, halfQv)
					wOut[5] = reduceModInToModOutSigned(vBuf[5], emb.modOut[j], qv, halfQv)
					wOut[6] = reduceModInToModOutSigned(vBuf[6], emb.modOut[j], qv, halfQv)
					wOut[7] = reduceModInToModOutSigned(vBuf[7], emb.modOut[j], qv, halfQv)
				}
			}
		}

		for i := M; i < len(v[0]); i++ {
			c := v[0][i]
			for j := 0; j < outLen; j++ {
				if emb.idx[j] == 0 {
					vOut[j][i] = c
				} else {
					vOut[j][i] = reduceModInToModOutSigned(c, emb.modOut[j], qv, halfQv)
				}
			}
		}

		return
	}

	fHi := emb.u64Pool.Get().(*[embedBatch]uint64)
	defer emb.u64Pool.Put(fHi)
	fLo := emb.u64Pool.Get().(*[embedBatch]uint64)
	defer emb.u64Pool.Put(fLo)

	vBuf := make([]*[embedBatch]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = emb.u64Pool.Get().(*[embedBatch]uint64)
		defer emb.u64Pool.Put(vBuf[i])
	}

	for k := 0; k < M; k += embedBatch {
		clear(fHi[:])
		clear(fLo[:])
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))
			wBuf := vBuf[i]

			compInv, compInvS := emb.compInv[i], emb.compInvS[i]
			modIn := emb.modIn[i]
			invHi, invLo := modIn.Div()

			vec.SMulScalarTo(wBuf[:], w[:], compInv, compInvS, modIn)

			for j := 0; j < embedBatch; j += 8 {
				fHi[j+0], fLo[j+0] = f128Acc(wBuf[j+0], invHi, invLo, fHi[j+0], fLo[j+0])
				fHi[j+1], fLo[j+1] = f128Acc(wBuf[j+1], invHi, invLo, fHi[j+1], fLo[j+1])
				fHi[j+2], fLo[j+2] = f128Acc(wBuf[j+2], invHi, invLo, fHi[j+2], fLo[j+2])
				fHi[j+3], fLo[j+3] = f128Acc(wBuf[j+3], invHi, invLo, fHi[j+3], fLo[j+3])

				fHi[j+4], fLo[j+4] = f128Acc(wBuf[j+4], invHi, invLo, fHi[j+4], fLo[j+4])
				fHi[j+5], fLo[j+5] = f128Acc(wBuf[j+5], invHi, invLo, fHi[j+5], fLo[j+5])
				fHi[j+6], fLo[j+6] = f128Acc(wBuf[j+6], invHi, invLo, fHi[j+6], fLo[j+6])
				fHi[j+7], fLo[j+7] = f128Acc(wBuf[j+7], invHi, invLo, fHi[j+7], fLo[j+7])
			}
		}

		for j := 0; j < embedBatch; j += 8 {
			fHi[j+0] += fLo[j+0] >> 63
			fHi[j+1] += fLo[j+1] >> 63
			fHi[j+2] += fLo[j+2] >> 63
			fHi[j+3] += fLo[j+3] >> 63

			fHi[j+4] += fLo[j+4] >> 63
			fHi[j+5] += fLo[j+5] >> 63
			fHi[j+6] += fLo[j+6] >> 63
			fHi[j+7] += fLo[j+7] >> 63
		}

		for i := 0; i < outLen; i++ {
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))

			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				copy(wOut[:], v[emb.idx[i]][k:k+embedBatch])
			} else {
				negMod, negModS := emb.negMod[i], emb.negModS[i]
				modOut := emb.modOut[i]

				vec.SMulScalarTo(wOut[:], fHi[:], negMod, negModS, modOut)

				comp, compS := emb.comp[i], emb.compS[i]
				for j := 0; j < inLen; j++ {
					wBuf := vBuf[j]
					comp, compS := comp[j], compS[j]

					vec.SMulAddScalarTo(wOut[:], wBuf[:], comp, compS, modOut)
				}
			}
		}
	}

	clear(fHi[:])
	clear(fLo[:])
	for i := 0; i < inLen; i++ {
		w := v[i][M:]
		wBuf := vBuf[i][:len(v[0])-M]

		compInv, compInvS := emb.compInv[i], emb.compInvS[i]
		modIn := emb.modIn[i]
		invHi, invLo := modIn.Div()

		vec.SMulScalarTo(wBuf[:], w[:], compInv, compInvS, modIn)

		for j := 0; j < len(v[0])-M; j++ {
			fHi[j], fLo[j] = f128Acc(wBuf[j], invHi, invLo, fHi[j], fLo[j])
		}
	}

	for j := 0; j < len(v[0])-M; j++ {
		fHi[j] += fLo[j] >> 63
	}

	for i := 0; i < outLen; i++ {
		wOut := vOut[i][M:]

		if 0 <= emb.idx[i] && emb.idx[i] < inLen {
			copy(wOut, v[emb.idx[i]][M:])
		} else {
			negMod, negModS := emb.negMod[i], emb.negModS[i]
			modOut := emb.modOut[i]

			vec.SMulScalarTo(wOut, fHi[:len(v[0])-M], negMod, negModS, modOut)

			comp, compS := emb.comp[i], emb.compS[i]
			for j := 0; j < inLen; j++ {
				wBuf := vBuf[j]
				comp, compS := comp[j], compS[j]

				vec.SMulAddScalarTo(wOut, wBuf[:len(v[0])-M], comp, compS, modOut)
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

	f128Pool *sync.Pool
	u64Pool  *sync.Pool
}

// NewScaler creates a new [Scaler].
func NewScaler(modOut []*num.Modulus, modIn []*num.Modulus) *Scaler {
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

		f128Pool: &sync.Pool{
			New: func() any {
				return &[8]float128.Float128{}
			},
		},
		u64Pool: &sync.Pool{
			New: func() any {
				return &[8]uint64{}
			},
		},
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
	if !(e.Rank() == eOut.Rank() && e.ModLen() == len(sc.modIn) && eOut.ModLen() == len(sc.modOut)) {
		panic("input(s) not consistent")
	} else if e.IsNTT {
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
	L := unsafe.Sizeof(uint64(0))

	f128 := sc.f128Pool.Get().(*[8]float128.Float128)
	defer sc.f128Pool.Put(f128)
	uHi := sc.u64Pool.Get().(*[8]uint64)
	defer sc.u64Pool.Put(uHi)
	uLo := sc.u64Pool.Get().(*[8]uint64)
	defer sc.u64Pool.Put(uLo)
	vBuf := make([]*[8]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = sc.u64Pool.Get().(*[8]uint64)
		defer sc.u64Pool.Put(vBuf[i])
	}

	for k := 0; k < M; k += 8 {
		clear(f128[:])
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[8]uint64)(unsafe.Add(r, uintptr(k)*L))
			wBuf := vBuf[i]

			compInv, compInvS := sc.compInv[i], sc.compInvS[i]
			scFrac := sc.scFrac[i]
			modIn := sc.modIn[i]

			wBuf[0] = num.SMul(w[0], compInv, compInvS, modIn)
			wBuf[1] = num.SMul(w[1], compInv, compInvS, modIn)
			wBuf[2] = num.SMul(w[2], compInv, compInvS, modIn)
			wBuf[3] = num.SMul(w[3], compInv, compInvS, modIn)

			wBuf[4] = num.SMul(w[4], compInv, compInvS, modIn)
			wBuf[5] = num.SMul(w[5], compInv, compInvS, modIn)
			wBuf[6] = num.SMul(w[6], compInv, compInvS, modIn)
			wBuf[7] = num.SMul(w[7], compInv, compInvS, modIn)

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(wBuf[0])), scFrac))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(wBuf[1])), scFrac))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(wBuf[2])), scFrac))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(wBuf[3])), scFrac))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(wBuf[4])), scFrac))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(wBuf[5])), scFrac))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(wBuf[6])), scFrac))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(wBuf[7])), scFrac))
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
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(k)*L))

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
				wBuf := vBuf[j]

				modScInt, modScIntS := modScInt[j], modScIntS[j]

				wOut[0] = num.Add(wOut[0], num.SMul(wBuf[0], modScInt, modScIntS, modOut), modOut)
				wOut[1] = num.Add(wOut[1], num.SMul(wBuf[1], modScInt, modScIntS, modOut), modOut)
				wOut[2] = num.Add(wOut[2], num.SMul(wBuf[2], modScInt, modScIntS, modOut), modOut)
				wOut[3] = num.Add(wOut[3], num.SMul(wBuf[3], modScInt, modScIntS, modOut), modOut)

				wOut[4] = num.Add(wOut[4], num.SMul(wBuf[4], modScInt, modScIntS, modOut), modOut)
				wOut[5] = num.Add(wOut[5], num.SMul(wBuf[5], modScInt, modScIntS, modOut), modOut)
				wOut[6] = num.Add(wOut[6], num.SMul(wBuf[6], modScInt, modScIntS, modOut), modOut)
				wOut[7] = num.Add(wOut[7], num.SMul(wBuf[7], modScInt, modScIntS, modOut), modOut)
			}
		}
	}

	for k := M; k < len(v[0]); k++ {
		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			vBuf[i][0] = num.SMul(v[i][k], sc.compInv[i], sc.compInvS[i], sc.modIn[i])
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(vBuf[i][0])), sc.scFrac[i]))
		}
		uHi[0], uLo[0] = float128.ToUint128(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.Reduce128(uHi[0], uLo[0], sc.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(vBuf[j][0], sc.scInt[i][j], sc.scIntS[i][j], sc.modOut[i]), sc.modOut[i])
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

	f128Pool *sync.Pool
	u64Pool  *sync.Pool
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

		f128Pool: &sync.Pool{
			New: func() any {
				return &[8]float128.Float128{}
			},
		},
		u64Pool: &sync.Pool{
			New: func() any {
				return &[8]uint64{}
			},
		},
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
	if !(e.Rank() == eOut.Rank() && e.ModLen() == len(sc.modIn) && eOut.ModLen() == len(sc.modOut)) {
		panic("input(s) not consistent")
	} else if e.IsNTT {
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
	L := unsafe.Sizeof(uint64(0))

	f128 := sc.f128Pool.Get().(*[8]float128.Float128)
	defer sc.f128Pool.Put(f128)
	uHi := sc.u64Pool.Get().(*[8]uint64)
	defer sc.u64Pool.Put(uHi)
	uLo := sc.u64Pool.Get().(*[8]uint64)
	defer sc.u64Pool.Put(uLo)
	vBuf := make([]*[8]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = sc.u64Pool.Get().(*[8]uint64)
		defer sc.u64Pool.Put(vBuf[i])
	}

	for k := 0; k < M; k += 8 {
		clear(f128[:])
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[8]uint64)(unsafe.Add(r, uintptr(k)*L))
			wBuf := vBuf[i]

			compInv, compInvS, modIn := sc.compInv[i], sc.compInvS[i], sc.modIn[i]
			inv := sc.inv[i]

			wBuf[0] = num.SMul(w[0], compInv, compInvS, modIn)
			wBuf[1] = num.SMul(w[1], compInv, compInvS, modIn)
			wBuf[2] = num.SMul(w[2], compInv, compInvS, modIn)
			wBuf[3] = num.SMul(w[3], compInv, compInvS, modIn)

			wBuf[4] = num.SMul(w[4], compInv, compInvS, modIn)
			wBuf[5] = num.SMul(w[5], compInv, compInvS, modIn)
			wBuf[6] = num.SMul(w[6], compInv, compInvS, modIn)
			wBuf[7] = num.SMul(w[7], compInv, compInvS, modIn)

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(wBuf[0])), inv))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(wBuf[1])), inv))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(wBuf[2])), inv))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(wBuf[3])), inv))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(wBuf[4])), inv))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(wBuf[5])), inv))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(wBuf[6])), inv))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(wBuf[7])), inv))
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
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(k)*L))

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
				wBuf := vBuf[j]

				modScInt, modScIntS := modScInt[j], modScIntS[j]

				wOut[0] = num.Add(wOut[0], num.SMul(wBuf[0], modScInt, modScIntS, modOut), modOut)
				wOut[1] = num.Add(wOut[1], num.SMul(wBuf[1], modScInt, modScIntS, modOut), modOut)
				wOut[2] = num.Add(wOut[2], num.SMul(wBuf[2], modScInt, modScIntS, modOut), modOut)
				wOut[3] = num.Add(wOut[3], num.SMul(wBuf[3], modScInt, modScIntS, modOut), modOut)

				wOut[4] = num.Add(wOut[4], num.SMul(wBuf[4], modScInt, modScIntS, modOut), modOut)
				wOut[5] = num.Add(wOut[5], num.SMul(wBuf[5], modScInt, modScIntS, modOut), modOut)
				wOut[6] = num.Add(wOut[6], num.SMul(wBuf[6], modScInt, modScIntS, modOut), modOut)
				wOut[7] = num.Add(wOut[7], num.SMul(wBuf[7], modScInt, modScIntS, modOut), modOut)
			}
		}

		clear(f128[:])
		for i := 0; i < inLen; i++ {
			wBuf := vBuf[i]

			modScFrac := sc.scFrac[i]

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(wBuf[0])), modScFrac))
			f128[1] = float128.Add(f128[1], float128.Mul(float128.FromInt64(int64(wBuf[1])), modScFrac))
			f128[2] = float128.Add(f128[2], float128.Mul(float128.FromInt64(int64(wBuf[2])), modScFrac))
			f128[3] = float128.Add(f128[3], float128.Mul(float128.FromInt64(int64(wBuf[3])), modScFrac))

			f128[4] = float128.Add(f128[4], float128.Mul(float128.FromInt64(int64(wBuf[4])), modScFrac))
			f128[5] = float128.Add(f128[5], float128.Mul(float128.FromInt64(int64(wBuf[5])), modScFrac))
			f128[6] = float128.Add(f128[6], float128.Mul(float128.FromInt64(int64(wBuf[6])), modScFrac))
			f128[7] = float128.Add(f128[7], float128.Mul(float128.FromInt64(int64(wBuf[7])), modScFrac))
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
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(k)*L))

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
			vBuf[i][0] = num.SMul(v[i][k], sc.compInv[i], sc.compInvS[i], sc.modIn[i])

			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(vBuf[i][0])), sc.inv[i]))
		}

		uLo[0] = float128.ToUint64(f128[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.SMul(sc.modOut[i].Value()-uLo[0], sc.ovfInt[i], sc.ovfIntS[i], sc.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(vBuf[j][0], sc.scInt[i][j], sc.scIntS[i][j], sc.modOut[i]), sc.modOut[i])
			}
		}

		f128[0] = float128.Float128{}
		for i := 0; i < inLen; i++ {
			f128[0] = float128.Add(f128[0], float128.Mul(float128.FromInt64(int64(vBuf[i][0])), sc.scFrac[i]))
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
