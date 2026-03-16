package crt

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// Embedder embeds a vector or polynomial into different modulus.
// In other words, it computes
//
//	[p]_modIn -> [p]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
type Embedder struct {
	// modIn is the input modulus.
	modIn []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus

	// modInv is the accumulated inverse of modIn.
	modInv [][]uint64
	// modInvS is the Shoup form of modIn.
	modInvS [][]uint64

	// base is the mixed-radix basis of modIn.
	base [][]uint64
	// baseS is the Shoup form of base.
	baseS [][]uint64

	// inModOut equals modIn modulo modOut.
	inModOut []uint64

	// idx holds the index of the input modulus limb if it overlaps with the output modulus limb.
	// For example, if modOut[i] = modIn[j], then idx[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	idx []int
}

// NewEmbedder creates a new [Embedder].
func NewEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *Embedder {
	if !isCoprime(modIn) || !isCoprime(modOut) {
		panic("modulus must be coprime")
	}

	modInv := make([][]uint64, len(modIn))
	modInvS := make([][]uint64, len(modIn))
	for i := 0; i < len(modIn); i++ {
		modInv[i] = make([]uint64, len(modIn)-i-1)
		modInvS[i] = make([]uint64, len(modIn)-i-1)
		for j := 0; j < len(modIn)-i-1; j++ {
			modInv[i][j] = num.Inv(modIn[i].Value(), modIn[i+j+1])
			modInvS[i][j] = num.SForm(modInv[i][j], modIn[i+j+1])
		}
	}

	base := make([][]uint64, len(modOut))
	baseS := make([][]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		base[i] = make([]uint64, len(modIn))
		baseS[i] = make([]uint64, len(modIn))

		base[i][0] = 1
		baseS[i][0] = num.SForm(1, modOut[i])
		for j := 1; j < len(modIn); j++ {
			base[i][j] = num.Mul(base[i][j-1], modIn[j-1].Value(), modOut[i])
			baseS[i][j] = num.SForm(base[i][j], modOut[i])
		}
	}

	inModOut := make([]uint64, len(modOut))
	idx := make([]int, len(modOut))
	for i := 0; i < len(modOut); i++ {
		inModOut[i] = 1
		idx[i] = -1
		for j := 0; j < len(modIn); j++ {
			inModOut[i] = num.Mul(inModOut[i], modIn[j].Value(), modOut[i])

			if modOut[i].Value() == modIn[j].Value() {
				idx[i] = j
			}
		}
	}

	return &Embedder{
		modIn:  modIn,
		modOut: modOut,

		modInv:  modInv,
		modInvS: modInvS,

		base:  base,
		baseS: baseS,

		inModOut: inModOut,

		idx: idx,
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

	if inLen == 1 {
		qv := emb.modIn[0].Value()
		halfQv := qv >> 1

		vBuf := embed64Pool.Get().(*[embedBatch]uint64)
		defer embed64Pool.Put(vBuf)

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

	vBuf := make([]*[embedBatch]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get().(*[embedBatch]uint64)
		defer embed64Pool.Put(vBuf[i])
	}

	vBool := embed64Pool.Get().(*[embedBatch]uint64)
	defer embed64Pool.Put(vBool)
	vCorr := embed64Pool.Get().(*[embedBatch]uint64)
	defer embed64Pool.Put(vCorr)

	qLastHalf := emb.modIn[inLen-1].Value() >> 1
	for k := 0; k < M; k += embedBatch {
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))
			copy(vBuf[i][:], w[:])
		}

		for i := 0; i < inLen; i++ {
			wBuf := vBuf[i]
			for j := i + 1; j < inLen; j++ {
				vec.SubTo(vBuf[j][:], vBuf[j][:], wBuf[:], emb.modIn[j])
				vec.SMulScalarTo(vBuf[j][:], vBuf[j][:], emb.modInv[i][j-i-1], emb.modInvS[i][j-i-1], emb.modIn[j])
			}
		}

		vLast := vBuf[inLen-1]
		for i := 0; i < embedBatch; i += 8 {
			vBool[i+0] = (qLastHalf - vLast[i+0]) >> 63
			vBool[i+1] = (qLastHalf - vLast[i+1]) >> 63
			vBool[i+2] = (qLastHalf - vLast[i+2]) >> 63
			vBool[i+3] = (qLastHalf - vLast[i+3]) >> 63

			vBool[i+4] = (qLastHalf - vLast[i+4]) >> 63
			vBool[i+5] = (qLastHalf - vLast[i+5]) >> 63
			vBool[i+6] = (qLastHalf - vLast[i+6]) >> 63
			vBool[i+7] = (qLastHalf - vLast[i+7]) >> 63
		}

		for i := 0; i < outLen; i++ {
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))
			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				copy(wOut[:], vBuf[emb.idx[i]][:])
			} else {
				base, baseS := emb.base[i], emb.baseS[i]
				inModOut := emb.inModOut[i]
				modOut := emb.modOut[i]

				vec.SMulScalarTo(wOut[:], vBuf[0][:], base[0], baseS[0], modOut)
				for j := 1; j < inLen; j++ {
					vec.SMulAddScalarTo(wOut[:], vBuf[j][:], base[j], baseS[j], modOut)
				}
				vec.MulScalarTo(vCorr[:], vBool[:], inModOut, nil)
				vec.SubTo(wOut[:], wOut[:], vCorr[:], modOut)
			}
		}
	}

	for i := 0; i < inLen; i++ {
		copy(vBuf[i][:], v[i][M:])
	}

	for i := 0; i < inLen; i++ {
		wBuf := vBuf[i][:len(v[0])-M]
		for j := i + 1; j < inLen; j++ {
			vec.SubTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], wBuf[:], emb.modIn[j])
			vec.SMulScalarTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], emb.modInv[i][j-i-1], emb.modInvS[i][j-i-1], emb.modIn[j])
		}
	}

	vLast := vBuf[inLen-1]
	for i := 0; i < len(v[0])-M; i++ {
		vBool[i] = (qLastHalf - vLast[i]) >> 63
	}

	for i := 0; i < outLen; i++ {
		wOut := vOut[i][M:]
		if 0 <= emb.idx[i] && emb.idx[i] < inLen {
			copy(wOut[:], vBuf[emb.idx[i]][:len(v[0])-M])
		} else {
			base, baseS := emb.base[i], emb.baseS[i]
			inModOut := emb.inModOut[i]
			modOut := emb.modOut[i]

			vec.SMulScalarTo(wOut[:], vBuf[0][:len(v[0])-M], base[0], baseS[0], modOut)
			for j := 1; j < inLen; j++ {
				vec.SMulAddScalarTo(wOut[:], vBuf[j][:len(v[0])-M], base[j], baseS[j], modOut)
			}
			vec.MulScalarTo(vCorr[:len(v[0])-M], vBool[:len(v[0])-M], inModOut, nil)
			vec.SubTo(wOut[:], wOut[:], vCorr[:len(v[0])-M], modOut)
		}
	}
}

// ModulusIn returns the input modulus.
func (emb *ApproxEmbedder) ModulusIn() []*num.Modulus {
	return emb.modIn
}

// ModulusOut returns the output modulus.
func (emb *ApproxEmbedder) ModulusOut() []*num.Modulus {
	return emb.modOut
}
