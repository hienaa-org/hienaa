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

var (
	embed64Pool = sync.Pool{
		New: func() any {
			return new([embedBatch]uint64)
		},
	}
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

	// modInv is the inverse of modIn.
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

	inLen, outLen := len(v), len(vOut)
	if inLen > len(emb.modIn) || len(vOut) > len(emb.modOut) {
		panic("input(s) not consistent")
	}

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

				modOut := emb.modOut[j]

				if emb.idx[j] == 0 {
					copy(wOut[:], vBuf[:])
				} else {
					wOut[0] = reduceModInToModOutSigned(vBuf[0], modOut, qv, halfQv)
					wOut[1] = reduceModInToModOutSigned(vBuf[1], modOut, qv, halfQv)
					wOut[2] = reduceModInToModOutSigned(vBuf[2], modOut, qv, halfQv)
					wOut[3] = reduceModInToModOutSigned(vBuf[3], modOut, qv, halfQv)

					wOut[4] = reduceModInToModOutSigned(vBuf[4], modOut, qv, halfQv)
					wOut[5] = reduceModInToModOutSigned(vBuf[5], modOut, qv, halfQv)
					wOut[6] = reduceModInToModOutSigned(vBuf[6], modOut, qv, halfQv)
					wOut[7] = reduceModInToModOutSigned(vBuf[7], modOut, qv, halfQv)
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

		for i := 0; i < outLen; i++ {
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))
			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				copy(wOut[:], vBuf[emb.idx[i]][:])
			}
		}

		for i := 0; i < inLen; i++ {
			wBuf := vBuf[i]
			for j := i + 1; j < inLen; j++ {
				modIn := emb.modIn[j]
				vec.SubTo(vBuf[j][:], vBuf[j][:], wBuf[:], modIn)
				vec.SMulScalarTo(vBuf[j][:], vBuf[j][:], emb.modInv[i][j-i-1], emb.modInvS[i][j-i-1], modIn)
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
			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				continue
			}

			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))

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

	for i := 0; i < inLen; i++ {
		copy(vBuf[i][:len(v[0])-M], v[i][M:])
	}

	for i := 0; i < outLen; i++ {
		if 0 <= emb.idx[i] && emb.idx[i] < inLen {
			copy(vOut[i][M:], v[emb.idx[i]][M:])
		}
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
		if 0 <= emb.idx[i] && emb.idx[i] < inLen {
			continue
		}

		wOut := vOut[i][M:]
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

// ModulusIn returns the input modulus.
func (emb *Embedder) ModulusIn() []*num.Modulus {
	return emb.modIn
}

// ModulusOut returns the output modulus.
func (emb *Embedder) ModulusOut() []*num.Modulus {
	return emb.modOut
}

// ApproxEmbedder embeds a vector or polynomial into different modulus,
// without removing the overflow term.
// In other words, it computes
//
//	[p]_modIn -> [p + modIn * I]_modOut
type ApproxEmbedder struct {
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

	// idx holds the index of the input modulus limb if it overlaps with the output modulus limb.
	// For example, if modOut[i] = modIn[j], then idx[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	idx []int
}

// NewApproxEmbedder creates a new [ApproxEmbedder].
func NewApproxEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *ApproxEmbedder {
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

	return &ApproxEmbedder{
		modIn:  modIn,
		modOut: modOut,

		compInv:  compInv,
		compInvS: compInvS,

		comp:  comp,
		compS: compS,

		idx: idx,
	}
}

// Embed returns the embedding of e to the output modulus.
// If e.ModLen() < len(emb.modIn), it only embeds the first e.ModLen() elements.
func (emb *ApproxEmbedder) Embed(e *Element) *Element {
	eOut := NewPoly(e.Rank(), e.ModLen())
	emb.EmbedTo(eOut, e)
	return eOut
}

// EmbedTo embeds e to eOut.
// If e.ModLen() < len(emb.modIn) or eOut.ModLen() < len(emb.modOut),
// it only embeds the first emb.ModLen() elements to eOut.ModLen() elements.
func (emb *ApproxEmbedder) EmbedTo(eOut, e *Element) {
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
func (emb *ApproxEmbedder) EmbedVec(v [][]uint64) [][]uint64 {
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
func (emb *ApproxEmbedder) EmbedVecTo(vOut, v [][]uint64) {
	M := (len(v[0]) >> logEmbedBatch) << logEmbedBatch
	L := unsafe.Sizeof(uint64(0))

	inLen, outLen := len(v), len(vOut)

	if inLen != len(emb.modIn) || len(vOut) > len(emb.modOut) {
		panic("input(s) not consistent")
	}

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

	for k := 0; k < M; k += embedBatch {
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))
			copy(vBuf[i][:], w[:])
		}

		for i := 0; i < outLen; i++ {
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))

			if 0 <= emb.idx[i] && emb.idx[i] < inLen {
				copy(wOut[:], v[emb.idx[i]][k:k+embedBatch])
			} else {
				modOut := emb.modOut[i]
				comp, compS := emb.comp[i], emb.compS[i]
				for j := 0; j < inLen; j++ {
					wBuf := vBuf[j]
					comp, compS := comp[j], compS[j]

					vec.SMulAddScalarTo(wOut[:], wBuf[:], comp, compS, modOut)
				}
			}
		}
	}

	for i := 0; i < inLen; i++ {
		copy(vBuf[i][:], v[i][M:])
	}

	for i := 0; i < outLen; i++ {
		wOut := vOut[i][M:]

		if 0 <= emb.idx[i] && emb.idx[i] < inLen {
			copy(wOut, v[emb.idx[i]][M:])
		} else {
			modOut := emb.modOut[i]
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
func (emb *ApproxEmbedder) ModulusIn() []*num.Modulus {
	return emb.modIn
}

// ModulusOut returns the output modulus.
func (emb *ApproxEmbedder) ModulusOut() []*num.Modulus {
	return emb.modOut
}

// Scaler scales a polynomial to different modulus.
// In other words, it computes
//
//	[p]_modIn -> [(modOut / modIn) * p]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
type Scaler struct {
	// modIn is the input modulus.
	modIn []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus
	// modInComp equals modIn/modGCD.
	modInComp []*num.Modulus

	// modInCompLastHalf is the half of modInComp[-1].
	modInCompLastHalf uint64
	// modGCDLen is the length of the GCD of modIn and modOut.
	modGCDLen int

	// idxInToGCD maps the index of modIn to the index of modGCD.
	// If modIn[i] is not in modGCD, idxInToGCD[i] is -1.
	idxInToGCD []int
	// idxInToComp maps the index of modIn to the index of modInComp.
	// If modIn[i] is in modGCD, idxInToComp[i] is -1.
	idxInToComp []int
	// idxOutToGCD maps the index of modOut to the index of modGCD.
	// If modOut[i] is not in modGCD, idxOutToGCD[i] is -1.
	idxOutToGCD []int

	// modInCompInv is the inverse of of modInComp.
	modInCompInv [][]uint64
	// modInCompInvS is the Shoup form of modInCompInv.
	modInCompInvS [][]uint64

	// base is the mixed-radix basis of modInComp.
	base [][]uint64
	// baseS is the Shoup form of base.
	baseS [][]uint64

	// inCompModOut equals modInComp modulo modOut.
	inCompModOut []uint64

	// emb is the internal embedder of this scaler.
	// Embeds modIn/modGCD -> modOut.
	emb *Embedder

	// outCompModInComp equals modOut/modGCD modulo modIn.
	outCompModIn  []uint64
	outCompModInS []uint64
	// negInCompInvModOut equals the inverse of modIn/modGCD modulo modOut.
	negInCompInvModOut  []uint64
	negInCompInvModOutS []uint64
}

func NewScaler(modOut, modIn []*num.Modulus) *Scaler {
	if !isCoprime(modIn) || !isCoprime(modOut) {
		panic("modulus must be coprime")
	}

	modGCD := make([]*num.Modulus, 0, min(len(modOut), len(modIn)))
	isGCDOut := make([]bool, len(modOut))
	isGCDIn := make([]bool, len(modIn))
	for i := 0; i < len(modOut); i++ {
		for j := 0; j < len(modIn); j++ {
			if modOut[i].Value() == modIn[j].Value() {
				modGCD = append(modGCD, modOut[i])
				isGCDOut[i] = true
				isGCDIn[j] = true
				break
			}
		}
	}

	var idxIn int
	idxInToGCD := make([]int, len(modIn))
	idxInToComp := make([]int, len(modIn))
	for i := 0; i < len(modIn); i++ {
		idxInToGCD[i] = -1
		idxInToComp[i] = -1
		for j := 0; j < len(modGCD); j++ {
			if modIn[i].Value() == modGCD[j].Value() {
				idxInToGCD[i] = j
				break
			}
		}
		if idxInToGCD[i] == -1 {
			idxInToComp[i] = idxIn
			idxIn++
		}
	}

	idxOutToGCD := make([]int, len(modOut))
	for i := 0; i < len(modOut); i++ {
		idxOutToGCD[i] = -1
		for j := 0; j < len(modGCD); j++ {
			if modOut[i].Value() == modGCD[j].Value() {
				idxOutToGCD[i] = j
				break
			}
		}
	}

	modInComp := make([]*num.Modulus, 0, len(modIn)-len(modGCD))
	for i := 0; i < len(modIn); i++ {
		if idxInToGCD[i] == -1 {
			modInComp = append(modInComp, modIn[i])
		}
	}

	modCompInv := make([][]uint64, len(modInComp))
	modCompInvS := make([][]uint64, len(modInComp))
	for i := 0; i < len(modInComp); i++ {
		modCompInv[i] = make([]uint64, len(modInComp)-i-1)
		modCompInvS[i] = make([]uint64, len(modInComp)-i-1)
		for j := 0; j < len(modInComp)-i-1; j++ {
			modCompInv[i][j] = num.Inv(modInComp[i].Value(), modInComp[i+j+1])
			modCompInvS[i][j] = num.SForm(modCompInv[i][j], modInComp[i+j+1])
		}
	}

	base := make([][]uint64, len(modOut))
	baseS := make([][]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		base[i] = make([]uint64, len(modInComp))
		baseS[i] = make([]uint64, len(modInComp))
		for j := 0; j < len(modInComp); j++ {
			if j == 0 {
				base[i][j] = 1
			} else {
				base[i][j] = num.Mul(base[i][j-1], modInComp[j-1].Value(), modOut[i])
			}
			baseS[i][j] = num.SForm(base[i][j], modOut[i])
		}
	}

	inCompModOut := make([]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		inCompModOut[i] = 1
		for j := 0; j < len(modInComp); j++ {
			inCompModOut[i] = num.Mul(inCompModOut[i], modInComp[j].Value(), modOut[i])
		}
	}

	outCompModIn := make([]uint64, len(modIn))
	outCompModInS := make([]uint64, len(modIn))
	for i := 0; i < len(modIn); i++ {
		outCompModIn[i] = 1
		for j := 0; j < len(modOut); j++ {
			if !isGCDOut[j] {
				outCompModIn[i] = num.Mul(outCompModIn[i], modOut[j].Value(), modIn[i])
			}
		}
		outCompModInS[i] = num.SForm(outCompModIn[i], modIn[i])
	}

	negInCompInvModOut := make([]uint64, len(modOut))
	negInCompInvModOutS := make([]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		negInCompInvModOut[i] = 1
		for j := 0; j < len(modIn); j++ {
			if !isGCDIn[j] {
				negInCompInvModOut[i] = num.Mul(negInCompInvModOut[i], num.Inv(modIn[j].Value(), modOut[i]), modOut[i])
			}
		}
		negInCompInvModOut[i] = num.Neg(negInCompInvModOut[i], modOut[i])
		negInCompInvModOutS[i] = num.SForm(negInCompInvModOut[i], modOut[i])
	}

	return &Scaler{
		modIn:     modIn,
		modOut:    modOut,
		modInComp: modInComp,

		modGCDLen: len(modIn) - len(modInComp),

		idxInToGCD:  idxInToGCD,
		idxInToComp: idxInToComp,
		idxOutToGCD: idxOutToGCD,

		modInCompInv:  modCompInv,
		modInCompInvS: modCompInvS,

		base:  base,
		baseS: baseS,

		inCompModOut: inCompModOut,

		outCompModIn:        outCompModIn,
		outCompModInS:       outCompModInS,
		negInCompInvModOut:  negInCompInvModOut,
		negInCompInvModOutS: negInCompInvModOutS,
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

	inLen, gcdLen, outLen := len(sc.modIn), sc.modGCDLen, len(sc.modOut)
	M := (len(v[0]) >> logEmbedBatch) << logEmbedBatch
	L := unsafe.Sizeof(uint64(0))

	if len(sc.modIn) == sc.modGCDLen {
		vBuf := make([]*[embedBatch]uint64, inLen)
		for i := range vBuf {
			vBuf[i] = embed64Pool.Get().(*[embedBatch]uint64)
			defer embed64Pool.Put(vBuf[i])
		}

		for k := 0; k < M; k += embedBatch {
			for i := 0; i < inLen; i++ {
				r := unsafe.Pointer(unsafe.SliceData(v[i]))
				w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))

				ii := sc.idxInToGCD[i]
				vec.SMulScalarTo(vBuf[ii][:], w[:], sc.outCompModIn[ii], sc.outCompModInS[ii], sc.modOut[i])
			}

			for i := 0; i < outLen; i++ {
				rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
				wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))

				ii := sc.idxOutToGCD[i]
				if ii >= 0 {
					copy(wOut[:], vBuf[ii][:])
				} else {
					clear(wOut[:])
				}
			}
		}

		for k := M; k < len(v[0]); k++ {
			for i := 0; i < inLen; i++ {
				w := v[i][M:]

				ii := sc.idxInToGCD[i]
				vec.SMulScalarTo(vBuf[ii][:len(v[0])-M], w[:], sc.outCompModIn[ii], sc.outCompModInS[ii], sc.modOut[i])
			}

			for i := 0; i < outLen; i++ {
				wOut := vOut[i][M:]

				ii := sc.idxOutToGCD[i]
				if ii >= 0 {
					copy(wOut[:], vBuf[ii][:len(v[0])-M])
				} else {
					clear(wOut[:])
				}
			}
		}

		return
	}

	vMul := make([]*[embedBatch]uint64, gcdLen)
	for i := range vMul {
		vMul[i] = embed64Pool.Get().(*[embedBatch]uint64)
		defer embed64Pool.Put(vMul[i])
	}
	vBuf := make([]*[embedBatch]uint64, inLen-gcdLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get().(*[embedBatch]uint64)
		defer embed64Pool.Put(vBuf[i])
	}

	vBool := embed64Pool.Get().(*[embedBatch]uint64)
	defer embed64Pool.Put(vBool)
	vCorr := embed64Pool.Get().(*[embedBatch]uint64)
	defer embed64Pool.Put(vCorr)

	var qLast uint64
	for i := inLen - 1; i >= 0; i-- {
		if sc.idxInToComp[i] >= 0 {
			qLast = sc.modIn[i].Value()
			break
		}
	}
	qLastHalf := qLast >> 1

	for k := 0; k < M; k += embedBatch {
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))

			modIn := sc.modIn[i]

			ii := sc.idxInToGCD[i]
			if ii >= 0 {
				vec.SMulScalarTo(vMul[ii][:], w[:], sc.outCompModIn[i], sc.outCompModInS[i], modIn)
			}
		}

		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[i]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))

			modIn := sc.modIn[i]

			ii := sc.idxInToComp[i]
			if ii >= 0 {
				vec.SMulScalarTo(vBuf[ii][:], w[:], sc.outCompModIn[i], sc.outCompModInS[i], modIn)
			}
		}

		if inLen-gcdLen == 1 {
			for i := 0; i < outLen; i++ {
				rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
				wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))
				wBuf := vBuf[0]

				modOut := sc.modOut[i]

				for j := 0; j < embedBatch; j += 8 {
					wOut[j+0] = reduceModInToModOutSigned(wBuf[j+0], modOut, qLast, qLastHalf)
					wOut[j+1] = reduceModInToModOutSigned(wBuf[j+1], modOut, qLast, qLastHalf)
					wOut[j+2] = reduceModInToModOutSigned(wBuf[j+2], modOut, qLast, qLastHalf)
					wOut[j+3] = reduceModInToModOutSigned(wBuf[j+3], modOut, qLast, qLastHalf)

					wOut[j+4] = reduceModInToModOutSigned(wBuf[j+4], modOut, qLast, qLastHalf)
					wOut[j+5] = reduceModInToModOutSigned(wBuf[j+5], modOut, qLast, qLastHalf)
					wOut[j+6] = reduceModInToModOutSigned(wBuf[j+6], modOut, qLast, qLastHalf)
					wOut[j+7] = reduceModInToModOutSigned(wBuf[j+7], modOut, qLast, qLastHalf)
				}
			}
		} else {
			for i := 0; i < inLen-gcdLen; i++ {
				wBuf := vBuf[i]
				for j := i + 1; j < inLen-gcdLen; j++ {
					modInComp := sc.modInComp[j]
					vec.SubTo(vBuf[j][:], vBuf[j][:], wBuf[:], modInComp)
					vec.SMulScalarTo(vBuf[j][:], vBuf[j][:], sc.modInCompInv[i][j-i-1], sc.modInCompInvS[i][j-i-1], modInComp)
				}
			}

			vLast := vBuf[inLen-gcdLen-1]
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

				base, baseS := sc.base[i], sc.baseS[i]
				inCompModOut := sc.inCompModOut[i]
				modOut := sc.modOut[i]

				vec.SMulScalarTo(wOut[:], vBuf[0][:], base[0], baseS[0], modOut)
				for j := 1; j < inLen-gcdLen; j++ {
					vec.SMulAddScalarTo(wOut[:], vBuf[j][:], base[j], baseS[j], modOut)
				}
				vec.MulScalarTo(vCorr[:], vBool[:], inCompModOut, nil)
				vec.SubTo(wOut[:], wOut[:], vCorr[:], modOut)
			}
		}

		for i := 0; i < outLen; i++ {
			rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
			wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))

			modOut := sc.modOut[i]

			ii := sc.idxOutToGCD[i]
			if ii >= 0 {
				vec.SubTo(wOut[:], wOut[:], vMul[ii][:], modOut)
			}
			vec.SMulScalarTo(wOut[:], wOut[:], sc.negInCompInvModOut[i], sc.negInCompInvModOutS[i], modOut)
		}
	}

	for i := 0; i < inLen; i++ {
		w := v[i][M:]

		ii := sc.idxInToGCD[i]
		if ii >= 0 {
			vec.SMulScalarTo(vMul[ii][:len(v[0])-M], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modIn[i])
		}
	}

	for i := 0; i < inLen; i++ {
		w := v[i][M:]

		ii := sc.idxInToComp[i]
		if ii >= 0 {
			vec.SMulScalarTo(vBuf[ii][:len(v[0])-M], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modIn[i])
		}
	}

	if inLen-gcdLen == 1 {
		for i := 0; i < outLen; i++ {
			for j := 0; j < len(v[0])-M; j++ {
				vOut[i][j+M] = reduceModInToModOutSigned(vBuf[0][j], sc.modOut[i], qLast, qLastHalf)
			}
		}
	} else {
		for i := 0; i < inLen-gcdLen; i++ {
			wBuf := vBuf[i][:len(v[0])-M]
			for j := i + 1; j < inLen-gcdLen; j++ {
				vec.SubTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], wBuf[:], sc.modInComp[j])
				vec.SMulScalarTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], sc.modInCompInv[i][j-i-1], sc.modInCompInvS[i][j-i-1], sc.modInComp[j])
			}
		}

		vLast := vBuf[inLen-gcdLen-1]
		for i := 0; i < len(v[0])-M; i++ {
			vBool[i] = (qLastHalf - vLast[i]) >> 63
		}

		for i := 0; i < outLen; i++ {
			wOut := vOut[i][M:]

			base, baseS := sc.base[i], sc.baseS[i]
			inCompModOut := sc.inCompModOut[i]
			modOut := sc.modOut[i]

			vec.SMulScalarTo(wOut[:], vBuf[0][:len(v[0])-M], base[0], baseS[0], modOut)
			for j := 1; j < inLen-gcdLen; j++ {
				vec.SMulAddScalarTo(wOut[:], vBuf[j][:len(v[0])-M], base[j], baseS[j], modOut)
			}
			vec.MulScalarTo(vCorr[:len(v[0])-M], vBool[:len(v[0])-M], inCompModOut, nil)
			vec.SubTo(wOut[:], wOut[:], vCorr[:len(v[0])-M], modOut)
		}
	}

	for i := 0; i < outLen; i++ {
		wOut := vOut[i][M:]

		modOut := sc.modOut[i]

		ii := sc.idxOutToGCD[i]
		if ii >= 0 {
			vec.SubTo(wOut[:], wOut[:], vMul[ii][:len(v[0])-M], modOut)
		}
		vec.SMulScalarTo(wOut[:], wOut[:], sc.negInCompInvModOut[i], sc.negInCompInvModOutS[i], modOut)
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
