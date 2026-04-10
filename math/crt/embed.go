package crt

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/internal/pool"
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
	embed64Pool = pool.NewPool(func() *[embedBatch]uint64 {
		return new([embedBatch]uint64)
	})
)

// embedToModOut returns sign(x) mod qOut for x in [0, qIn).
func embedToModOut(x uint64, qOut *num.Modulus, qIn, halfQIn uint64) uint64 {
	if x <= halfQIn {
		return num.Reduce(x, qOut)
	}
	return num.Neg(num.Reduce(qIn-x, qOut), qOut)
}

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

		vBuf := embed64Pool.Get()
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
					wOut[0] = embedToModOut(vBuf[0], modOut, qv, halfQv)
					wOut[1] = embedToModOut(vBuf[1], modOut, qv, halfQv)
					wOut[2] = embedToModOut(vBuf[2], modOut, qv, halfQv)
					wOut[3] = embedToModOut(vBuf[3], modOut, qv, halfQv)

					wOut[4] = embedToModOut(vBuf[4], modOut, qv, halfQv)
					wOut[5] = embedToModOut(vBuf[5], modOut, qv, halfQv)
					wOut[6] = embedToModOut(vBuf[6], modOut, qv, halfQv)
					wOut[7] = embedToModOut(vBuf[7], modOut, qv, halfQv)
				}
			}
		}

		for i := M; i < len(v[0]); i++ {
			c := v[0][i]
			for j := 0; j < outLen; j++ {
				if emb.idx[j] == 0 {
					vOut[j][i] = c
				} else {
					vOut[j][i] = embedToModOut(c, emb.modOut[j], qv, halfQv)
				}
			}
		}

		return
	}

	vBuf := make([]*[embedBatch]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get()
		defer embed64Pool.Put(vBuf[i])
	}

	vBool := embed64Pool.Get()
	defer embed64Pool.Put(vBool)
	vCorr := embed64Pool.Get()
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

		vBuf := embed64Pool.Get()
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
					wOut[0] = embedToModOut(vBuf[0], emb.modOut[j], qv, halfQv)
					wOut[1] = embedToModOut(vBuf[1], emb.modOut[j], qv, halfQv)
					wOut[2] = embedToModOut(vBuf[2], emb.modOut[j], qv, halfQv)
					wOut[3] = embedToModOut(vBuf[3], emb.modOut[j], qv, halfQv)

					wOut[4] = embedToModOut(vBuf[4], emb.modOut[j], qv, halfQv)
					wOut[5] = embedToModOut(vBuf[5], emb.modOut[j], qv, halfQv)
					wOut[6] = embedToModOut(vBuf[6], emb.modOut[j], qv, halfQv)
					wOut[7] = embedToModOut(vBuf[7], emb.modOut[j], qv, halfQv)
				}
			}
		}

		for i := M; i < len(v[0]); i++ {
			c := v[0][i]
			for j := 0; j < outLen; j++ {
				if emb.idx[j] == 0 {
					vOut[j][i] = c
				} else {
					vOut[j][i] = embedToModOut(c, emb.modOut[j], qv, halfQv)
				}
			}
		}

		return
	}

	vBuf := make([]*[embedBatch]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get()
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
			vBuf[i] = embed64Pool.Get()
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
		vMul[i] = embed64Pool.Get()
		defer embed64Pool.Put(vMul[i])
	}
	vBuf := make([]*[embedBatch]uint64, inLen-gcdLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get()
		defer embed64Pool.Put(vBuf[i])
	}

	vBool := embed64Pool.Get()
	defer embed64Pool.Put(vBool)
	vCorr := embed64Pool.Get()
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
					wOut[j+0] = embedToModOut(wBuf[j+0], modOut, qLast, qLastHalf)
					wOut[j+1] = embedToModOut(wBuf[j+1], modOut, qLast, qLastHalf)
					wOut[j+2] = embedToModOut(wBuf[j+2], modOut, qLast, qLastHalf)
					wOut[j+3] = embedToModOut(wBuf[j+3], modOut, qLast, qLastHalf)

					wOut[j+4] = embedToModOut(wBuf[j+4], modOut, qLast, qLastHalf)
					wOut[j+5] = embedToModOut(wBuf[j+5], modOut, qLast, qLastHalf)
					wOut[j+6] = embedToModOut(wBuf[j+6], modOut, qLast, qLastHalf)
					wOut[j+7] = embedToModOut(wBuf[j+7], modOut, qLast, qLastHalf)
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
				vOut[i][j+M] = embedToModOut(vBuf[0][j], sc.modOut[i], qLast, qLastHalf)
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
