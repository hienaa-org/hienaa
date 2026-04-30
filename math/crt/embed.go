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

// Embedder embeds a polynomial into different modulus.
type Embedder struct {
	opIn  *Operator
	opOut *Operator

	vecEmb *VecEmbedder
	pool   *pool.Pool[*[]uint64]
}

// NewEmbedder creates a new [Embedder].
func NewEmbedder(opOut *Operator, opIn *Operator) *Embedder {
	return &Embedder{
		opIn:  opIn,
		opOut: opOut,

		vecEmb: NewVecEmbedder(opOut.mod, opIn.mod),
		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, opIn.params.Rank())
			return &v
		}),
	}
}

// WithPool sets the pool of [Embedder] and returns it.
func (emb *Embedder) WithPool(p *pool.Pool[*[]uint64]) *Embedder {
	emb.pool = p
	return emb
}

// Embed embeds e to the output modulus.
func (emb *Embedder) Embed(e *Element, isNTT bool) *Element {
	eOut := NewPoly(e.Rank(), len(emb.opOut.mod))
	emb.EmbedTo(eOut, e, isNTT)
	return eOut
}

// EmbedTo embeds e to eOut.
func (emb *Embedder) EmbedTo(eOut, e *Element, isNTT bool) {
	if eOut.Rank() != e.Rank() {
		panic("input(s) not consistent")
	}

	switch e.Type() {
	case TypeScalar:
		emb.vecEmb.EmbedTo(eOut.Coeffs, e.Coeffs)

	case TypePoly:
		if e.Rank() != emb.opIn.params.Rank() {
			panic("input(s) not consistent")
		}

		if e.IsNTT {
			if isNTT {
				emb.embedNTTTo(eOut, e)
				eOut.IsNTT = true
			} else {
				eInvNTT, put := getElementFromPool(emb.pool, e.ModLen(), false)
				defer put()

				emb.opIn.InvNTTTo(eInvNTT, e)
				emb.vecEmb.EmbedTo(eOut.Coeffs, eInvNTT.Coeffs)
				eOut.IsNTT = false
			}
		} else {
			emb.vecEmb.EmbedTo(eOut.Coeffs, e.Coeffs)
			eOut.IsNTT = false
			if isNTT {
				emb.opOut.Slice(0, eOut.ModLen()).FwdNTTTo(eOut, eOut)
			}
		}
	}
}

// embedNTTTo embeds eNTT to eOutNTT.
func (emb *Embedder) embedNTTTo(eOutNTT, eNTT *Element) {
	M := (len(eNTT.Coeffs[0]) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	inLen, outLen := len(eNTT.Coeffs), len(eOutNTT.Coeffs)
	if inLen != len(emb.vecEmb.modIn) || outLen > len(emb.vecEmb.modOut) {
		panic("input(s) not consistent")
	}

	eNTTCopy, put := getElementFromPool(emb.pool, eNTT.ModLen(), false)
	defer put()
	eNTTCopy.CopyFrom(eNTT)

	eInvNTT, put := getElementFromPool(emb.pool, eNTT.ModLen(), false)
	defer put()
	emb.opIn.InvNTTTo(eInvNTT, eNTT)

	if inLen == 1 {
		qv := emb.vecEmb.modIn[0].Value()
		halfQv := qv >> 1

		r := unsafe.Pointer(unsafe.SliceData(eInvNTT.Coeffs[0]))

		for i := 0; i < M; i += 8 {
			wIn := (*[8]uint64)(unsafe.Add(r, uintptr(i)*L))

			for j := 0; j < outLen; j++ {
				rOut := unsafe.Pointer(unsafe.SliceData(eOutNTT.Coeffs[j]))
				wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))

				modOut := emb.vecEmb.modOut[j]

				if emb.vecEmb.idx[j] != 0 {
					wOut[0] = embedToModOut(wIn[0], modOut, qv, halfQv)
					wOut[1] = embedToModOut(wIn[1], modOut, qv, halfQv)
					wOut[2] = embedToModOut(wIn[2], modOut, qv, halfQv)
					wOut[3] = embedToModOut(wIn[3], modOut, qv, halfQv)

					wOut[4] = embedToModOut(wIn[4], modOut, qv, halfQv)
					wOut[5] = embedToModOut(wIn[5], modOut, qv, halfQv)
					wOut[6] = embedToModOut(wIn[6], modOut, qv, halfQv)
					wOut[7] = embedToModOut(wIn[7], modOut, qv, halfQv)
				}
			}
		}

		for i := M; i < len(eInvNTT.Coeffs[0]); i++ {
			c := eInvNTT.Coeffs[0][i]
			for j := 0; j < outLen; j++ {
				if emb.vecEmb.idx[j] != 0 {
					eOutNTT.Coeffs[j][i] = embedToModOut(c, emb.vecEmb.modOut[j], qv, halfQv)
				}
			}
		}

		for i := 0; i < outLen; i++ {
			if emb.vecEmb.idx[i] == 0 {
				copy(eOutNTT.Coeffs[i][:], eNTTCopy.Coeffs[0][:])
			} else if emb.opOut.ntt[i] != nil {
				emb.opOut.ntt[i].ForwardTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i])
			}
		}

		return
	}

	vBoolPtr := emb.pool.Get()
	vBool := *vBoolPtr
	defer emb.pool.Put(vBoolPtr)

	vCorrPtr := emb.pool.Get()
	vCorr := *vCorrPtr
	defer emb.pool.Put(vCorrPtr)

	for i := 0; i < inLen; i++ {
		for j := i + 1; j < inLen; j++ {
			modIn := emb.vecEmb.modIn[j]
			vec.SubTo(eInvNTT.Coeffs[j], eInvNTT.Coeffs[j], eInvNTT.Coeffs[i], modIn)
			vec.SMulScalarTo(eInvNTT.Coeffs[j], eInvNTT.Coeffs[j], emb.vecEmb.modInv[i][j-i-1], emb.vecEmb.modInvS[i][j-i-1], modIn)
		}
	}

	clear(vBool[:])
	for i := 0; i < M; i += 8 {
		vBool[i+0] = isMixedRadixNegative(eInvNTT.Coeffs, i+0, emb.vecEmb.modInHalf)
		vBool[i+1] = isMixedRadixNegative(eInvNTT.Coeffs, i+1, emb.vecEmb.modInHalf)
		vBool[i+2] = isMixedRadixNegative(eInvNTT.Coeffs, i+2, emb.vecEmb.modInHalf)
		vBool[i+3] = isMixedRadixNegative(eInvNTT.Coeffs, i+3, emb.vecEmb.modInHalf)

		vBool[i+4] = isMixedRadixNegative(eInvNTT.Coeffs, i+4, emb.vecEmb.modInHalf)
		vBool[i+5] = isMixedRadixNegative(eInvNTT.Coeffs, i+5, emb.vecEmb.modInHalf)
		vBool[i+6] = isMixedRadixNegative(eInvNTT.Coeffs, i+6, emb.vecEmb.modInHalf)
		vBool[i+7] = isMixedRadixNegative(eInvNTT.Coeffs, i+7, emb.vecEmb.modInHalf)
	}
	for i := M; i < len(eNTT.Coeffs[0]); i++ {
		vBool[i] = isMixedRadixNegative(eInvNTT.Coeffs, i, emb.vecEmb.modInHalf)
	}

	for i := 0; i < outLen; i++ {
		if 0 <= emb.vecEmb.idx[i] && emb.vecEmb.idx[i] < inLen {
			continue
		}

		base, baseS := emb.vecEmb.base[i], emb.vecEmb.baseS[i]
		inModOut := emb.vecEmb.inModOut[i]
		modOut := emb.vecEmb.modOut[i]

		vec.SMulScalarTo(eOutNTT.Coeffs[i], eInvNTT.Coeffs[0], base[0], baseS[0], modOut)
		for j := 1; j < inLen; j++ {
			vec.SMulAddScalarTo(eOutNTT.Coeffs[i], eInvNTT.Coeffs[j], base[j], baseS[j], modOut)
		}
		vec.MulScalarTo(vCorr, vBool, inModOut, nil)
		vec.SubTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i], vCorr, modOut)
	}

	for i := 0; i < outLen; i++ {
		if 0 <= emb.vecEmb.idx[i] && emb.vecEmb.idx[i] < inLen {
			copy(eOutNTT.Coeffs[i], eNTTCopy.Coeffs[emb.vecEmb.idx[i]])
		} else if emb.opOut.ntt[i] != nil {
			emb.opOut.ntt[i].ForwardTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i])
		}
	}

	eOutNTT.IsNTT = true
}

// OperatorIn returns the input operator.
func (emb *Embedder) OperatorIn() *Operator {
	return emb.opIn
}

// OperatorOut returns the output operator.
func (emb *Embedder) OperatorOut() *Operator {
	return emb.opOut
}

// VecEmbedder embeds a vector or polynomial into different modulus.
// In other words, it computes
//
//	[p]_modIn -> [p]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
type VecEmbedder struct {
	// modIn is the input modulus.
	modIn []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus

	// modInHalf is half of modIn.
	modInHalf []uint64

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

// NewVecEmbedder creates a new [VecEmbedder].
func NewVecEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *VecEmbedder {
	if !isCoprime(modIn) || !isCoprime(modOut) {
		panic("modulus must be coprime")
	}

	modInHalf := make([]uint64, len(modIn))
	for i := 0; i < len(modIn); i++ {
		modInHalf[i] = modIn[i].Value() >> 1
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

	return &VecEmbedder{
		modIn:  modIn,
		modOut: modOut,

		modInHalf: modInHalf,

		modInv:  modInv,
		modInvS: modInvS,

		base:  base,
		baseS: baseS,

		inModOut: inModOut,

		idx: idx,
	}
}

// EmbedVec returns the embedding of v to the output modulus.
func (emb *VecEmbedder) Embed(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(emb.modOut))
	for i := 0; i < len(emb.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	emb.EmbedTo(vOut, v)
	return vOut
}

// EmbedTo embeds v to vOut.
// If len(vOut) < len(emb.modOut), it only embeds to len(vOut) elements.
func (emb *VecEmbedder) EmbedTo(vOut, v [][]uint64) {
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

		clear(vBool[:])
		for i := 0; i < embedBatch; i += 8 {
			vBool[i+0] = isMixedRadixNegative(vBuf, i+0, emb.modInHalf)
			vBool[i+1] = isMixedRadixNegative(vBuf, i+1, emb.modInHalf)
			vBool[i+2] = isMixedRadixNegative(vBuf, i+2, emb.modInHalf)
			vBool[i+3] = isMixedRadixNegative(vBuf, i+3, emb.modInHalf)

			vBool[i+4] = isMixedRadixNegative(vBuf, i+4, emb.modInHalf)
			vBool[i+5] = isMixedRadixNegative(vBuf, i+5, emb.modInHalf)
			vBool[i+6] = isMixedRadixNegative(vBuf, i+6, emb.modInHalf)
			vBool[i+7] = isMixedRadixNegative(vBuf, i+7, emb.modInHalf)
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

	clear(vBool[:])
	for i := 0; i < len(v[0])-M; i++ {
		vBool[i] = isMixedRadixNegative(vBuf, i, emb.modInHalf)
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
func (emb *VecEmbedder) ModulusIn() []*num.Modulus {
	return emb.modIn
}

// ModulusOut returns the output modulus.
func (emb *VecEmbedder) ModulusOut() []*num.Modulus {
	return emb.modOut
}

// Scaler scales a polynomial from the input modulus to the output modulus.
type Scaler struct {
	opIn  *Operator
	opOut *Operator

	vecSc *VecScaler
	pool  *pool.Pool[*[]uint64]
}

// NewScaler creates a new [Scaler].
func NewScaler(opOut *Operator, opIn *Operator) *Scaler {
	return &Scaler{
		opIn:  opIn,
		opOut: opOut,

		vecSc: NewVecScaler(opOut.mod, opIn.mod),
		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, opIn.params.Rank())
			return &v
		}),
	}
}

// WithPool sets the pool of [Scaler] and returns it.
func (sc *Scaler) WithPool(p *pool.Pool[*[]uint64]) *Scaler {
	sc.pool = p
	return sc
}

// Scale scales e to the output modulus.
func (sc *Scaler) Scale(e *Element, isNTT bool) *Element {
	pOut := NewPoly(e.Rank(), len(sc.opOut.mod))
	sc.ScaleTo(pOut, e, isNTT)
	return pOut
}

// ScaleTo scales e to eOut.
func (sc *Scaler) ScaleTo(eOut, e *Element, isNTT bool) {
	if !(e.Rank() == eOut.Rank() && e.ModLen() == len(sc.opIn.mod) && eOut.ModLen() == len(sc.opOut.mod)) {
		panic("input(s) not consistent")
	}

	if e.IsNTT {
		if isNTT {
			sc.scaleNTTTo(eOut, e)
			eOut.IsNTT = true
		} else {
			eInvNTT, put := getElementFromPool(sc.pool, e.ModLen(), false)
			defer put()

			sc.opIn.InvNTTTo(eInvNTT, e)
			sc.vecSc.ScaleTo(eOut.Coeffs, eInvNTT.Coeffs)
			eOut.IsNTT = false
		}
	} else {
		sc.vecSc.ScaleTo(eOut.Coeffs, e.Coeffs)
		eOut.IsNTT = false
		if isNTT {
			sc.opOut.Slice(0, eOut.ModLen()).FwdNTTTo(eOut, eOut)
		}
	}
}

// scaleNTTTo scales an element in the NTT form.
func (sc *Scaler) scaleNTTTo(eOutNTT, eNTT *Element) {
	M := (len(eNTT.Coeffs[0]) >> 3) << 3

	inLen, gcdLen, outLen := len(eNTT.Coeffs), sc.vecSc.modGCDLen, len(eOutNTT.Coeffs)
	if inLen != len(sc.vecSc.modIn) || outLen != len(sc.vecSc.modOut) {
		panic("input(s) not consistent")
	}

	if len(sc.vecSc.modIn) == sc.vecSc.modGCDLen {
		eNTTCopy, put := getElementFromPool(sc.pool, eNTT.ModLen(), true)
		defer put()

		for i := 0; i < inLen; i++ {
			ii := sc.vecSc.idxInToGCD[i]
			vec.SMulScalarTo(eNTTCopy.Coeffs[ii], eNTT.Coeffs[i], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], sc.vecSc.modIn[i])
		}

		for i := 0; i < outLen; i++ {
			ii := sc.vecSc.idxOutToGCD[i]
			if ii >= 0 {
				copy(eOutNTT.Coeffs[i], eNTTCopy.Coeffs[ii])
			} else {
				clear(eOutNTT.Coeffs[i])
			}
		}

		return
	}

	eNTTCopy, put := getElementFromPool(sc.pool, inLen-gcdLen, false)
	defer put()

	vMul, put := getElementFromPool(sc.pool, gcdLen, false)
	defer put()

	vBoolPtr := sc.pool.Get()
	vBool := *vBoolPtr
	defer sc.pool.Put(vBoolPtr)

	vCorrPtr := sc.pool.Get()
	vCorr := *vCorrPtr
	defer sc.pool.Put(vCorrPtr)

	var qLast uint64
	for i := inLen - 1; i >= 0; i-- {
		if sc.vecSc.idxInToComp[i] >= 0 {
			qLast = sc.vecSc.modIn[i].Value()
			break
		}
	}
	qLastHalf := qLast >> 1

	for i := 0; i < inLen; i++ {
		modIn := sc.vecSc.modIn[i]
		ii := sc.vecSc.idxInToGCD[i]
		if ii >= 0 {
			vec.SMulScalarTo(vMul.Coeffs[ii][:], eNTT.Coeffs[i][:], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], modIn)
		}
	}

	for i := 0; i < inLen; i++ {
		modIn := sc.vecSc.modIn[i]
		ii := sc.vecSc.idxInToComp[i]
		if ii >= 0 {
			vec.SMulScalarTo(eNTTCopy.Coeffs[ii][:], eNTT.Coeffs[i][:], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], modIn)
			if sc.opIn.ntt[i] != nil {
				sc.opIn.ntt[i].InverseTo(eNTTCopy.Coeffs[ii], eNTTCopy.Coeffs[ii])
			}
		}
	}

	if inLen-gcdLen == 1 {
		for i := 0; i < outLen; i++ {
			modOut := sc.vecSc.modOut[i]

			for j := 0; j < M; j += 8 {
				eOutNTT.Coeffs[i][j+0] = embedToModOut(eNTTCopy.Coeffs[0][j+0], modOut, qLast, qLastHalf)
				eOutNTT.Coeffs[i][j+1] = embedToModOut(eNTTCopy.Coeffs[0][j+1], modOut, qLast, qLastHalf)
				eOutNTT.Coeffs[i][j+2] = embedToModOut(eNTTCopy.Coeffs[0][j+2], modOut, qLast, qLastHalf)
				eOutNTT.Coeffs[i][j+3] = embedToModOut(eNTTCopy.Coeffs[0][j+3], modOut, qLast, qLastHalf)

				eOutNTT.Coeffs[i][j+4] = embedToModOut(eNTTCopy.Coeffs[0][j+4], modOut, qLast, qLastHalf)
				eOutNTT.Coeffs[i][j+5] = embedToModOut(eNTTCopy.Coeffs[0][j+5], modOut, qLast, qLastHalf)
				eOutNTT.Coeffs[i][j+6] = embedToModOut(eNTTCopy.Coeffs[0][j+6], modOut, qLast, qLastHalf)
				eOutNTT.Coeffs[i][j+7] = embedToModOut(eNTTCopy.Coeffs[0][j+7], modOut, qLast, qLastHalf)
			}

			for j := M; j < len(eOutNTT.Coeffs[i]); j++ {
				eOutNTT.Coeffs[i][j] = embedToModOut(eNTTCopy.Coeffs[0][j], modOut, qLast, qLastHalf)
			}
		}
	} else {
		for i := 0; i < inLen-gcdLen; i++ {
			for j := i + 1; j < inLen-gcdLen; j++ {
				modInComp := sc.vecSc.modInComp[j]
				vec.SubTo(eNTTCopy.Coeffs[j], eNTTCopy.Coeffs[j], eNTTCopy.Coeffs[i], modInComp)
				vec.SMulScalarTo(eNTTCopy.Coeffs[j], eNTTCopy.Coeffs[j], sc.vecSc.modInCompInv[i][j-i-1], sc.vecSc.modInCompInvS[i][j-i-1], modInComp)
			}
		}

		clear(vBool[:])
		for i := 0; i < M; i += 8 {
			vBool[i+0] = isMixedRadixNegative(eNTTCopy.Coeffs, i+0, sc.vecSc.modInCompHalf)
			vBool[i+1] = isMixedRadixNegative(eNTTCopy.Coeffs, i+1, sc.vecSc.modInCompHalf)
			vBool[i+2] = isMixedRadixNegative(eNTTCopy.Coeffs, i+2, sc.vecSc.modInCompHalf)
			vBool[i+3] = isMixedRadixNegative(eNTTCopy.Coeffs, i+3, sc.vecSc.modInCompHalf)

			vBool[i+4] = isMixedRadixNegative(eNTTCopy.Coeffs, i+4, sc.vecSc.modInCompHalf)
			vBool[i+5] = isMixedRadixNegative(eNTTCopy.Coeffs, i+5, sc.vecSc.modInCompHalf)
			vBool[i+6] = isMixedRadixNegative(eNTTCopy.Coeffs, i+6, sc.vecSc.modInCompHalf)
			vBool[i+7] = isMixedRadixNegative(eNTTCopy.Coeffs, i+7, sc.vecSc.modInCompHalf)
		}
		for i := M; i < len(eNTTCopy.Coeffs[0]); i++ {
			vBool[i] = isMixedRadixNegative(eNTTCopy.Coeffs, i, sc.vecSc.modInCompHalf)
		}

		for i := 0; i < outLen; i++ {
			base, baseS := sc.vecSc.base[i], sc.vecSc.baseS[i]
			inCompModOut := sc.vecSc.inCompModOut[i]
			modOut := sc.vecSc.modOut[i]

			vec.SMulScalarTo(eOutNTT.Coeffs[i], eNTTCopy.Coeffs[0], base[0], baseS[0], modOut)
			for j := 1; j < inLen-gcdLen; j++ {
				vec.SMulAddScalarTo(eOutNTT.Coeffs[i], eNTTCopy.Coeffs[j], base[j], baseS[j], modOut)
			}
			vec.MulScalarTo(vCorr[:], vBool[:], inCompModOut, nil)
			vec.SubTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i], vCorr[:], modOut)
		}
	}

	for i := 0; i < outLen; i++ {
		if sc.opOut.ntt[i] != nil {
			sc.opOut.ntt[i].ForwardTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i])
		}

		modOut := sc.vecSc.modOut[i]

		ii := sc.vecSc.idxOutToGCD[i]
		if ii >= 0 {
			vec.SubTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i], vMul.Coeffs[ii], modOut)
		}
		vec.SMulScalarTo(eOutNTT.Coeffs[i], eOutNTT.Coeffs[i], sc.vecSc.negInCompInvModOut[i], sc.vecSc.negInCompInvModOutS[i], modOut)
	}

	eOutNTT.IsNTT = true
}

// OperatorIn returns the input operator.
func (emb *Scaler) OperatorIn() *Operator {
	return emb.opIn
}

// OperatorOut returns the output operator.
func (emb *Scaler) OperatorOut() *Operator {
	return emb.opOut
}

// VecScaler scales a polynomial to different modulus.
// In other words, it computes
//
//	[p]_modIn -> [(modOut / modIn) * p]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
type VecScaler struct {
	// modIn is the input modulus.
	modIn []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus
	// modInComp equals modIn/modGCD.
	modInComp []*num.Modulus

	// modInCompHalf is half of modInComp.
	modInCompHalf []uint64

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
	emb *VecEmbedder

	// outCompModInComp equals modOut/modGCD modulo modIn.
	outCompModIn []uint64
	// outCompModInS is the Shoup form of outCompModInComp.
	outCompModInS []uint64
	// negInCompInvModOut equals the inverse of modIn/modGCD modulo modOut.
	negInCompInvModOut []uint64
	// negInCompInvModOutS is the Shoup form of negInCompInvModOut.
	negInCompInvModOutS []uint64
}

func NewVecScaler(modOut, modIn []*num.Modulus) *VecScaler {
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
	modInCompHalf := make([]uint64, 0, len(modIn)-len(modGCD))
	for i := 0; i < len(modIn); i++ {
		if idxInToGCD[i] == -1 {
			modInComp = append(modInComp, modIn[i])
			modInCompHalf = append(modInCompHalf, modIn[i].Value()>>1)
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

	return &VecScaler{
		modIn:     modIn,
		modOut:    modOut,
		modInComp: modInComp,

		modInCompHalf: modInCompHalf,

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

// ScaleVec returns the scaled vector of v.
func (sc *VecScaler) Scale(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(sc.modOut))
	for i := 0; i < len(sc.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	sc.ScaleTo(vOut, v)
	return vOut
}

// ScaleVecTo scales v to vOut.
func (sc *VecScaler) ScaleTo(vOut, v [][]uint64) {
	M := (len(v[0]) >> logEmbedBatch) << logEmbedBatch
	L := unsafe.Sizeof(uint64(0))

	inLen, gcdLen, outLen := len(sc.modIn), sc.modGCDLen, len(sc.modOut)
	if len(v) != len(sc.modIn) || len(vOut) != len(sc.modOut) {
		panic("input(s) not consistent")
	}

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
				vec.SMulScalarTo(vBuf[ii][:], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modIn[i])
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

		for i := 0; i < inLen; i++ {
			w := v[i][M:]

			ii := sc.idxInToGCD[i]
			vec.SMulScalarTo(vBuf[ii][:len(v[0])-M], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modIn[i])
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

			clear(vBool[:])
			for i := 0; i < embedBatch; i += 8 {
				vBool[i+0] = isMixedRadixNegative(vBuf, i+0, sc.modInCompHalf)
				vBool[i+1] = isMixedRadixNegative(vBuf, i+1, sc.modInCompHalf)
				vBool[i+2] = isMixedRadixNegative(vBuf, i+2, sc.modInCompHalf)
				vBool[i+3] = isMixedRadixNegative(vBuf, i+3, sc.modInCompHalf)

				vBool[i+4] = isMixedRadixNegative(vBuf, i+4, sc.modInCompHalf)
				vBool[i+5] = isMixedRadixNegative(vBuf, i+5, sc.modInCompHalf)
				vBool[i+6] = isMixedRadixNegative(vBuf, i+6, sc.modInCompHalf)
				vBool[i+7] = isMixedRadixNegative(vBuf, i+7, sc.modInCompHalf)
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

		clear(vBool[:])
		for i := 0; i < len(v[0])-M; i++ {
			vBool[i] = isMixedRadixNegative(vBuf, i, sc.modInCompHalf)
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
func (sc *VecScaler) ModulusIn() []*num.Modulus {
	return sc.modIn
}

// ModulusOut returns the output modulus.
func (sc *VecScaler) ModulusOut() []*num.Modulus {
	return sc.modOut
}
