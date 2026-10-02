package crt

import (
	"slices"
	"unsafe"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
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

// Embedder embeds *[Scalar] and *[Poly] into different modulus.
// In other words, it computes
//
//	[e]_modIn -> [e]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
type Embedder struct {
	OpIn  *Operator
	OpOut *Operator

	vecEmb *VecEmbedder
	pool   *pool.Pool[*[]uint64]
}

// NewEmbedder creates a new [Embedder].
func NewEmbedder(opOut *Operator, opIn *Operator) *Embedder {
	if !opOut.params.IsEqual(opIn.params) {
		panic("inconsistent parameters")
	}

	return &Embedder{
		OpIn:  opIn,
		OpOut: opOut,

		vecEmb: NewVecEmbedder(opOut.mod, opIn.mod),
		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, max(opIn.params.Rank(), len(opIn.mod)))
			return &v
		}),
	}
}

// WithPool sets the pool of [Embedder] and returns it.
func (emb *Embedder) WithPool(p *pool.Pool[*[]uint64]) *Embedder {
	emb.pool = p
	return emb
}

// Embed returns the embedding of e to the output modulus.
func (emb *Embedder) Embed[T *Scalar | *Poly](e T, outForm dft.Form) T {
	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	if eIsScalar {
		eOut := NewScalar(len(emb.OpOut.mod))
		emb.EmbedTo(eOut, eScalar, outForm)
		return any(eOut).(T)
	}
	eOut := emb.OpOut.NewPoly()
	emb.EmbedTo(eOut, ePoly, outForm)
	return any(eOut).(T)
}

// EmbedTo embeds e to eOut.
func (emb *Embedder) EmbedTo[T *Scalar | *Poly](eOut, e T, outForm dft.Form) {
	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	if eIsScalar {
		emb.embedScalarTo(eOutScalar, eScalar)
		return
	}

	switch {
	case ePoly.Form == dft.FormNTT && outForm == dft.FormNTT:
		emb.embedNTTTo(eOutPoly, ePoly)
		eOutPoly.Form = dft.FormNTT

	case ePoly.Form == dft.FormNTT && outForm == dft.FormCoeff:
		pInvNTT, put := getPolyFromPool(emb.pool, ePoly.Rank(), ePoly.ModLen(), dft.FormCoeff)
		defer put()

		emb.OpIn.InvNTTTo(pInvNTT, ePoly)
		emb.vecEmb.EmbedTo(eOutPoly.Coeffs, pInvNTT.Coeffs)
		eOutPoly.Form = dft.FormCoeff

	case ePoly.Form == dft.FormCoeff:
		emb.vecEmb.EmbedTo(eOutPoly.Coeffs, ePoly.Coeffs)
		eOutPoly.Form = dft.FormCoeff
		if outForm == dft.FormNTT {
			emb.OpOut.Slice(0, eOutPoly.ModLen()).FwdNTTTo(eOutPoly, eOutPoly)
		}
	}
}

func (emb *Embedder) embedNTTTo(pOutNTT, pNTT *Poly) {
	M := (len(pNTT.Coeffs[0]) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	inLen, outLen := len(pNTT.Coeffs), len(pOutNTT.Coeffs)
	if inLen != len(emb.vecEmb.modIn) || outLen > len(emb.vecEmb.modOut) {
		panic("input(s) not consistent")
	}

	N := len(pNTT.Coeffs[0])
	for i := 1; i < inLen; i++ {
		if len(pNTT.Coeffs[i]) != N {
			panic("inconsistent input(s)")
		}
	}
	for i := 0; i < outLen; i++ {
		if len(pOutNTT.Coeffs[i]) != N {
			panic("inconsistent input(s)")
		}
	}

	pNTTCopy, put := getPolyFromPool(emb.pool, N, inLen, dft.FormCoeff)
	defer put()
	pNTTCopy.CopyFrom(pNTT)

	pInvNTT, put := getPolyFromPool(emb.pool, N, inLen, dft.FormCoeff)
	defer put()
	emb.OpIn.InvNTTTo(pInvNTT, pNTT)

	if inLen == 1 {
		qv := emb.vecEmb.modInSorted[0].Value()
		halfQv := qv >> 1

		r := unsafe.Pointer(unsafe.SliceData(pInvNTT.Coeffs[0]))

		for i := 0; i < M; i += 8 {
			wIn := (*[8]uint64)(unsafe.Add(r, uintptr(i)*L))

			for j := 0; j < outLen; j++ {
				rOut := unsafe.Pointer(unsafe.SliceData(pOutNTT.Coeffs[j]))
				wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(i)*L))

				qOut := emb.vecEmb.modOut[j].Value()
				qOutDiv, _ := emb.vecEmb.modOut[j].Div()

				if emb.vecEmb.idx[j] != 0 {
					wOut[0] = embedToModOut(wIn[0], qOut, qOutDiv, qv, halfQv)
					wOut[1] = embedToModOut(wIn[1], qOut, qOutDiv, qv, halfQv)
					wOut[2] = embedToModOut(wIn[2], qOut, qOutDiv, qv, halfQv)
					wOut[3] = embedToModOut(wIn[3], qOut, qOutDiv, qv, halfQv)

					wOut[4] = embedToModOut(wIn[4], qOut, qOutDiv, qv, halfQv)
					wOut[5] = embedToModOut(wIn[5], qOut, qOutDiv, qv, halfQv)
					wOut[6] = embedToModOut(wIn[6], qOut, qOutDiv, qv, halfQv)
					wOut[7] = embedToModOut(wIn[7], qOut, qOutDiv, qv, halfQv)
				}
			}
		}

		for i := M; i < len(pInvNTT.Coeffs[0]); i++ {
			c := pInvNTT.Coeffs[0][i]
			for j := 0; j < outLen; j++ {
				if emb.vecEmb.idx[j] != 0 {
					qOut := emb.vecEmb.modOut[j].Value()
					qOutDiv, _ := emb.vecEmb.modOut[j].Div()

					pOutNTT.Coeffs[j][i] = embedToModOut(c, qOut, qOutDiv, qv, halfQv)
				}
			}
		}

		for i := 0; i < outLen; i++ {
			if emb.vecEmb.idx[i] == 0 {
				copy(pOutNTT.Coeffs[i][:], pNTTCopy.Coeffs[0][:])
			} else if emb.OpOut.ntt[i] != nil {
				emb.OpOut.ntt[i].ForwardTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i])
			}
		}

		return
	}

	vBoolPtr := emb.pool.Get()
	vBool := (*vBoolPtr)[:N]
	defer emb.pool.Put(vBoolPtr)

	vCorrPtr := emb.pool.Get()
	vCorr := (*vCorrPtr)[:N]
	defer emb.pool.Put(vCorrPtr)

	pInvNTTCoeffs := make([][]uint64, inLen)
	for i := 0; i < inLen; i++ {
		pInvNTTCoeffs[i] = pInvNTT.Coeffs[emb.vecEmb.modInMap[i]]
	}

	for i := 0; i < inLen; i++ {
		for j := i + 1; j < inLen; j++ {
			vec.SubTo(pInvNTTCoeffs[j], pInvNTTCoeffs[j], pInvNTTCoeffs[i], emb.vecEmb.modInSorted[j])
			vec.SMulTo(pInvNTTCoeffs[j], pInvNTTCoeffs[j], emb.vecEmb.modInv[i][j-i-1], emb.vecEmb.modInvS[i][j-i-1], emb.vecEmb.modInSorted[j])
		}
	}

	clear(vBool[:])
	for i := 0; i < M; i += 8 {
		vBool[i+0] = isMixedRadixNegative(pInvNTTCoeffs, i+0, emb.vecEmb.modInHalf)
		vBool[i+1] = isMixedRadixNegative(pInvNTTCoeffs, i+1, emb.vecEmb.modInHalf)
		vBool[i+2] = isMixedRadixNegative(pInvNTTCoeffs, i+2, emb.vecEmb.modInHalf)
		vBool[i+3] = isMixedRadixNegative(pInvNTTCoeffs, i+3, emb.vecEmb.modInHalf)

		vBool[i+4] = isMixedRadixNegative(pInvNTTCoeffs, i+4, emb.vecEmb.modInHalf)
		vBool[i+5] = isMixedRadixNegative(pInvNTTCoeffs, i+5, emb.vecEmb.modInHalf)
		vBool[i+6] = isMixedRadixNegative(pInvNTTCoeffs, i+6, emb.vecEmb.modInHalf)
		vBool[i+7] = isMixedRadixNegative(pInvNTTCoeffs, i+7, emb.vecEmb.modInHalf)
	}
	for i := M; i < len(pNTT.Coeffs[0]); i++ {
		vBool[i] = isMixedRadixNegative(pInvNTTCoeffs, i, emb.vecEmb.modInHalf)
	}

	for i := 0; i < outLen; i++ {
		if 0 <= emb.vecEmb.idx[i] && emb.vecEmb.idx[i] < inLen {
			continue
		}

		base, baseS := emb.vecEmb.base[i], emb.vecEmb.baseS[i]
		inModOut := emb.vecEmb.inModOut[i]
		modOut := emb.vecEmb.modOut[i]

		vec.SMulTo(pOutNTT.Coeffs[i], pInvNTTCoeffs[0], base[0], baseS[0], modOut)
		for j := 1; j < inLen; j++ {
			vec.SMulAddTo(pOutNTT.Coeffs[i], pInvNTTCoeffs[j], base[j], baseS[j], modOut)
		}
		vec.MulTo(vCorr, vBool, inModOut, nil)
		vec.SubTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i], vCorr, modOut)
	}

	for i := 0; i < outLen; i++ {
		if 0 <= emb.vecEmb.idx[i] && emb.vecEmb.idx[i] < inLen {
			copy(pOutNTT.Coeffs[i], pNTTCopy.Coeffs[emb.vecEmb.modInMap[emb.vecEmb.idx[i]]])
		} else if emb.OpOut.ntt[i] != nil {
			emb.OpOut.ntt[i].ForwardTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i])
		}
	}

	pOutNTT.Form = dft.FormNTT
}

func (emb *Embedder) embedScalarTo(cOut, c *Scalar) {
	vecEmb := emb.vecEmb
	inLen, outLen := c.ModLen(), cOut.ModLen()
	if inLen == 0 || inLen != len(vecEmb.modIn) || outLen > len(vecEmb.modOut) {
		panic("inconsistent input(s)")
	}

	if inLen == 1 {
		x := c.Value[0]
		qv := vecEmb.modInSorted[0].Value()
		halfQv := qv >> 1
		for i := 0; i < outLen; i++ {
			if vecEmb.idx[i] == 0 {
				cOut.Value[i] = x
			} else {
				qOut := vecEmb.modOut[i].Value()
				qOutDivHi, _ := vecEmb.modOut[i].Div()
				cOut.Value[i] = embedToModOut(x, qOut, qOutDivHi, qv, halfQv)
			}
		}
		return
	}

	cBufPtr := emb.pool.Get()
	cBuf := (*cBufPtr)[:inLen]
	defer emb.pool.Put(cBufPtr)

	for i := 0; i < inLen; i++ {
		cBuf[i] = c.Value[vecEmb.modInMap[i]]
	}
	for i := 0; i < outLen; i++ {
		if 0 <= vecEmb.idx[i] && vecEmb.idx[i] < inLen {
			cOut.Value[i] = cBuf[vecEmb.idx[i]]
		}
	}

	for i := 0; i < inLen; i++ {
		for j := i + 1; j < inLen; j++ {
			cBuf[j] = num.Sub(cBuf[j], cBuf[i], vecEmb.modInSorted[j])
			cBuf[j] = num.SMul(cBuf[j], vecEmb.modInv[i][j-i-1], vecEmb.modInvS[i][j-i-1], vecEmb.modInSorted[j])
		}
	}

	isNegative := false
	for i := inLen - 1; i >= 0; i-- {
		if cBuf[i] > vecEmb.modInHalf[i] {
			isNegative = true
			break
		}
		if cBuf[i] < vecEmb.modInHalf[i] {
			break
		}
	}

	for i := 0; i < outLen; i++ {
		if 0 <= vecEmb.idx[i] && vecEmb.idx[i] < inLen {
			continue
		}

		base, baseS := vecEmb.base[i], vecEmb.baseS[i]
		modOut := vecEmb.modOut[i]
		x := num.SMul(cBuf[0], base[0], baseS[0], modOut)
		for j := 1; j < inLen; j++ {
			x = num.Add(x, num.SMul(cBuf[j], base[j], baseS[j], modOut), modOut)
		}
		if isNegative {
			x = num.Sub(x, vecEmb.inModOut[i], modOut)
		}
		cOut.Value[i] = x
	}
}

// VecEmbedder embeds a vector into different modulus.
// In other words, it computes
//
//	[v]_modIn -> [v]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
type VecEmbedder struct {
	// modIn is the input modulus in the caller's order.
	modIn []*num.Modulus
	// modInSorted is the input modulus sorted for mixed-radix conversion.
	modInSorted []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus

	// modInMap maps the index of modInSorted to modIn.
	modInMap []int

	// modInHalf is modInSorted/2 in the mixed-radix basis of modInSorted.
	modInHalf []uint64

	// modInv is the inverse of modInSorted.
	modInv [][]uint64
	// modInvS is the Shoup form of modInv.
	modInvS [][]uint64

	// base is the mixed-radix basis of modInSorted.
	base [][]uint64
	// baseS is the Shoup form of base.
	baseS [][]uint64

	// inModOut equals modInSorted modulo modOut.
	inModOut []uint64

	// idx holds the index of the sorted input modulus limb if it overlaps with the output modulus limb.
	// For example, if modOut[i] = modInSorted[j], then idx[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	idx []int
}

// NewVecEmbedder creates a new [VecEmbedder].
func NewVecEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *VecEmbedder {
	modInSorted := make([]*num.Modulus, len(modIn))
	copy(modInSorted, modIn)
	slices.SortFunc(modInSorted, num.CmpModulus)

	modInMap := make([]int, len(modInSorted))
	for i := range modInSorted {
		for j := range modIn {
			if modInSorted[i].Value() == modIn[j].Value() {
				modInMap[i] = j
				break
			}
		}
	}

	modInHalf := halfProductMixedRadix(modInSorted)

	modInv := make([][]uint64, len(modInSorted))
	modInvS := make([][]uint64, len(modInSorted))
	for i := 0; i < len(modInSorted); i++ {
		modInv[i] = make([]uint64, len(modInSorted)-i-1)
		modInvS[i] = make([]uint64, len(modInSorted)-i-1)
		for j := 0; j < len(modInSorted)-i-1; j++ {
			modInv[i][j] = num.Inv(modInSorted[i].Value(), modInSorted[i+j+1])
			modInvS[i][j] = num.SForm(modInv[i][j], modInSorted[i+j+1])
		}
	}

	base := make([][]uint64, len(modOut))
	baseS := make([][]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		base[i] = make([]uint64, len(modInSorted))
		base[i][0] = 1
		for j := 1; j < len(modInSorted); j++ {
			base[i][j] = num.Mul(base[i][j-1], modInSorted[j-1].Value(), modOut[i])
		}
		baseS[i] = vec.SForm(base[i], modOut[i])
	}

	inModOut := make([]uint64, len(modOut))
	idx := make([]int, len(modOut))
	for i := 0; i < len(modOut); i++ {
		inModOut[i] = 1
		idx[i] = -1
		for j := 0; j < len(modInSorted); j++ {
			inModOut[i] = num.Mul(inModOut[i], modInSorted[j].Value(), modOut[i])
			if modOut[i].Value() == modInSorted[j].Value() {
				idx[i] = j
			}
		}
	}

	return &VecEmbedder{
		modIn:       modIn,
		modInSorted: modInSorted,
		modOut:      modOut,

		modInMap: modInMap,

		modInHalf: modInHalf,

		modInv:  modInv,
		modInvS: modInvS,

		base:  base,
		baseS: baseS,

		inModOut: inModOut,

		idx: idx,
	}
}

// Embed returns the embedding of v to the output modulus.
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
//
// Panics when v and vOut is inconsistent.
func (emb *VecEmbedder) EmbedTo(vOut, v [][]uint64) {
	inLen, outLen := len(v), len(vOut)
	if inLen == 0 || inLen != len(emb.modIn) || outLen > len(emb.modOut) {
		panic("inconsistent input(s)")
	}

	N := len(v[0])
	for i := 1; i < inLen; i++ {
		if len(v[i]) != N {
			panic("inconsistent input(s)")
		}
	}
	for i := 0; i < outLen; i++ {
		if len(vOut[i]) != N {
			panic("inconsistent input(s)")
		}
	}

	M := (N >> logEmbedBatch) << logEmbedBatch
	L := unsafe.Sizeof(uint64(0))

	if inLen == 1 {
		qv := emb.modInSorted[0].Value()
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

				qOut := emb.modOut[j].Value()
				qOutDivHi, _ := emb.modOut[j].Div()

				if emb.idx[j] == 0 {
					copy(wOut[:], vBuf[:])
				} else {
					wOut[0] = embedToModOut(vBuf[0], qOut, qOutDivHi, qv, halfQv)
					wOut[1] = embedToModOut(vBuf[1], qOut, qOutDivHi, qv, halfQv)
					wOut[2] = embedToModOut(vBuf[2], qOut, qOutDivHi, qv, halfQv)
					wOut[3] = embedToModOut(vBuf[3], qOut, qOutDivHi, qv, halfQv)

					wOut[4] = embedToModOut(vBuf[4], qOut, qOutDivHi, qv, halfQv)
					wOut[5] = embedToModOut(vBuf[5], qOut, qOutDivHi, qv, halfQv)
					wOut[6] = embedToModOut(vBuf[6], qOut, qOutDivHi, qv, halfQv)
					wOut[7] = embedToModOut(vBuf[7], qOut, qOutDivHi, qv, halfQv)
				}
			}
		}

		for i := M; i < len(v[0]); i++ {
			c := v[0][i]
			for j := 0; j < outLen; j++ {
				if emb.idx[j] == 0 {
					vOut[j][i] = c
				} else {
					qOut := emb.modOut[j].Value()
					qOutDivHi, _ := emb.modOut[j].Div()
					vOut[j][i] = embedToModOut(c, qOut, qOutDivHi, qv, halfQv)
				}
			}
		}

		return
	}

	vBuf := make([]*[embedBatch]uint64, inLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get()
	}
	defer func() {
		for i := range vBuf {
			embed64Pool.Put(vBuf[i])
		}
	}()

	vBool := embed64Pool.Get()
	defer embed64Pool.Put(vBool)
	vCorr := embed64Pool.Get()
	defer embed64Pool.Put(vCorr)

	for k := 0; k < M; k += embedBatch {
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[emb.modInMap[i]]))
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
				vec.SubTo(vBuf[j][:], vBuf[j][:], wBuf[:], emb.modInSorted[j])
				vec.SMulTo(vBuf[j][:], vBuf[j][:], emb.modInv[i][j-i-1], emb.modInvS[i][j-i-1], emb.modInSorted[j])
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

			vec.SMulTo(wOut[:], vBuf[0][:], base[0], baseS[0], modOut)
			for j := 1; j < inLen; j++ {
				vec.SMulAddTo(wOut[:], vBuf[j][:], base[j], baseS[j], modOut)
			}
			vec.MulTo(vCorr[:], vBool[:], inModOut, nil)
			vec.SubTo(wOut[:], wOut[:], vCorr[:], modOut)
		}
	}

	for i := 0; i < inLen; i++ {
		copy(vBuf[i][:len(v[0])-M], v[emb.modInMap[i]][M:])
	}

	for i := 0; i < outLen; i++ {
		if 0 <= emb.idx[i] && emb.idx[i] < inLen {
			copy(vOut[i][M:], vBuf[emb.idx[i]][:len(v[0])-M])
		}
	}

	for i := 0; i < inLen; i++ {
		wBuf := vBuf[i][:len(v[0])-M]
		for j := i + 1; j < inLen; j++ {
			vec.SubTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], wBuf[:], emb.modInSorted[j])
			vec.SMulTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], emb.modInv[i][j-i-1], emb.modInvS[i][j-i-1], emb.modInSorted[j])
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

		vec.SMulTo(wOut[:], vBuf[0][:len(v[0])-M], base[0], baseS[0], modOut)
		for j := 1; j < inLen; j++ {
			vec.SMulAddTo(wOut[:], vBuf[j][:len(v[0])-M], base[j], baseS[j], modOut)
		}
		vec.MulTo(vCorr[:len(v[0])-M], vBool[:len(v[0])-M], inModOut, nil)
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

// Scaler scales *[Scalar] and *[Poly] from the input modulus to the output modulus.
// In other words, it computes
//
//	[e]_modIn -> [(modOut / modIn) * e]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
// Each modulus from modIn and modOut should be equal or coprime.
type Scaler struct {
	OpIn  *Operator
	OpOut *Operator

	vecSc *VecScaler
	pool  *pool.Pool[*[]uint64]
}

// NewScaler creates a new [Scaler].
func NewScaler(opOut *Operator, opIn *Operator) *Scaler {
	if !opOut.params.IsEqual(opIn.params) {
		panic("inconsistent parameters")
	}

	return &Scaler{
		OpIn:  opIn,
		OpOut: opOut,

		vecSc: NewVecScaler(opOut.mod, opIn.mod),
		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, max(opIn.params.Rank(), len(opIn.mod)))
			return &v
		}),
	}
}

// WithPool sets the pool of [Scaler] and returns it.
func (sc *Scaler) WithPool(p *pool.Pool[*[]uint64]) *Scaler {
	sc.pool = p
	return sc
}

// Scale returns the scaling of e to the output modulus.
func (sc *Scaler) Scale[T *Scalar | *Poly](e T, outForm dft.Form) T {
	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	if eIsScalar {
		eOut := NewScalar(len(sc.OpOut.mod))
		sc.ScaleTo(eOut, eScalar, outForm)
		return any(eOut).(T)
	}
	eOut := sc.OpOut.NewPoly()
	sc.ScaleTo(eOut, ePoly, outForm)
	return any(eOut).(T)
}

// ScaleTo scales e to eOut.
func (sc *Scaler) ScaleTo[T *Scalar | *Poly](eOut, e T, outForm dft.Form) {
	eOutScalar, _ := any(eOut).(*Scalar)
	eOutPoly, _ := any(eOut).(*Poly)

	eScalar, eIsScalar := any(e).(*Scalar)
	ePoly, _ := any(e).(*Poly)

	if eIsScalar {
		sc.scaleScalarTo(eOutScalar, eScalar)
		return
	}

	switch {
	case ePoly.Form == dft.FormNTT && outForm == dft.FormNTT:
		sc.scaleNTTTo(eOutPoly, ePoly)
		eOutPoly.Form = dft.FormNTT

	case ePoly.Form == dft.FormNTT && outForm == dft.FormCoeff:
		pInvNTT, put := getPolyFromPool(sc.pool, ePoly.Rank(), ePoly.ModLen(), dft.FormCoeff)
		defer put()

		sc.OpIn.InvNTTTo(pInvNTT, ePoly)
		sc.vecSc.ScaleTo(eOutPoly.Coeffs, pInvNTT.Coeffs)
		eOutPoly.Form = dft.FormCoeff

	case ePoly.Form == dft.FormCoeff:
		sc.vecSc.ScaleTo(eOutPoly.Coeffs, ePoly.Coeffs)
		eOutPoly.Form = dft.FormCoeff
		if outForm == dft.FormNTT {
			sc.OpOut.Slice(0, eOutPoly.ModLen()).FwdNTTTo(eOutPoly, eOutPoly)
		}
	}
}

func (sc *Scaler) scaleNTTTo(pOutNTT, pNTT *Poly) {
	M := (len(pNTT.Coeffs[0]) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	inLen, gcdLen, outLen := len(pNTT.Coeffs), sc.vecSc.modGCDLen, len(pOutNTT.Coeffs)
	if inLen != len(sc.vecSc.modIn) || outLen != len(sc.vecSc.modOut) {
		panic("input(s) not consistent")
	}

	N := len(pNTT.Coeffs[0])
	for i := 1; i < inLen; i++ {
		if len(pNTT.Coeffs[i]) != N {
			panic("inconsistent input(s)")
		}
	}
	for i := 0; i < outLen; i++ {
		if len(pOutNTT.Coeffs[i]) != N {
			panic("inconsistent input(s)")
		}
	}

	if inLen == outLen && inLen == gcdLen {
		pNTTCopy, put := getPolyFromPool(sc.pool, N, pNTT.ModLen(), dft.FormNTT)
		defer put()

		for j := 0; j < outLen; j++ {
			copy(pNTTCopy.Coeffs[j], pNTT.Coeffs[sc.vecSc.idxGCDToIn[sc.vecSc.idxOutToGCD[j]]])
		}
		pOutNTT.CopyFrom(pNTTCopy)
		pOutNTT.Form = dft.FormNTT
		return
	}

	if len(sc.vecSc.modInSorted) == sc.vecSc.modGCDLen {
		eNTTCopy, put := getPolyFromPool(sc.pool, N, pNTT.ModLen(), dft.FormNTT)
		defer put()

		for i := 0; i < inLen; i++ {
			ii := sc.vecSc.idxInToGCD[i]
			vec.SMulTo(eNTTCopy.Coeffs[ii], pNTT.Coeffs[sc.vecSc.modInMap[i]], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], sc.vecSc.modInSorted[i])
		}

		for i := 0; i < outLen; i++ {
			ii := sc.vecSc.idxOutToGCD[i]
			if ii >= 0 {
				copy(pOutNTT.Coeffs[i], eNTTCopy.Coeffs[ii])
			} else {
				clear(pOutNTT.Coeffs[i])
			}
		}

		return
	}

	pNTTCopy, put := getPolyFromPool(sc.pool, N, inLen-gcdLen, dft.FormCoeff)
	defer put()

	pMul, put := getPolyFromPool(sc.pool, N, gcdLen, dft.FormCoeff)
	defer put()

	vBoolPtr := sc.pool.Get()
	vBool := (*vBoolPtr)[:N]
	defer sc.pool.Put(vBoolPtr)

	vCorrPtr := sc.pool.Get()
	vCorr := (*vCorrPtr)[:N]
	defer sc.pool.Put(vCorrPtr)

	var qLast uint64
	for i := inLen - 1; i >= 0; i-- {
		if sc.vecSc.idxInToComp[i] >= 0 {
			qLast = sc.vecSc.modInSorted[i].Value()
			break
		}
	}
	qLastHalf := qLast >> 1

	for i := 0; i < inLen; i++ {
		modIn := sc.vecSc.modInSorted[i]
		ii := sc.vecSc.idxInToGCD[i]
		if ii >= 0 {
			vec.SMulTo(pMul.Coeffs[ii][:], pNTT.Coeffs[sc.vecSc.modInMap[i]][:], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], modIn)
		}
	}

	for i := 0; i < inLen; i++ {
		modIn := sc.vecSc.modInSorted[i]
		ii := sc.vecSc.idxInToComp[i]
		if ii >= 0 {
			srcIdx := sc.vecSc.modInMap[i]
			vec.SMulTo(pNTTCopy.Coeffs[ii][:], pNTT.Coeffs[srcIdx][:], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], modIn)
			if sc.OpIn.ntt[srcIdx] != nil {
				sc.OpIn.ntt[srcIdx].InverseTo(pNTTCopy.Coeffs[ii], pNTTCopy.Coeffs[ii])
			}
		}
	}

	if inLen-gcdLen == 1 {
		r := unsafe.Pointer(unsafe.SliceData(pNTTCopy.Coeffs[0]))
		for i := 0; i < outLen; i++ {
			rOut := unsafe.Pointer(unsafe.SliceData(pOutNTT.Coeffs[i]))

			qOut := sc.vecSc.modOut[i].Value()
			qOutDivHi, _ := sc.vecSc.modOut[i].Div()

			for j := 0; j < M; j += 8 {
				wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(j)*L))
				w := (*[8]uint64)(unsafe.Add(r, uintptr(j)*L))

				wOut[0] = embedToModOut(w[0], qOut, qOutDivHi, qLast, qLastHalf)
				wOut[1] = embedToModOut(w[1], qOut, qOutDivHi, qLast, qLastHalf)
				wOut[2] = embedToModOut(w[2], qOut, qOutDivHi, qLast, qLastHalf)
				wOut[3] = embedToModOut(w[3], qOut, qOutDivHi, qLast, qLastHalf)

				wOut[4] = embedToModOut(w[4], qOut, qOutDivHi, qLast, qLastHalf)
				wOut[5] = embedToModOut(w[5], qOut, qOutDivHi, qLast, qLastHalf)
				wOut[6] = embedToModOut(w[6], qOut, qOutDivHi, qLast, qLastHalf)
				wOut[7] = embedToModOut(w[7], qOut, qOutDivHi, qLast, qLastHalf)
			}

			for j := M; j < len(pOutNTT.Coeffs[i]); j++ {
				pOutNTT.Coeffs[i][j] = embedToModOut(pNTTCopy.Coeffs[0][j], qOut, qOutDivHi, qLast, qLastHalf)
			}
		}
	} else {
		for i := 0; i < inLen-gcdLen; i++ {
			for j := i + 1; j < inLen-gcdLen; j++ {
				modInComp := sc.vecSc.modInComp[j]
				vec.SubTo(pNTTCopy.Coeffs[j], pNTTCopy.Coeffs[j], pNTTCopy.Coeffs[i], modInComp)
				vec.SMulTo(pNTTCopy.Coeffs[j], pNTTCopy.Coeffs[j], sc.vecSc.modInCompInv[i][j-i-1], sc.vecSc.modInCompInvS[i][j-i-1], modInComp)
			}
		}

		clear(vBool[:])
		for i := 0; i < M; i += 8 {
			vBool[i+0] = isMixedRadixNegative(pNTTCopy.Coeffs, i+0, sc.vecSc.modInCompHalf)
			vBool[i+1] = isMixedRadixNegative(pNTTCopy.Coeffs, i+1, sc.vecSc.modInCompHalf)
			vBool[i+2] = isMixedRadixNegative(pNTTCopy.Coeffs, i+2, sc.vecSc.modInCompHalf)
			vBool[i+3] = isMixedRadixNegative(pNTTCopy.Coeffs, i+3, sc.vecSc.modInCompHalf)

			vBool[i+4] = isMixedRadixNegative(pNTTCopy.Coeffs, i+4, sc.vecSc.modInCompHalf)
			vBool[i+5] = isMixedRadixNegative(pNTTCopy.Coeffs, i+5, sc.vecSc.modInCompHalf)
			vBool[i+6] = isMixedRadixNegative(pNTTCopy.Coeffs, i+6, sc.vecSc.modInCompHalf)
			vBool[i+7] = isMixedRadixNegative(pNTTCopy.Coeffs, i+7, sc.vecSc.modInCompHalf)
		}
		for i := M; i < len(pNTTCopy.Coeffs[0]); i++ {
			vBool[i] = isMixedRadixNegative(pNTTCopy.Coeffs, i, sc.vecSc.modInCompHalf)
		}

		for i := 0; i < outLen; i++ {
			base, baseS := sc.vecSc.base[i], sc.vecSc.baseS[i]
			inCompModOut := sc.vecSc.inCompModOut[i]
			modOut := sc.vecSc.modOut[i]

			vec.SMulTo(pOutNTT.Coeffs[i], pNTTCopy.Coeffs[0], base[0], baseS[0], modOut)
			for j := 1; j < inLen-gcdLen; j++ {
				vec.SMulAddTo(pOutNTT.Coeffs[i], pNTTCopy.Coeffs[j], base[j], baseS[j], modOut)
			}
			vec.MulTo(vCorr[:], vBool[:], inCompModOut, nil)
			vec.SubTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i], vCorr[:], modOut)
		}
	}

	for i := 0; i < outLen; i++ {
		if sc.OpOut.ntt[i] != nil {
			sc.OpOut.ntt[i].ForwardTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i])
		}

		modOut := sc.vecSc.modOut[i]

		ii := sc.vecSc.idxOutToGCD[i]
		if ii >= 0 {
			vec.SubTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i], pMul.Coeffs[ii], modOut)
		}
		vec.SMulTo(pOutNTT.Coeffs[i], pOutNTT.Coeffs[i], sc.vecSc.negInCompInvModOut[i], sc.vecSc.negInCompInvModOutS[i], modOut)
	}

	pOutNTT.Form = dft.FormNTT
}

func (sc *Scaler) scaleScalarTo(cOut, c *Scalar) {
	inLen, gcdLen, outLen := c.ModLen(), sc.vecSc.modGCDLen, cOut.ModLen()
	if inLen == 0 || inLen != len(sc.vecSc.modIn) || outLen != len(sc.vecSc.modOut) {
		panic("inconsistent input(s)")
	}

	cBufPtr := sc.pool.Get()
	cBuf := (*cBufPtr)[:inLen]
	defer sc.pool.Put(cBufPtr)

	if inLen == outLen && inLen == gcdLen {
		for i := 0; i < inLen; i++ {
			cBuf[i] = c.Value[sc.vecSc.idxGCDToIn[sc.vecSc.idxOutToGCD[i]]]
		}
		copy(cOut.Value, cBuf)
		return
	}

	if inLen == gcdLen {
		for i := 0; i < inLen; i++ {
			ii := sc.vecSc.idxInToGCD[i]
			cBuf[ii] = num.SMul(c.Value[sc.vecSc.modInMap[i]], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], sc.vecSc.modInSorted[i])
		}
		for i := 0; i < outLen; i++ {
			ii := sc.vecSc.idxOutToGCD[i]
			if ii >= 0 {
				cOut.Value[i] = cBuf[ii]
			} else {
				cOut.Value[i] = 0
			}
		}
		return
	}

	cMulPtr := sc.pool.Get()
	cMul := (*cMulPtr)[:gcdLen]
	defer sc.pool.Put(cMulPtr)

	cBuf = cBuf[:inLen-gcdLen]
	for i := 0; i < inLen; i++ {
		ii := sc.vecSc.idxInToGCD[i]
		if ii >= 0 {
			cMul[ii] = num.SMul(c.Value[sc.vecSc.modInMap[i]], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], sc.vecSc.modInSorted[i])
		}
	}
	for i := 0; i < inLen; i++ {
		ii := sc.vecSc.idxInToComp[i]
		if ii >= 0 {
			cBuf[ii] = num.SMul(c.Value[sc.vecSc.modInMap[i]], sc.vecSc.outCompModIn[i], sc.vecSc.outCompModInS[i], sc.vecSc.modInSorted[i])
		}
	}

	if inLen-gcdLen == 1 {
		x := cBuf[0]
		qv := sc.vecSc.modInComp[0].Value()
		halfQv := qv >> 1
		for i := 0; i < outLen; i++ {
			qOut := sc.vecSc.modOut[i].Value()
			qOutDivHi, _ := sc.vecSc.modOut[i].Div()
			cOut.Value[i] = embedToModOut(x, qOut, qOutDivHi, qv, halfQv)
		}
	} else {
		for i := 0; i < inLen-gcdLen; i++ {
			for j := i + 1; j < inLen-gcdLen; j++ {
				cBuf[j] = num.Sub(cBuf[j], cBuf[i], sc.vecSc.modInComp[j])
				cBuf[j] = num.SMul(cBuf[j], sc.vecSc.modInCompInv[i][j-i-1], sc.vecSc.modInCompInvS[i][j-i-1], sc.vecSc.modInComp[j])
			}
		}

		isNegative := false
		for i := inLen - gcdLen - 1; i >= 0; i-- {
			if cBuf[i] > sc.vecSc.modInCompHalf[i] {
				isNegative = true
				break
			}
			if cBuf[i] < sc.vecSc.modInCompHalf[i] {
				break
			}
		}

		for i := 0; i < outLen; i++ {
			base, baseS := sc.vecSc.base[i], sc.vecSc.baseS[i]
			modOut := sc.vecSc.modOut[i]
			x := num.SMul(cBuf[0], base[0], baseS[0], modOut)
			for j := 1; j < inLen-gcdLen; j++ {
				x = num.Add(x, num.SMul(cBuf[j], base[j], baseS[j], modOut), modOut)
			}
			if isNegative {
				x = num.Sub(x, sc.vecSc.inCompModOut[i], modOut)
			}
			cOut.Value[i] = x
		}
	}

	for i := 0; i < outLen; i++ {
		modOut := sc.vecSc.modOut[i]
		ii := sc.vecSc.idxOutToGCD[i]
		if ii >= 0 {
			cOut.Value[i] = num.Sub(cOut.Value[i], cMul[ii], modOut)
		}
		cOut.Value[i] = num.SMul(cOut.Value[i], sc.vecSc.negInCompInvModOut[i], sc.vecSc.negInCompInvModOutS[i], modOut)
	}
}

// VecScaler scales a vector to different modulus.
// In other words, it computes
//
//	[v]_modIn -> [(modOut / modIn) * v]_modOut
//
// It uses mixed-radix representation conversion, so the computation is exact.
// Each modulus from modIn and modOut should be equal or coprime.
type VecScaler struct {
	// modIn is the input modulus order.
	modIn []*num.Modulus
	// modInSorted is the input modulus sorted for mixed-radix conversion.
	modInSorted []*num.Modulus
	// modOut is the output modulus.
	modOut []*num.Modulus
	// modInComp equals modInSorted/modGCD.
	modInComp []*num.Modulus

	// modInMap maps the index of modInSorted to modIn.
	modInMap []int

	// modInCompHalf is modInComp/2 in the mixed-radix basis of modInComp.
	modInCompHalf []uint64

	// modGCDLen is the length of the GCD of modIn and modOut.
	modGCDLen int

	// idxInToGCD maps the index of modInSorted to the index of modGCD.
	// If modInSorted[i] is not in modGCD, idxInToGCD[i] is -1.
	idxInToGCD []int
	// idxInToComp maps the index of modInSorted to the index of modInComp.
	// If modInSorted[i] is in modGCD, idxInToComp[i] is -1.
	idxInToComp []int
	// idxOutToGCD maps the index of modOut to the index of modGCD.
	// If modOut[i] is not in modGCD, idxOutToGCD[i] is -1.
	idxOutToGCD []int
	// idxGCDToIn maps the index of modGCD to the index of modIn.
	idxGCDToIn []int

	// modInCompInv is the inverse of modInComp.
	modInCompInv [][]uint64
	// modInCompInvS is the Shoup form of modInCompInv.
	modInCompInvS [][]uint64

	// base is the mixed-radix basis of modInComp.
	base [][]uint64
	// baseS is the Shoup form of base.
	baseS [][]uint64

	// inCompModOut equals modInComp modulo modOut.
	inCompModOut []uint64

	// outCompModIn equals modOut/modGCD modulo modInSorted.
	outCompModIn []uint64
	// outCompModInS is the Shoup form of outCompModIn.
	outCompModInS []uint64
	// negInCompInvModOut equals the negative inverse of modIn/modGCD modulo modOut.
	negInCompInvModOut []uint64
	// negInCompInvModOutS is the Shoup form of negInCompInvModOut.
	negInCompInvModOutS []uint64
}

// NewVecScaler creates a new [VecScaler].
func NewVecScaler(modOut, modIn []*num.Modulus) *VecScaler {
	modInSorted := make([]*num.Modulus, len(modIn))
	copy(modInSorted, modIn)
	slices.SortFunc(modInSorted, num.CmpModulus)

	modInMap := make([]int, len(modInSorted))
	for i := range modInSorted {
		for j := range modIn {
			if modInSorted[i].Value() == modIn[j].Value() {
				modInMap[i] = j
				break
			}
		}
	}

	modGCD := make([]*num.Modulus, 0, min(len(modOut), len(modIn)))
	isGCDOut := make([]bool, len(modOut))
	isGCDIn := make([]bool, len(modInSorted))
	for i := 0; i < len(modOut); i++ {
		for j := 0; j < len(modInSorted); j++ {
			if modOut[i].Value() == modInSorted[j].Value() {
				modGCD = append(modGCD, modOut[i])
				isGCDOut[i] = true
				isGCDIn[j] = true
				break
			}
		}
	}

	var idxIn int
	idxInToGCD := make([]int, len(modInSorted))
	idxInToComp := make([]int, len(modInSorted))
	for i := 0; i < len(modInSorted); i++ {
		idxInToGCD[i] = -1
		idxInToComp[i] = -1
		for j := 0; j < len(modGCD); j++ {
			if modInSorted[i].Value() == modGCD[j].Value() {
				idxInToGCD[i] = j
				break
			}
		}
		if idxInToGCD[i] == -1 {
			idxInToComp[i] = idxIn
			idxIn++
		}
	}

	idxGCDToIn := make([]int, len(modGCD))
	for i := range idxGCDToIn {
		idxGCDToIn[i] = -1
	}
	for i, idx := range idxInToGCD {
		if idx >= 0 {
			idxGCDToIn[idx] = modInMap[i]
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

	modInComp := make([]*num.Modulus, 0, len(modInSorted)-len(modGCD))
	for i := 0; i < len(modInSorted); i++ {
		if idxInToGCD[i] == -1 {
			modInComp = append(modInComp, modInSorted[i])
		}
	}
	modInCompHalf := halfProductMixedRadix(modInComp)

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
		for j := 0; j < len(modInComp); j++ {
			if j == 0 {
				base[i][j] = 1
			} else {
				base[i][j] = num.Mul(base[i][j-1], modInComp[j-1].Value(), modOut[i])
			}
		}
		baseS[i] = vec.SForm(base[i], modOut[i])
	}

	inCompModOut := make([]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		inCompModOut[i] = 1
		for j := 0; j < len(modInComp); j++ {
			inCompModOut[i] = num.Mul(inCompModOut[i], modInComp[j].Value(), modOut[i])
		}
	}

	outCompModIn := make([]uint64, len(modInSorted))
	outCompModInS := make([]uint64, len(modInSorted))
	for i := 0; i < len(modInSorted); i++ {
		outCompModIn[i] = 1
		for j := 0; j < len(modOut); j++ {
			if !isGCDOut[j] {
				outCompModIn[i] = num.Mul(outCompModIn[i], modOut[j].Value(), modInSorted[i])
			}
		}
		outCompModInS[i] = num.SForm(outCompModIn[i], modInSorted[i])
	}

	negInCompInvModOut := make([]uint64, len(modOut))
	negInCompInvModOutS := make([]uint64, len(modOut))
	for i := 0; i < len(modOut); i++ {
		negInCompInvModOut[i] = 1
		for j := 0; j < len(modInSorted); j++ {
			if !isGCDIn[j] {
				negInCompInvModOut[i] = num.Mul(negInCompInvModOut[i], num.Inv(modInSorted[j].Value(), modOut[i]), modOut[i])
			}
		}
		negInCompInvModOut[i] = num.Neg(negInCompInvModOut[i], modOut[i])
		negInCompInvModOutS[i] = num.SForm(negInCompInvModOut[i], modOut[i])
	}

	return &VecScaler{
		modIn:       modIn,
		modInSorted: modInSorted,
		modOut:      modOut,
		modInComp:   modInComp,

		modInMap: modInMap,

		modInCompHalf: modInCompHalf,

		modGCDLen: len(modInSorted) - len(modInComp),

		idxInToGCD:  idxInToGCD,
		idxInToComp: idxInToComp,
		idxOutToGCD: idxOutToGCD,
		idxGCDToIn:  idxGCDToIn,

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

// Scale returns the scaled vector of v.
func (sc *VecScaler) Scale(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(sc.modOut))
	for i := 0; i < len(sc.modOut); i++ {
		vOut[i] = make([]uint64, len(v[0]))
	}
	sc.ScaleTo(vOut, v)
	return vOut
}

// ScaleTo scales v to vOut.
//
// Panics when v and vOut is inconsistent.
func (sc *VecScaler) ScaleTo(vOut, v [][]uint64) {
	M := (len(v[0]) >> logEmbedBatch) << logEmbedBatch
	L := unsafe.Sizeof(uint64(0))

	inLen, gcdLen, outLen := len(sc.modInSorted), sc.modGCDLen, len(sc.modOut)
	if len(v) != len(sc.modIn) || len(vOut) != len(sc.modOut) {
		panic("inconsistent input(s)")
	}

	N := len(v[0])
	for i := 1; i < inLen; i++ {
		if len(v[i]) != N {
			panic("inconsistent input(s)")
		}
	}
	for i := 0; i < outLen; i++ {
		if len(vOut[i]) != N {
			panic("inconsistent input(s)")
		}
	}

	if inLen == outLen && inLen == gcdLen {
		vBuf := make([]*[embedBatch]uint64, inLen)
		for i := range vBuf {
			vBuf[i] = embed64Pool.Get()
		}
		defer func() {
			for i := range vBuf {
				embed64Pool.Put(vBuf[i])
			}
		}()

		for k := 0; k < M; k += embedBatch {
			for i := 0; i < inLen; i++ {
				r := unsafe.Pointer(unsafe.SliceData(v[sc.idxGCDToIn[sc.idxOutToGCD[i]]]))
				w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))
				copy(vBuf[i][:], w[:])
			}

			for i := 0; i < outLen; i++ {
				rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
				wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))

				copy(wOut[:], vBuf[i][:])
			}
		}

		for i := 0; i < inLen; i++ {
			copy(vBuf[i][:len(v[0])-M], v[sc.idxGCDToIn[sc.idxOutToGCD[i]]][M:])
		}

		for i := 0; i < outLen; i++ {
			copy(vOut[i][M:], vBuf[i][:len(v[0])-M])
		}

		return
	}

	if len(sc.modInSorted) == sc.modGCDLen {
		vBuf := make([]*[embedBatch]uint64, inLen)
		for i := range vBuf {
			vBuf[i] = embed64Pool.Get()
		}
		defer func() {
			for i := range vBuf {
				embed64Pool.Put(vBuf[i])
			}
		}()

		for k := 0; k < M; k += embedBatch {
			for i := 0; i < inLen; i++ {
				r := unsafe.Pointer(unsafe.SliceData(v[sc.modInMap[i]]))
				w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))

				ii := sc.idxInToGCD[i]
				vec.SMulTo(vBuf[ii][:], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modInSorted[i])
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
			w := v[sc.modInMap[i]][M:]

			ii := sc.idxInToGCD[i]
			vec.SMulTo(vBuf[ii][:len(v[0])-M], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modInSorted[i])
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
	}
	defer func() {
		for i := range vMul {
			embed64Pool.Put(vMul[i])
		}
	}()

	vBuf := make([]*[embedBatch]uint64, inLen-gcdLen)
	for i := range vBuf {
		vBuf[i] = embed64Pool.Get()
	}
	defer func() {
		for i := range vBuf {
			embed64Pool.Put(vBuf[i])
		}
	}()

	vBool := embed64Pool.Get()
	defer embed64Pool.Put(vBool)
	vCorr := embed64Pool.Get()
	defer embed64Pool.Put(vCorr)

	var qLast uint64
	for i := inLen - 1; i >= 0; i-- {
		if sc.idxInToComp[i] >= 0 {
			qLast = sc.modInSorted[i].Value()
			break
		}
	}
	qLastHalf := qLast >> 1

	for k := 0; k < M; k += embedBatch {
		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[sc.modInMap[i]]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))

			modIn := sc.modInSorted[i]

			ii := sc.idxInToGCD[i]
			if ii >= 0 {
				vec.SMulTo(vMul[ii][:], w[:], sc.outCompModIn[i], sc.outCompModInS[i], modIn)
			}
		}

		for i := 0; i < inLen; i++ {
			r := unsafe.Pointer(unsafe.SliceData(v[sc.modInMap[i]]))
			w := (*[embedBatch]uint64)(unsafe.Add(r, uintptr(k)*L))

			modIn := sc.modInSorted[i]

			ii := sc.idxInToComp[i]
			if ii >= 0 {
				vec.SMulTo(vBuf[ii][:], w[:], sc.outCompModIn[i], sc.outCompModInS[i], modIn)
			}
		}

		if inLen-gcdLen == 1 {
			for i := 0; i < outLen; i++ {
				rOut := unsafe.Pointer(unsafe.SliceData(vOut[i]))
				wOut := (*[embedBatch]uint64)(unsafe.Add(rOut, uintptr(k)*L))
				wBuf := vBuf[0]

				qOut := sc.modOut[i].Value()
				qOutDivHi, _ := sc.modOut[i].Div()

				for j := 0; j < embedBatch; j += 8 {
					wOut[j+0] = embedToModOut(wBuf[j+0], qOut, qOutDivHi, qLast, qLastHalf)
					wOut[j+1] = embedToModOut(wBuf[j+1], qOut, qOutDivHi, qLast, qLastHalf)
					wOut[j+2] = embedToModOut(wBuf[j+2], qOut, qOutDivHi, qLast, qLastHalf)
					wOut[j+3] = embedToModOut(wBuf[j+3], qOut, qOutDivHi, qLast, qLastHalf)

					wOut[j+4] = embedToModOut(wBuf[j+4], qOut, qOutDivHi, qLast, qLastHalf)
					wOut[j+5] = embedToModOut(wBuf[j+5], qOut, qOutDivHi, qLast, qLastHalf)
					wOut[j+6] = embedToModOut(wBuf[j+6], qOut, qOutDivHi, qLast, qLastHalf)
					wOut[j+7] = embedToModOut(wBuf[j+7], qOut, qOutDivHi, qLast, qLastHalf)
				}
			}
		} else {
			for i := 0; i < inLen-gcdLen; i++ {
				wBuf := vBuf[i]
				for j := i + 1; j < inLen-gcdLen; j++ {
					modInComp := sc.modInComp[j]
					vec.SubTo(vBuf[j][:], vBuf[j][:], wBuf[:], modInComp)
					vec.SMulTo(vBuf[j][:], vBuf[j][:], sc.modInCompInv[i][j-i-1], sc.modInCompInvS[i][j-i-1], modInComp)
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

				vec.SMulTo(wOut[:], vBuf[0][:], base[0], baseS[0], modOut)
				for j := 1; j < inLen-gcdLen; j++ {
					vec.SMulAddTo(wOut[:], vBuf[j][:], base[j], baseS[j], modOut)
				}
				vec.MulTo(vCorr[:], vBool[:], inCompModOut, nil)
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
			vec.SMulTo(wOut[:], wOut[:], sc.negInCompInvModOut[i], sc.negInCompInvModOutS[i], modOut)
		}
	}

	for i := 0; i < inLen; i++ {
		w := v[sc.modInMap[i]][M:]

		ii := sc.idxInToGCD[i]
		if ii >= 0 {
			vec.SMulTo(vMul[ii][:len(v[0])-M], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modInSorted[i])
		}
	}

	for i := 0; i < inLen; i++ {
		w := v[sc.modInMap[i]][M:]

		ii := sc.idxInToComp[i]
		if ii >= 0 {
			vec.SMulTo(vBuf[ii][:len(v[0])-M], w[:], sc.outCompModIn[i], sc.outCompModInS[i], sc.modInSorted[i])
		}
	}

	if inLen-gcdLen == 1 {
		for i := 0; i < outLen; i++ {
			for j := 0; j < len(v[0])-M; j++ {
				qOut := sc.modOut[i].Value()
				qOutDivHi, _ := sc.modOut[i].Div()
				vOut[i][j+M] = embedToModOut(vBuf[0][j], qOut, qOutDivHi, qLast, qLastHalf)
			}
		}
	} else {
		for i := 0; i < inLen-gcdLen; i++ {
			wBuf := vBuf[i][:len(v[0])-M]
			for j := i + 1; j < inLen-gcdLen; j++ {
				vec.SubTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], wBuf[:], sc.modInComp[j])
				vec.SMulTo(vBuf[j][:len(v[0])-M], vBuf[j][:len(v[0])-M], sc.modInCompInv[i][j-i-1], sc.modInCompInvS[i][j-i-1], sc.modInComp[j])
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

			vec.SMulTo(wOut[:], vBuf[0][:len(v[0])-M], base[0], baseS[0], modOut)
			for j := 1; j < inLen-gcdLen; j++ {
				vec.SMulAddTo(wOut[:], vBuf[j][:len(v[0])-M], base[j], baseS[j], modOut)
			}
			vec.MulTo(vCorr[:len(v[0])-M], vBool[:len(v[0])-M], inCompModOut, nil)
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
		vec.SMulTo(wOut[:], wOut[:], sc.negInCompInvModOut[i], sc.negInCompInvModOutS[i], modOut)
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
