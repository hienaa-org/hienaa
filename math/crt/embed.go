package crt

import (
	"math/big"
	"unsafe"

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

	// invHi is the high 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invHi []uint64
	// invLo is the low 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invLo []uint64

	// idx holds the index of the input modulus limb if it overlaps with the output modulus limb.
	// For example, if modOut[i] = modIn[j], then idx[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	idx []int

	buf embedderBuffer
}

// embedderBuffer is a buffer for [Embedder].
type embedderBuffer struct {
	// fHi is a buffer for the high 64 bits of the fixed-point number.
	// Always has length 8.
	fHi []uint64
	// fLo is a buffer for the low 64 bits of the fixed-point number.
	// Always has length 8.
	fLo []uint64
	// f64 is a buffer for the 64-bit fixed-point number.
	// Always has length 8.
	f64 []uint64
	// in is a buffer for input coefficient.
	// Always has length [len(modIn)][8].
	in [][]uint64
}

// NewEmbedder creates a new [Embedder].
func NewEmbedder(modOut []*num.Modulus, modIn []*num.Modulus) *Embedder {
	switch {
	case len(modIn) == 0:
		panic("NewEmbedder: modIn cannot be empty")
	case len(modOut) == 0:
		panic("NewEmbedder: modOut cannot be empty")
	}

	compInv := make([]uint64, len(modIn))
	compInvS := make([]uint64, len(modIn))

	invHi := make([]uint64, len(modIn))
	invLo := make([]uint64, len(modIn))

	tmpInt := big.NewInt(0)
	tmpRat := big.NewRat(0, 1)
	for i := 0; i < len(modIn); i++ {
		compInv[i] = 1
		for j := 0; j < len(modIn); j++ {
			if i != j {
				compInv[i] = num.Mul(compInv[i], num.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
			compInvS[i] = num.SForm(compInv[i], modIn[i])
		}

		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Num().Lsh(tmpRat.Num().SetInt64(1), fixedPrec-64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		invHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		invLo[i] = tmpInt.Uint64()
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

		invHi: invHi,
		invLo: invLo,

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
		fHi: make([]uint64, 8),
		fLo: make([]uint64, 8),
		f64: make([]uint64, 8),
		in:  in,
	}
}

// Embed returns the embedding of p to the output modulus.
// If p.ModLen() < len(e.modIn), it only embeds the first p.ModLen() elements.
func (e *Embedder) Embed(p *Poly) *Poly {
	pOut := NewPoly(p.Rank(), p.ModLen())
	e.EmbedTo(pOut, p)
	return pOut
}

// EmbedTo embeds the p to pOut.
// If p.ModLen() < len(e.modIn) or pOut.ModLen() < len(e.modOut),
// it only embeds the first p.ModLen() elements to pOut.ModLen() elements.
func (e *Embedder) EmbedTo(pOut, p *Poly) {
	if p.isNTT || pOut.isNTT {
		panic("Embed: cannot embed NTT polynomials")
	}
	e.EmbedVecTo(p.Coeffs, pOut.Coeffs)
}

// EmbedVec returns the embedding of v to the output modulus.
// If len(v) < len(e.modIn), it only embeds the first len(v) elements.
func (e *Embedder) EmbedVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(e.modOut))
	for i := 0; i < len(e.modOut); i++ {
		vOut[i] = make([]uint64, len(v[i]))
	}
	e.EmbedVecTo(vOut, v)
	return vOut
}

// EmbedVecTo embeds v to vOut.
// If len(v) < len(e.modIn) or len(vOut) < len(e.modOut),
// it only embeds the first len(v) elements to len(vOut) elements.
func (e *Embedder) EmbedVecTo(vOut, v [][]uint64) {
	M := (len(v[0]) >> 3) << 3

	inLen, outLen := min(len(v), len(e.modIn)), min(len(vOut), len(e.modOut))

	if inLen == 1 {
		qv := e.modIn[0].Value()
		halfQv := qv >> 1
		for i := 0; i < outLen; i++ {
			if e.idx[i] == 0 {
				copy(vOut[i], v[0])
			} else {
				for j := 0; j < M; j += 8 {
					wIn := (*[8]uint64)(unsafe.Pointer(&v[0][j]))
					wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][j]))

					wOut[0] = reduceModInToModOutSigned(wIn[0], e.modOut[i], qv, halfQv)
					wOut[1] = reduceModInToModOutSigned(wIn[1], e.modOut[i], qv, halfQv)
					wOut[2] = reduceModInToModOutSigned(wIn[2], e.modOut[i], qv, halfQv)
					wOut[3] = reduceModInToModOutSigned(wIn[3], e.modOut[i], qv, halfQv)

					wOut[4] = reduceModInToModOutSigned(wIn[4], e.modOut[i], qv, halfQv)
					wOut[5] = reduceModInToModOutSigned(wIn[5], e.modOut[i], qv, halfQv)
					wOut[6] = reduceModInToModOutSigned(wIn[6], e.modOut[i], qv, halfQv)
					wOut[7] = reduceModInToModOutSigned(wIn[7], e.modOut[i], qv, halfQv)
				}

				for j := M; j < len(v[0]); j++ {
					vOut[i][j] = reduceModInToModOutSigned(v[0][j], e.modOut[i], qv, halfQv)
				}
			}
		}
		return
	}

	fHi := (*[8]uint64)(unsafe.Pointer(&e.buf.fHi[0]))
	fLo := (*[8]uint64)(unsafe.Pointer(&e.buf.fLo[0]))
	var hi, lo uint64

	for k := 0; k < M; k += 8 {
		clear(fHi[:])
		clear(fLo[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&e.buf.in[i][0]))

			compInv, compInvS := e.compInv[i], e.compInvS[i]
			invHi, invLo := e.invHi[i], e.invLo[i]
			modIn := e.modIn[i]

			bufIn[0] = num.SMul(wIn[0], compInv, compInvS, modIn)
			bufIn[1] = num.SMul(wIn[1], compInv, compInvS, modIn)
			bufIn[2] = num.SMul(wIn[2], compInv, compInvS, modIn)
			bufIn[3] = num.SMul(wIn[3], compInv, compInvS, modIn)

			bufIn[4] = num.SMul(wIn[4], compInv, compInvS, modIn)
			bufIn[5] = num.SMul(wIn[5], compInv, compInvS, modIn)
			bufIn[6] = num.SMul(wIn[6], compInv, compInvS, modIn)
			bufIn[7] = num.SMul(wIn[7], compInv, compInvS, modIn)

			hi, lo = mulAndFloor(bufIn[0], invHi, invLo)
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)

			hi, lo = mulAndFloor(bufIn[1], invHi, invLo)
			fHi[1], fLo[1] = add128(fHi[1], fLo[1], hi, lo)

			hi, lo = mulAndFloor(bufIn[2], invHi, invLo)
			fHi[2], fLo[2] = add128(fHi[2], fLo[2], hi, lo)

			hi, lo = mulAndFloor(bufIn[3], invHi, invLo)
			fHi[3], fLo[3] = add128(fHi[3], fLo[3], hi, lo)

			hi, lo = mulAndFloor(bufIn[4], invHi, invLo)
			fHi[4], fLo[4] = add128(fHi[4], fLo[4], hi, lo)

			hi, lo = mulAndFloor(bufIn[5], invHi, invLo)
			fHi[5], fLo[5] = add128(fHi[5], fLo[5], hi, lo)

			hi, lo = mulAndFloor(bufIn[6], invHi, invLo)
			fHi[6], fLo[6] = add128(fHi[6], fLo[6], hi, lo)

			hi, lo = mulAndFloor(bufIn[7], invHi, invLo)
			fHi[7], fLo[7] = add128(fHi[7], fLo[7], hi, lo)
		}

		fLo[0] = roundTo64(fHi[0], fLo[0])
		fLo[1] = roundTo64(fHi[1], fLo[1])
		fLo[2] = roundTo64(fHi[2], fLo[2])
		fLo[3] = roundTo64(fHi[3], fLo[3])

		fLo[4] = roundTo64(fHi[4], fLo[4])
		fLo[5] = roundTo64(fHi[5], fLo[5])
		fLo[6] = roundTo64(fHi[6], fLo[6])
		fLo[7] = roundTo64(fHi[7], fLo[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			if 0 <= e.idx[i] && e.idx[i] < inLen {
				copy(wOut[:], v[e.idx[i]][k:k+8])
			} else {
				negMod, negModS := e.negMod[i], e.negModS[i]
				modOut := e.modOut[i]

				wOut[0] = num.SMul(fLo[0], negMod, negModS, modOut)
				wOut[1] = num.SMul(fLo[1], negMod, negModS, modOut)
				wOut[2] = num.SMul(fLo[2], negMod, negModS, modOut)
				wOut[3] = num.SMul(fLo[3], negMod, negModS, modOut)

				wOut[4] = num.SMul(fLo[4], negMod, negModS, modOut)
				wOut[5] = num.SMul(fLo[5], negMod, negModS, modOut)
				wOut[6] = num.SMul(fLo[6], negMod, negModS, modOut)
				wOut[7] = num.SMul(fLo[7], negMod, negModS, modOut)

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
		fHi[0], fLo[0] = 0, 0
		for i := 0; i < inLen; i++ {
			e.buf.in[i][0] = num.SMul(v[i][k], e.compInv[i], e.compInvS[i], e.modIn[i])

			hi, lo := mulAndFloor(e.buf.in[i][0], e.invHi[i], e.invLo[i])
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)
		}
		fLo[0] = roundTo64(fHi[0], fLo[0])

		for i := 0; i < outLen; i++ {
			if 0 <= e.idx[i] && e.idx[i] < inLen {
				vOut[i][k] = v[e.idx[i]][k]
			} else {
				vOut[i][k] = num.SMul(fLo[0], e.negMod[i], e.negModS[i], e.modOut[i])
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

		invHi: e.invHi,
		invLo: e.invLo,

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

	// scFracHi is the high 64 bits of the fixed-point approximation of the fractional part of the modulus scaling factor.
	scFracHi []uint64
	// scFracLo is the low 64 bits of the fixed-point approximation of the fractional part of the modulus scaling factor.
	scFracLo []uint64

	buf embedderBuffer
}

// NewScaler creates a new [Scaler].
func NewScaler(modOut []*num.Modulus, modIn []*num.Modulus) *Scaler {
	switch {
	case len(modIn) == 0:
		panic("NewScaler: modIn cannot be empty")
	case len(modOut) == 0:
		panic("NewScaler: modOut cannot be empty")
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

	scFracHi := make([]uint64, len(modIn))
	scFracLo := make([]uint64, len(modIn))

	modOutBig := big.NewInt(1)
	tmpInt := big.NewInt(0)
	for i := 0; i < len(modOut); i++ {
		modOutBig.Mul(modOutBig, tmpInt.SetUint64(modOut[i].Value()))
	}

	tmpRat := big.NewRat(0, 1)
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

		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Num().Set(tmpInt.Lsh(tmpInt, fixedPrec-64))

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		scFracHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		scFracLo[i] = tmpInt.Uint64()
	}

	return &Scaler{
		modIn:  modIn,
		modOut: modOut,

		compInv:  compInv,
		compInvS: compInvS,

		scInt:  scInt,
		scIntS: scIntS,

		scFracHi: scFracHi,
		scFracLo: scFracLo,

		buf: newEmbedderBuffer(modIn),
	}
}

// Scale returns the scaled polynomial of p.
func (s *Scaler) Scale(p *Poly) *Poly {
	pOut := NewPoly(p.Rank(), p.ModLen())
	s.ScaleTo(pOut, p)
	return pOut
}

// ScaleTo scales p to pOut.
func (s *Scaler) ScaleTo(pOut, p *Poly) {
	switch {
	case p.isNTT || pOut.isNTT:
		panic("ScaleTo: cannot scale NTT polynomials")
	case p.ModLen() != len(s.modIn):
		panic("ScaleTo: len(p.modulus) != len(s.ModIn)")
	}

	s.ScaleVecTo(p.Coeffs, pOut.Coeffs)
}

// ScaleVec returns the scaled vector of v.
func (s *Scaler) ScaleVec(v [][]uint64) [][]uint64 {
	vOut := make([][]uint64, len(s.modOut))
	for i := 0; i < len(s.modOut); i++ {
		vOut[i] = make([]uint64, len(v[i]))
	}
	s.ScaleVecTo(vOut, v)
	return vOut
}

// ScaleVecTo scales v to vOut.
func (s *Scaler) ScaleVecTo(vOut, v [][]uint64) {
	inLen, outLen := len(s.modIn), len(s.modOut)
	M := (len(v[0]) >> 3) << 3

	switch {
	case inLen != len(v):
		panic("ScaleTo: len(in) != len(s.modIn)")
	case outLen != len(vOut):
		panic("ScaleTo: len(out) != len(s.modOut)")
	}

	var hi, lo uint64
	fHi := (*[8]uint64)(unsafe.Pointer(&s.buf.fHi[0]))
	fLo := (*[8]uint64)(unsafe.Pointer(&s.buf.fLo[0]))

	for k := 0; k < M; k += 8 {
		clear(s.buf.fHi[:])
		clear(s.buf.fLo[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[i][0]))

			compInv, compInvS := s.compInv[i], s.compInvS[i]
			scFracHi, scFracLo := s.scFracHi[i], s.scFracLo[i]
			modIn := s.modIn[i]

			bufIn[0] = num.SMul(wIn[0], compInv, compInvS, modIn)
			bufIn[1] = num.SMul(wIn[1], compInv, compInvS, modIn)
			bufIn[2] = num.SMul(wIn[2], compInv, compInvS, modIn)
			bufIn[3] = num.SMul(wIn[3], compInv, compInvS, modIn)

			bufIn[4] = num.SMul(wIn[4], compInv, compInvS, modIn)
			bufIn[5] = num.SMul(wIn[5], compInv, compInvS, modIn)
			bufIn[6] = num.SMul(wIn[6], compInv, compInvS, modIn)
			bufIn[7] = num.SMul(wIn[7], compInv, compInvS, modIn)

			hi, lo = mulAndFloor(bufIn[0], scFracHi, scFracLo)
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)

			hi, lo = mulAndFloor(bufIn[1], scFracHi, scFracLo)
			fHi[1], fLo[1] = add128(fHi[1], fLo[1], hi, lo)

			hi, lo = mulAndFloor(bufIn[2], scFracHi, scFracLo)
			fHi[2], fLo[2] = add128(fHi[2], fLo[2], hi, lo)

			hi, lo = mulAndFloor(bufIn[3], scFracHi, scFracLo)
			fHi[3], fLo[3] = add128(fHi[3], fLo[3], hi, lo)

			hi, lo = mulAndFloor(bufIn[4], scFracHi, scFracLo)
			fHi[4], fLo[4] = add128(fHi[4], fLo[4], hi, lo)

			hi, lo = mulAndFloor(bufIn[5], scFracHi, scFracLo)
			fHi[5], fLo[5] = add128(fHi[5], fLo[5], hi, lo)

			hi, lo = mulAndFloor(bufIn[6], scFracHi, scFracLo)
			fHi[6], fLo[6] = add128(fHi[6], fLo[6], hi, lo)

			hi, lo = mulAndFloor(bufIn[7], scFracHi, scFracLo)
			fHi[7], fLo[7] = add128(fHi[7], fLo[7], hi, lo)
		}

		fHi[0], fLo[0] = roundTo128(fHi[0], fLo[0])
		fHi[1], fLo[1] = roundTo128(fHi[1], fLo[1])
		fHi[2], fLo[2] = roundTo128(fHi[2], fLo[2])
		fHi[3], fLo[3] = roundTo128(fHi[3], fLo[3])

		fHi[4], fLo[4] = roundTo128(fHi[4], fLo[4])
		fHi[5], fLo[5] = roundTo128(fHi[5], fLo[5])
		fHi[6], fLo[6] = roundTo128(fHi[6], fLo[6])
		fHi[7], fLo[7] = roundTo128(fHi[7], fLo[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			modOut := s.modOut[i]

			wOut[0] = num.Reduce128(fHi[0], fLo[0], modOut)
			wOut[1] = num.Reduce128(fHi[1], fLo[1], modOut)
			wOut[2] = num.Reduce128(fHi[2], fLo[2], modOut)
			wOut[3] = num.Reduce128(fHi[3], fLo[3], modOut)

			wOut[4] = num.Reduce128(fHi[4], fLo[4], modOut)
			wOut[5] = num.Reduce128(fHi[5], fLo[5], modOut)
			wOut[6] = num.Reduce128(fHi[6], fLo[6], modOut)
			wOut[7] = num.Reduce128(fHi[7], fLo[7], modOut)

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
		fHi[0], fLo[0] = 0, 0
		for i := 0; i < inLen; i++ {
			s.buf.in[i][0] = num.SMul(v[i][k], s.compInv[i], s.compInvS[i], s.modIn[i])
			hi, lo = mulAndFloor(s.buf.in[i][0], s.scFracHi[i], s.scFracLo[i])
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)
		}
		fHi[0], fLo[0] = roundTo128(fHi[0], fLo[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.Reduce128(fHi[0], fLo[0], s.modOut[i])
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

		scFracHi: s.scFracHi,
		scFracLo: s.scFracLo,

		buf: newEmbedderBuffer(s.modIn),
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

	// invHi is the high 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invHi []uint64
	// invLo is the low 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invLo []uint64

	// scInt is the integer part of the modulus scaling factor.
	scInt [][]uint64
	// scIntS is the Shoup form of modScInt.
	scIntS [][]uint64

	// scFracHi is the high 64 bits of the fixed-point approximation of the fractional part of the scaling factor.
	scFracHi []uint64
	// scFracLo is the low 64 bits of the fixed-point approximation of the fractional part of the scaling factor.
	scFracLo []uint64

	// ovfInt is the integer part for the constant to compute the overflow multiplied by the scaling factor.
	ovfInt []uint64
	// ovfIntS is the Shoup form of ovfInt.
	ovfIntS []uint64

	// ovfFracHi is the high 64 bits of the fixed-point approximation of the fractional part
	// for the constant to compute the overflow multiplied by the scaling factor.
	ovfFracHi uint64
	// ovfFracLo is the low 64 bits of the fixed-point approximation of the fractional part
	// for the constant to compute the overflow multiplied by the scaling factor.
	ovfFracLo uint64

	buf embedderBuffer
}

// NewScaleEmbedder creates a new [ScaleEmbedder].
func NewScaleEmbedder(scale *big.Rat, modOut []*num.Modulus, modIn []*num.Modulus) *ScaleEmbedder {
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

	invHi := make([]uint64, inLen)
	invLo := make([]uint64, inLen)

	modInBig := big.NewInt(1)
	tmpInt := big.NewInt(0)
	tmpRat := big.NewRat(0, 1)
	modScIntBig := big.NewInt(0)

	for i := 0; i < inLen; i++ {
		modInBig.Mul(modInBig, tmpInt.SetUint64(modIn[i].Value()))

		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Num().Lsh(tmpRat.Num().SetInt64(1), fixedPrec-64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		invHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		invLo[i] = tmpInt.Uint64()
	}

	scInt := make([][]uint64, outLen)
	scIntS := make([][]uint64, outLen)
	for i := 0; i < outLen; i++ {
		scInt[i] = make([]uint64, inLen)
		scIntS[i] = make([]uint64, inLen)
	}

	scFracHi := make([]uint64, inLen)
	scFracLo := make([]uint64, inLen)

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

		tmpRat.Num().Lsh(tmpRat.Num(), fixedPrec-64)
		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		scFracHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)
		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		scFracLo[i] = tmpInt.Uint64()
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
	tmpRat.Num().Lsh(tmpRat.Num(), fixedPrec-64)
	tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
	ovfFracHi := tmpInt.Uint64()

	tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
	tmpRat.Num().Lsh(tmpRat.Num(), 64)
	tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
	ovfFracLo := tmpInt.Uint64()

	return &ScaleEmbedder{
		modIn:  modIn,
		modOut: modOut,
		scale:  scale,

		compInv:  compInv,
		compInvS: compInvS,

		invHi: invHi,
		invLo: invLo,

		scInt:  scInt,
		scIntS: scIntS,

		scFracHi: scFracHi,
		scFracLo: scFracLo,

		ovfInt:  ovfInt,
		ovfIntS: ovfIntS,

		ovfFracHi: ovfFracHi,
		ovfFracLo: ovfFracLo,

		buf: newEmbedderBuffer(modIn),
	}
}

// ScaleEmbed scales and embeds p to the output modulus.
func (s *ScaleEmbedder) ScaleEmbed(p *Poly) *Poly {
	pOut := NewPoly(p.Rank(), p.ModLen())
	s.ScaleEmbedTo(pOut, p)
	return pOut
}

// ScaleEmbedTo scales and embeds p to pOut.
func (s *ScaleEmbedder) ScaleEmbedTo(pOut, p *Poly) {
	switch {
	case p.isNTT || pOut.isNTT:
		panic("ScaleEmbedTo: cannot scale NTT polynomials")
	case p.ModLen() != len(s.modIn):
		panic("ScaleEmbedTo: len(p.modulus) != len(s.ModIn)")
	}

	s.ScaleEmbedVecTo(p.Coeffs, pOut.Coeffs)
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
	inLen, outLen := len(s.modIn), len(s.modOut)
	M := (len(v[0]) >> 3) << 3

	switch {
	case inLen != len(v):
		panic("ScaleTo: len(in) != len(s.modIn)")
	case outLen != len(vOut):
		panic("ScaleTo: len(out) != len(s.modOut)")
	}

	var hi, lo uint64
	fHi := (*[8]uint64)(unsafe.Pointer(&s.buf.fHi[0]))
	fLo := (*[8]uint64)(unsafe.Pointer(&s.buf.fLo[0]))
	f64 := (*[8]uint64)(unsafe.Pointer(&s.buf.f64[0]))

	for k := 0; k < M; k += 8 {
		clear(fHi[:])
		clear(fLo[:])
		for i := 0; i < inLen; i++ {
			wIn := (*[8]uint64)(unsafe.Pointer(&v[i][k]))
			bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[i][0]))

			compInv, compInvS, modIn := s.compInv[i], s.compInvS[i], s.modIn[i]
			invLo, invHi := s.invLo[i], s.invHi[i]

			bufIn[0] = num.SMul(wIn[0], compInv, compInvS, modIn)
			bufIn[1] = num.SMul(wIn[1], compInv, compInvS, modIn)
			bufIn[2] = num.SMul(wIn[2], compInv, compInvS, modIn)
			bufIn[3] = num.SMul(wIn[3], compInv, compInvS, modIn)

			bufIn[4] = num.SMul(wIn[4], compInv, compInvS, modIn)
			bufIn[5] = num.SMul(wIn[5], compInv, compInvS, modIn)
			bufIn[6] = num.SMul(wIn[6], compInv, compInvS, modIn)
			bufIn[7] = num.SMul(wIn[7], compInv, compInvS, modIn)

			hi, lo = mulAndFloor(bufIn[0], invHi, invLo)
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)

			hi, lo = mulAndFloor(bufIn[1], invHi, invLo)
			fHi[1], fLo[1] = add128(fHi[1], fLo[1], hi, lo)

			hi, lo = mulAndFloor(bufIn[2], invHi, invLo)
			fHi[2], fLo[2] = add128(fHi[2], fLo[2], hi, lo)

			hi, lo = mulAndFloor(bufIn[3], invHi, invLo)
			fHi[3], fLo[3] = add128(fHi[3], fLo[3], hi, lo)

			hi, lo = mulAndFloor(bufIn[4], invHi, invLo)
			fHi[4], fLo[4] = add128(fHi[4], fLo[4], hi, lo)

			hi, lo = mulAndFloor(bufIn[5], invHi, invLo)
			fHi[5], fLo[5] = add128(fHi[5], fLo[5], hi, lo)

			hi, lo = mulAndFloor(bufIn[6], invHi, invLo)
			fHi[6], fLo[6] = add128(fHi[6], fLo[6], hi, lo)

			hi, lo = mulAndFloor(bufIn[7], invHi, invLo)
			fHi[7], fLo[7] = add128(fHi[7], fLo[7], hi, lo)
		}

		f64[0] = roundTo64(fHi[0], fLo[0])
		f64[1] = roundTo64(fHi[1], fLo[1])
		f64[2] = roundTo64(fHi[2], fLo[2])
		f64[3] = roundTo64(fHi[3], fLo[3])

		f64[4] = roundTo64(fHi[4], fLo[4])
		f64[5] = roundTo64(fHi[5], fLo[5])
		f64[6] = roundTo64(fHi[6], fLo[6])
		f64[7] = roundTo64(fHi[7], fLo[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			intOv, intOvS := s.ovfInt[i], s.ovfIntS[i]
			modOut, modOutv := s.modOut[i], s.modOut[i].Value()

			wOut[0] = num.SMul(modOutv-f64[0], intOv, intOvS, modOut)
			wOut[1] = num.SMul(modOutv-f64[1], intOv, intOvS, modOut)
			wOut[2] = num.SMul(modOutv-f64[2], intOv, intOvS, modOut)
			wOut[3] = num.SMul(modOutv-f64[3], intOv, intOvS, modOut)

			wOut[4] = num.SMul(modOutv-f64[4], intOv, intOvS, modOut)
			wOut[5] = num.SMul(modOutv-f64[5], intOv, intOvS, modOut)
			wOut[6] = num.SMul(modOutv-f64[6], intOv, intOvS, modOut)
			wOut[7] = num.SMul(modOutv-f64[7], intOv, intOvS, modOut)

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

		clear(s.buf.fHi[:])
		clear(s.buf.fLo[:])
		for i := 0; i < inLen; i++ {
			bufIn := (*[8]uint64)(unsafe.Pointer(&s.buf.in[i][0]))

			modscFracHi, modscFracLo := s.scFracHi[i], s.scFracLo[i]

			hi, lo = mulAndFloor(bufIn[0], modscFracHi, modscFracLo)
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)

			hi, lo = mulAndFloor(bufIn[1], modscFracHi, modscFracLo)
			fHi[1], fLo[1] = add128(fHi[1], fLo[1], hi, lo)

			hi, lo = mulAndFloor(bufIn[2], modscFracHi, modscFracLo)
			fHi[2], fLo[2] = add128(fHi[2], fLo[2], hi, lo)

			hi, lo = mulAndFloor(bufIn[3], modscFracHi, modscFracLo)
			fHi[3], fLo[3] = add128(fHi[3], fLo[3], hi, lo)

			hi, lo = mulAndFloor(bufIn[4], modscFracHi, modscFracLo)
			fHi[4], fLo[4] = add128(fHi[4], fLo[4], hi, lo)

			hi, lo = mulAndFloor(bufIn[5], modscFracHi, modscFracLo)
			fHi[5], fLo[5] = add128(fHi[5], fLo[5], hi, lo)

			hi, lo = mulAndFloor(bufIn[6], modscFracHi, modscFracLo)
			fHi[6], fLo[6] = add128(fHi[6], fLo[6], hi, lo)

			hi, lo = mulAndFloor(bufIn[7], modscFracHi, modscFracLo)
			fHi[7], fLo[7] = add128(fHi[7], fLo[7], hi, lo)
		}

		hi, lo = mulAndFloor(f64[0], s.ovfFracHi, s.ovfFracLo)
		fHi[0], fLo[0] = sub128(fHi[0], fLo[0], hi, lo)

		hi, lo = mulAndFloor(f64[1], s.ovfFracHi, s.ovfFracLo)
		fHi[1], fLo[1] = sub128(fHi[1], fLo[1], hi, lo)

		hi, lo = mulAndFloor(f64[2], s.ovfFracHi, s.ovfFracLo)
		fHi[2], fLo[2] = sub128(fHi[2], fLo[2], hi, lo)

		hi, lo = mulAndFloor(f64[3], s.ovfFracHi, s.ovfFracLo)
		fHi[3], fLo[3] = sub128(fHi[3], fLo[3], hi, lo)

		hi, lo = mulAndFloor(f64[4], s.ovfFracHi, s.ovfFracLo)
		fHi[4], fLo[4] = sub128(fHi[4], fLo[4], hi, lo)

		hi, lo = mulAndFloor(f64[5], s.ovfFracHi, s.ovfFracLo)
		fHi[5], fLo[5] = sub128(fHi[5], fLo[5], hi, lo)

		hi, lo = mulAndFloor(f64[6], s.ovfFracHi, s.ovfFracLo)
		fHi[6], fLo[6] = sub128(fHi[6], fLo[6], hi, lo)

		hi, lo = mulAndFloor(f64[7], s.ovfFracHi, s.ovfFracLo)
		fHi[7], fLo[7] = sub128(fHi[7], fLo[7], hi, lo)

		fHi[0], fLo[0] = roundTo128Signed(fHi[0], fLo[0])
		fHi[1], fLo[1] = roundTo128Signed(fHi[1], fLo[1])
		fHi[2], fLo[2] = roundTo128Signed(fHi[2], fLo[2])
		fHi[3], fLo[3] = roundTo128Signed(fHi[3], fLo[3])

		fHi[4], fLo[4] = roundTo128Signed(fHi[4], fLo[4])
		fHi[5], fLo[5] = roundTo128Signed(fHi[5], fLo[5])
		fHi[6], fLo[6] = roundTo128Signed(fHi[6], fLo[6])
		fHi[7], fLo[7] = roundTo128Signed(fHi[7], fLo[7])

		for i := 0; i < outLen; i++ {
			wOut := (*[8]uint64)(unsafe.Pointer(&vOut[i][k]))

			modOut := s.modOut[i]

			wOut[0] = add64To128Signed(wOut[0], fHi[0], fLo[0], modOut)
			wOut[1] = add64To128Signed(wOut[1], fHi[1], fLo[1], modOut)
			wOut[2] = add64To128Signed(wOut[2], fHi[2], fLo[2], modOut)
			wOut[3] = add64To128Signed(wOut[3], fHi[3], fLo[3], modOut)

			wOut[4] = add64To128Signed(wOut[4], fHi[4], fLo[4], modOut)
			wOut[5] = add64To128Signed(wOut[5], fHi[5], fLo[5], modOut)
			wOut[6] = add64To128Signed(wOut[6], fHi[6], fLo[6], modOut)
			wOut[7] = add64To128Signed(wOut[7], fHi[7], fLo[7], modOut)
		}
	}

	for k := M; k < len(v[0]); k++ {
		fHi[0], fLo[0] = 0, 0
		for i := 0; i < inLen; i++ {
			s.buf.in[i][0] = num.SMul(v[i][k], s.compInv[i], s.compInvS[i], s.modIn[i])

			hi, lo = mulAndFloor(s.buf.in[i][0], s.invHi[i], s.invLo[i])
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)
		}

		f64[0] = roundTo64(fHi[0], fLo[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = num.SMul(s.modOut[i].Value()-f64[0], s.ovfInt[i], s.ovfIntS[i], s.modOut[i])
			for j := 0; j < inLen; j++ {
				vOut[i][k] = num.Add(vOut[i][k], num.SMul(s.buf.in[j][0], s.scInt[i][j], s.scIntS[i][j], s.modOut[i]), s.modOut[i])
			}
		}

		fHi[0], fLo[0] = 0, 0
		for i := 0; i < inLen; i++ {
			hi, lo = mulAndFloor(s.buf.in[i][0], s.scFracHi[i], s.scFracLo[i])
			fHi[0], fLo[0] = add128(fHi[0], fLo[0], hi, lo)
		}

		hi, lo = mulAndFloor(f64[0], s.ovfFracHi, s.ovfFracLo)
		fHi[0], fLo[0] = sub128(fHi[0], fLo[0], hi, lo)

		fHi[0], fLo[0] = roundTo128Signed(fHi[0], fLo[0])

		for i := 0; i < outLen; i++ {
			vOut[i][k] = add64To128Signed(vOut[i][k], fHi[0], fLo[0], s.modOut[i])
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

		invHi: s.invHi,
		invLo: s.invLo,

		scInt:  s.scInt,
		scIntS: s.scIntS,

		scFracHi: s.scFracHi,
		scFracLo: s.scFracLo,

		ovfInt:  s.ovfInt,
		ovfIntS: s.ovfIntS,

		ovfFracHi: s.ovfFracHi,
		ovfFracLo: s.ovfFracLo,

		buf: newEmbedderBuffer(s.modIn),
	}
}
