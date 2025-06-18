package poly

import (
	"fmt"
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/hienaa-org/hienaa/math/mod"
)

const (
	fixedPrec = 124
	floatPrec = 62
	roundMask = 1<<(fixedPrec-floatPrec) - 1
)

// MultAndFloor returns floor((x * (yHi * 2^64 + yLo)) / 2^float_prec).
// We assure that the output is in [0, 2^128).
func MultAndFloor(x, yHi, yLo uint64) (uint64, uint64) {
	var resHi, resLo, tmpHi, tmpLo uint64

	resHi, resLo = bits.Mul64(x, yHi)
	resHi <<= (64 - floatPrec)
	resHi += (resLo >> floatPrec)
	resLo <<= (64 - floatPrec)

	tmpHi, tmpLo = bits.Mul64(x, yLo)
	tmpLo >>= floatPrec
	tmpLo += (tmpHi << (64 - floatPrec))
	tmpHi >>= floatPrec

	resLo, carry := bits.Add64(resLo, tmpLo, 0)
	resHi, _ = bits.Add64(resHi, tmpHi, carry)

	return resHi, resLo
}

// RoundTo64 returns round((xHi * 2^64 + xLo) / 2^fixed_prec).
// We assure that the output is in [0, 2^64).
func RoundTo64(xHi, xLo uint64) uint64 {
	var res, carry uint64

	res = (xLo >> (fixedPrec - floatPrec)) + (xHi << (64 - fixedPrec + floatPrec))
	carry = (xLo & roundMask) >> (fixedPrec - floatPrec - 1)
	res, _ = bits.Add64(res, carry, 0)

	return res
}

// RoundTo128 returns round((xHi * 2^64 + xLo) / 2^fixed_prec).
// We assure that the output is in [0, 2^128).
func RoundTo128(xHi, xLo uint64) (uint64, uint64) {
	var resHi, resLo, carry uint64

	resLo = (xLo >> (fixedPrec - floatPrec)) + (xHi << (64 - fixedPrec + floatPrec))
	resHi = xHi >> (fixedPrec - floatPrec)
	carry = (xLo & roundMask) >> (fixedPrec - floatPrec - 1)

	resLo, carry = bits.Add64(resLo, carry, 0)
	resHi, _ = bits.Add64(resHi, carry, 0)

	return resHi, resLo
}

func RoundTo128Signed(xHi, xLo uint64) (uint64, uint64) {
	if xHi&0x8000000000000000 != 0 {
		resHi, resLo := RoundTo128(-xHi, -xLo)
		return -resHi, -resLo
	} else {
		return RoundTo128(xHi, xLo)
	}
}

// Embedder holds the precomputed constants for the fixed-point HPS style basis embedding algorithm.
type Embedder struct {
	// modIn is the input modulus of the Embedder.
	modIn []*mod.Modulus
	// modOut is the output modulus of the Embedder.
	modOut []*mod.Modulus

	// compInv is the modular inverse of the compliment of the input modulus limb.
	compInv []uint64
	// compInvS is the Shoup multiplication form of the modular inverse of the compliment of the input modulus limb.
	compInvS []uint64

	// comp is the compliment of the input modulus limb.
	comp [][]uint64
	// compS is the Shoup multiplication form of the compliment of the input modulus limb.
	compS [][]uint64

	// negMod is the negative of the input modulus modulo the output modulus limb.
	negMod []uint64
	// negModS is the Shoup multiplication form of the negative of the input modulus modulo the output modulus limb.
	negModS []uint64

	// invLo is the low 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invLo []uint64
	// invHi is the high 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invHi []uint64

	// index holds the index of the input modulus limb if it overlaps with the output modulus limb.
	// for example, if modOut[i] = modIn[j], then index[i] = j.
	// -1 if the input modulus limb does not overlap with the output modulus limb.
	index []int

	// vLo is a temporary buffer used for the high-precision fixed-point arithmetic.
	vLo [8]uint64
	// vHi is a temporary buffer used for the high-precision fixed-point arithmetic.
	vHi [8]uint64

	// buff is a temporary buffer used for the embedding process.
	buff [][8]uint64
}

// NewEmbedder creates a new Embedder.
// modIn is the input modulus of the Embedder.
// modOut is the output modulus of the Embedder.
func NewEmbedder(modIn []*mod.Modulus, modOut []*mod.Modulus) *Embedder {
	switch {
	case len(modIn) == 0:
		panic("NewEmbedder: modIn cannot be empty")
	case len(modOut) == 0:
		panic("NewEmbedder: modOut cannot be empty")
	}

	lenIn := len(modIn)
	lenOut := len(modOut)

	compInv := make([]uint64, lenIn)
	compInvS := make([]uint64, lenIn)

	comp := make([][]uint64, lenOut)
	compS := make([][]uint64, lenOut)
	for i := 0; i < lenOut; i++ {
		comp[i] = make([]uint64, lenIn)
		compS[i] = make([]uint64, lenIn)
	}

	negMod := make([]uint64, lenOut)
	negModS := make([]uint64, lenOut)

	invLo := make([]uint64, lenIn)
	invHi := make([]uint64, lenIn)

	index := make([]int, lenOut)

	vLo := [8]uint64{}
	vHi := [8]uint64{}

	buff := make([][8]uint64, lenIn)
	for i := 0; i < lenIn; i++ {
		buff[i] = [8]uint64{}
	}

	tmpInt := big.NewInt(0)
	tmpRat := big.NewRat(0, 1)

	for i := 0; i < lenIn; i++ {
		compInv[i] = 1
		for j := 0; j < lenIn; j++ {
			if i != j {
				compInv[i] = mod.Mul(compInv[i], mod.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
			compInvS[i] = mod.SForm(compInv[i], modIn[i])
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

	for i := 0; i < lenOut; i++ {
		negMod[i] = modOut[i].Value() - 1
		index[i] = -1

		for j := 0; j < lenIn; j++ {
			comp[i][j] = 1
			for k := 0; k < lenIn; k++ {
				if j != k {
					comp[i][j] = mod.Mul(comp[i][j], modIn[k].Value(), modOut[i])
				}
			}
			compS[i][j] = mod.SForm(comp[i][j], modOut[i])

			negMod[i] = mod.Mul(negMod[i], modIn[j].Value(), modOut[i])

			if modOut[i].Value() == modIn[j].Value() {
				index[i] = j
			}
		}

		negModS[i] = mod.SForm(negMod[i], modOut[i])
	}

	return &Embedder{
		modIn:    modIn,
		modOut:   modOut,
		compInv:  compInv,
		compInvS: compInvS,
		comp:     comp,
		compS:    compS,
		negMod:   negMod,
		negModS:  negModS,
		invLo:    invLo,
		invHi:    invHi,
		index:    index,
		vLo:      vLo,
		vHi:      vHi,
		buff:     buff,
	}
}

// EmbedTo embeds the input coefficient vector into the output coefficient vector.
// All vectors are in signed representation.
func (e *Embedder) EmbedTo(in, out [][]uint64) {
	degree, lenIn, lenOut := len(in[0]), len(in), len(out)
	M := (degree >> 3) << 3

	lenIn, lenOut = min(lenIn, len(e.modIn)), min(lenOut, len(e.modOut))

	if lenIn == 1 {
		Q := e.modIn[0].Value()
		halfQ := Q >> 1

		for i := 0; i < lenOut; i++ {
			if e.index[i] == 0 {
				copy(out[i], in[0])
			} else {
				outi, in0, modOuti := out[i], in[0], e.modOut[i]

				for j := 0; j < M; j += 8 {
					if in0[j] <= halfQ {
						outi[j] = mod.Reduce(in0[j], modOuti)
					} else {
						outi[j] = mod.Reduce(Q-in0[j], modOuti)
					}

					if in0[j+1] <= halfQ {
						outi[j+1] = mod.Reduce(in0[j+1], modOuti)
					} else {
						outi[j+1] = mod.Reduce(Q-in0[j+1], modOuti)
					}

					if in0[j+2] <= halfQ {
						outi[j+2] = mod.Reduce(in0[j+2], modOuti)
					} else {
						outi[j+2] = mod.Reduce(Q-in0[j+2], modOuti)
					}

					if in0[j+3] <= halfQ {
						outi[j+3] = mod.Reduce(in0[j+3], modOuti)
					} else {
						outi[j+3] = mod.Reduce(Q-in0[j+3], modOuti)
					}

					if in0[j+4] <= halfQ {
						outi[j+4] = mod.Reduce(in0[j+4], modOuti)
					} else {
						outi[j+4] = mod.Reduce(Q-in0[j+4], modOuti)
					}

					if in0[j+5] <= halfQ {
						outi[j+5] = mod.Reduce(in0[j+5], modOuti)
					} else {
						outi[j+5] = mod.Reduce(Q-in0[j+5], modOuti)
					}

					if in0[j+6] <= halfQ {
						outi[j+6] = mod.Reduce(in0[j+6], modOuti)
					} else {
						outi[j+6] = mod.Reduce(Q-in0[j+6], modOuti)
					}

					if in0[j+7] <= halfQ {
						outi[j+7] = mod.Reduce(in0[j+7], modOuti)
					} else {
						outi[j+7] = mod.Reduce(Q-in0[j+7], modOuti)
					}
				}

				for j := M; j < degree; j++ {
					if in0[j] <= halfQ {
						outi[j] = mod.Reduce(in0[j], modOuti)
					} else {
						outi[j] = mod.Reduce(Q-in0[j], modOuti)
					}
				}
			}
		}
	} else {
		var tmpLo, tmpHi, carry uint64

		for k := 0; k < M; k += 8 {
			clear(e.vLo[:])
			clear(e.vHi[:])

			for i := 0; i < lenIn; i++ {
				ini := (*[8]uint64)(unsafe.Pointer(&in[i][k]))
				compInvi, compInvSi, modIni := e.compInv[i], e.compInvS[i], e.modIn[i]
				invLoi, invHii := e.invLo[i], e.invHi[i]

				e.buff[i][0] = mod.SMul(ini[0], compInvi, compInvSi, modIni)
				e.buff[i][1] = mod.SMul(ini[1], compInvi, compInvSi, modIni)
				e.buff[i][2] = mod.SMul(ini[2], compInvi, compInvSi, modIni)
				e.buff[i][3] = mod.SMul(ini[3], compInvi, compInvSi, modIni)
				e.buff[i][4] = mod.SMul(ini[4], compInvi, compInvSi, modIni)
				e.buff[i][5] = mod.SMul(ini[5], compInvi, compInvSi, modIni)
				e.buff[i][6] = mod.SMul(ini[6], compInvi, compInvSi, modIni)
				e.buff[i][7] = mod.SMul(ini[7], compInvi, compInvSi, modIni)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][0], invHii, invLoi)
				e.vLo[0], carry = bits.Add64(e.vLo[0], tmpLo, 0)
				e.vHi[0], _ = bits.Add64(e.vHi[0], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][1], invHii, invLoi)
				e.vLo[1], carry = bits.Add64(e.vLo[1], tmpLo, 0)
				e.vHi[1], _ = bits.Add64(e.vHi[1], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][2], invHii, invLoi)
				e.vLo[2], carry = bits.Add64(e.vLo[2], tmpLo, 0)
				e.vHi[2], _ = bits.Add64(e.vHi[2], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][3], invHii, invLoi)
				e.vLo[3], carry = bits.Add64(e.vLo[3], tmpLo, 0)
				e.vHi[3], _ = bits.Add64(e.vHi[3], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][4], invHii, invLoi)
				e.vLo[4], carry = bits.Add64(e.vLo[4], tmpLo, 0)
				e.vHi[4], _ = bits.Add64(e.vHi[4], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][5], invHii, invLoi)
				e.vLo[5], carry = bits.Add64(e.vLo[5], tmpLo, 0)
				e.vHi[5], _ = bits.Add64(e.vHi[5], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][6], invHii, invLoi)
				e.vLo[6], carry = bits.Add64(e.vLo[6], tmpLo, 0)
				e.vHi[6], _ = bits.Add64(e.vHi[6], tmpHi, carry)

				tmpHi, tmpLo = MultAndFloor(e.buff[i][7], invHii, invLoi)
				e.vLo[7], carry = bits.Add64(e.vLo[7], tmpLo, 0)
				e.vHi[7], _ = bits.Add64(e.vHi[7], tmpHi, carry)
			}

			e.vLo[0] = RoundTo64(e.vHi[0], e.vLo[0])
			e.vLo[1] = RoundTo64(e.vHi[1], e.vLo[1])
			e.vLo[2] = RoundTo64(e.vHi[2], e.vLo[2])
			e.vLo[3] = RoundTo64(e.vHi[3], e.vLo[3])
			e.vLo[4] = RoundTo64(e.vHi[4], e.vLo[4])
			e.vLo[5] = RoundTo64(e.vHi[5], e.vLo[5])
			e.vLo[6] = RoundTo64(e.vHi[6], e.vLo[6])
			e.vLo[7] = RoundTo64(e.vHi[7], e.vLo[7])

			for i := 0; i < lenOut; i++ {
				outi := (*[8]uint64)(unsafe.Pointer(&out[i][k]))

				if 0 <= e.index[i] && e.index[i] < lenIn {
					ini := (*[8]uint64)(unsafe.Pointer(&in[e.index[i]][k]))
					copy(outi[:], ini[:])
				} else {
					compi, compSi := e.comp[i], e.compS[i]
					negModi, negModSi := e.negMod[i], e.negModS[i]
					modOuti := e.modOut[i]

					outi[0] = mod.SMul(e.vLo[0], negModi, negModSi, modOuti)
					outi[1] = mod.SMul(e.vLo[1], negModi, negModSi, modOuti)
					outi[2] = mod.SMul(e.vLo[2], negModi, negModSi, modOuti)
					outi[3] = mod.SMul(e.vLo[3], negModi, negModSi, modOuti)
					outi[4] = mod.SMul(e.vLo[4], negModi, negModSi, modOuti)
					outi[5] = mod.SMul(e.vLo[5], negModi, negModSi, modOuti)
					outi[6] = mod.SMul(e.vLo[6], negModi, negModSi, modOuti)
					outi[7] = mod.SMul(e.vLo[7], negModi, negModSi, modOuti)

					for j := 0; j < lenIn; j++ {
						buffj := e.buff[j]
						compij, compSij := compi[j], compSi[j]

						outi[0] = mod.Add(outi[0], mod.SMul(buffj[0], compij, compSij, modOuti), modOuti)
						outi[1] = mod.Add(outi[1], mod.SMul(buffj[1], compij, compSij, modOuti), modOuti)
						outi[2] = mod.Add(outi[2], mod.SMul(buffj[2], compij, compSij, modOuti), modOuti)
						outi[3] = mod.Add(outi[3], mod.SMul(buffj[3], compij, compSij, modOuti), modOuti)
						outi[4] = mod.Add(outi[4], mod.SMul(buffj[4], compij, compSij, modOuti), modOuti)
						outi[5] = mod.Add(outi[5], mod.SMul(buffj[5], compij, compSij, modOuti), modOuti)
						outi[6] = mod.Add(outi[6], mod.SMul(buffj[6], compij, compSij, modOuti), modOuti)
						outi[7] = mod.Add(outi[7], mod.SMul(buffj[7], compij, compSij, modOuti), modOuti)
					}
				}
			}
		}

		for k := M; k < degree; k++ {
			for i := 0; i < lenIn; i++ {
				e.buff[i][0] = mod.SMul(in[i][k], e.compInv[i], e.compInvS[i], e.modIn[i])
			}

			vHi0, vLo0 := uint64(0), uint64(0)
			for i := 0; i < lenIn; i++ {
				tmpHi, tmpLo = MultAndFloor(e.buff[i][0], e.invHi[i], e.invLo[i])
				vLo0, carry = bits.Add64(vLo0, tmpLo, 0)
				vHi0, _ = bits.Add64(vHi0, tmpHi, carry)
			}
			vLo0 = RoundTo64(vHi0, vLo0)

			for i := 0; i < lenOut; i++ {
				if 0 <= e.index[i] && e.index[i] < lenIn {
					out[i][k] = in[e.index[i]][k]
				} else {
					out[i][k] = mod.SMul(vLo0, e.negMod[i], e.negModS[i], e.modOut[i])

					for j := 0; j < lenIn; j++ {
						out[i][k] = mod.Add(out[i][k], mod.SMul(e.buff[j][0], e.comp[i][j], e.compS[i][j], e.modOut[i]), e.modOut[i])
					}
				}
			}
		}
	}
}

// Scaler holds the precomputed constants for the fixed-point HPS style scaling algorithm.
type Scaler struct {
	// modIn is the input modulus of the Scaler.
	modIn []*mod.Modulus
	// modOut is the output modulus of the Scaler.
	modOut []*mod.Modulus

	// compInv is the modular inverse of the compliment of the input modulus limb.
	compInv []uint64
	// compInvS is the Shoup multiplication form of the modular inverse of the compliment of the input modulus limb.
	compInvS []uint64

	// intSc is the integer part of the scaling factor.
	intSc [][]uint64
	// intScS is the Shoup multiplication form of the integer part of the scaling factor.
	intScS [][]uint64

	// decScLo is the low 64 bits of the fixed-point approximation of the decimal part of the scaling factor.
	decScLo []uint64
	// decScHi is the high 64 bits of the fixed-point approximation of the decimal part of the scaling factor.
	decScHi []uint64

	// vLo is a temporary buffer used for the high-precision fixed-point arithmetic.
	vLo [8]uint64
	// vHi is a temporary buffer used for the high-precision fixed-point arithmetic.
	vHi [8]uint64

	// buff is a temporary buffer used for the scaling process.
	buff [][8]uint64
}

func NewScaler(modIn []*mod.Modulus, modOut []*mod.Modulus) *Scaler {
	switch {
	case len(modIn) == 0:
		panic("NewScaler: modIn cannot be empty")
	case len(modOut) == 0:
		panic("NewScaler: modOut cannot be empty")
	}

	lenIn := len(modIn)
	lenOut := len(modOut)

	compInv := make([]uint64, lenIn)
	compInvS := make([]uint64, lenIn)

	intSc := make([][]uint64, lenOut)
	intScS := make([][]uint64, lenOut)
	for i := 0; i < lenOut; i++ {
		intSc[i] = make([]uint64, lenIn)
		intScS[i] = make([]uint64, lenIn)
	}

	decScLo := make([]uint64, lenIn)
	decScHi := make([]uint64, lenIn)

	vLo := [8]uint64{}
	vHi := [8]uint64{}

	buff := make([][8]uint64, lenIn)
	for i := 0; i < lenIn; i++ {
		buff[i] = [8]uint64{}
	}

	for i := 0; i < lenIn; i++ {
		compInv[i] = 1
		for j := 0; j < lenIn; j++ {
			if i != j {
				compInv[i] = mod.Mul(compInv[i], mod.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
		}
		compInvS[i] = mod.SForm(compInv[i], modIn[i])
	}

	Pbig := big.NewInt(1)
	tmpInt := big.NewInt(0)
	tmpRat := big.NewRat(0, 1)
	intScBig := big.NewInt(0)

	for i := 0; i < lenOut; i++ {
		Pbig.Mul(Pbig, tmpInt.SetUint64(modOut[i].Value()))
	}

	for i := 0; i < lenIn; i++ {
		intScBig.Div(Pbig, tmpInt.SetUint64(modIn[i].Value()))
		for j := 0; j < lenOut; j++ {
			tmpInt.Mod(intScBig, tmpInt.SetUint64(modOut[j].Value()))
			intSc[j][i] = tmpInt.Uint64()
			intScS[j][i] = mod.SForm(intSc[j][i], modOut[j])
		}

		tmpInt.Mul(intScBig, tmpInt.SetUint64(modIn[i].Value()))
		tmpInt.Sub(Pbig, tmpInt)

		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Num().Set(tmpInt.Lsh(tmpInt, fixedPrec-64))

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		decScHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		decScLo[i] = tmpInt.Uint64()
	}

	return &Scaler{
		modIn:    modIn,
		modOut:   modOut,
		compInv:  compInv,
		compInvS: compInvS,
		intSc:    intSc,
		intScS:   intScS,
		decScLo:  decScLo,
		decScHi:  decScHi,
		vLo:      vLo,
		vHi:      vHi,
		buff:     buff,
	}
}

// ScaleTo scales the input polynomial to the output polynomial.
func (s *Scaler) ScaleTo(in, out [][]uint64) {
	degree, lenIn, lenOut := len(in[0]), len(s.modIn), len(s.modOut)
	M := (degree >> 3) << 3

	switch {
	case lenIn != len(in):
		panic("ScaleTo: len(in) != len(s.modIn)")
	case lenOut != len(out):
		panic("ScaleTo: len(out) != len(s.modOut)")
	}

	var tmpLo, tmpHi, carry uint64

	for k := 0; k < M; k += 8 {
		clear(s.vLo[:])
		clear(s.vHi[:])

		for i := 0; i < lenIn; i++ {
			ini := (*[8]uint64)(unsafe.Pointer(&in[i][k]))
			compInv, compInvS, modIn := s.compInv[i], s.compInvS[i], s.modIn[i]
			decScHi, decScLo := s.decScHi[i], s.decScLo[i]

			s.buff[i][0] = mod.SMul(ini[0], compInv, compInvS, modIn)
			s.buff[i][1] = mod.SMul(ini[1], compInv, compInvS, modIn)
			s.buff[i][2] = mod.SMul(ini[2], compInv, compInvS, modIn)
			s.buff[i][3] = mod.SMul(ini[3], compInv, compInvS, modIn)
			s.buff[i][4] = mod.SMul(ini[4], compInv, compInvS, modIn)
			s.buff[i][5] = mod.SMul(ini[5], compInv, compInvS, modIn)
			s.buff[i][6] = mod.SMul(ini[6], compInv, compInvS, modIn)
			s.buff[i][7] = mod.SMul(ini[7], compInv, compInvS, modIn)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][0], decScHi, decScLo)
			s.vLo[0], carry = bits.Add64(s.vLo[0], tmpLo, 0)
			s.vHi[0], _ = bits.Add64(s.vHi[0], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][1], decScHi, decScLo)
			s.vLo[1], carry = bits.Add64(s.vLo[1], tmpLo, 0)
			s.vHi[1], _ = bits.Add64(s.vHi[1], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][2], decScHi, decScLo)
			s.vLo[2], carry = bits.Add64(s.vLo[2], tmpLo, 0)
			s.vHi[2], _ = bits.Add64(s.vHi[2], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][3], decScHi, decScLo)
			s.vLo[3], carry = bits.Add64(s.vLo[3], tmpLo, 0)
			s.vHi[3], _ = bits.Add64(s.vHi[3], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][4], decScHi, decScLo)
			s.vLo[4], carry = bits.Add64(s.vLo[4], tmpLo, 0)
			s.vHi[4], _ = bits.Add64(s.vHi[4], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][5], decScHi, decScLo)
			s.vLo[5], carry = bits.Add64(s.vLo[5], tmpLo, 0)
			s.vHi[5], _ = bits.Add64(s.vHi[5], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][6], decScHi, decScLo)
			s.vLo[6], carry = bits.Add64(s.vLo[6], tmpLo, 0)
			s.vHi[6], _ = bits.Add64(s.vHi[6], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][7], decScHi, decScLo)
			s.vLo[7], carry = bits.Add64(s.vLo[7], tmpLo, 0)
			s.vHi[7], _ = bits.Add64(s.vHi[7], tmpHi, carry)
		}

		s.vHi[0], s.vLo[0] = RoundTo128(s.vHi[0], s.vLo[0])
		s.vHi[1], s.vLo[1] = RoundTo128(s.vHi[1], s.vLo[1])
		s.vHi[2], s.vLo[2] = RoundTo128(s.vHi[2], s.vLo[2])
		s.vHi[3], s.vLo[3] = RoundTo128(s.vHi[3], s.vLo[3])
		s.vHi[4], s.vLo[4] = RoundTo128(s.vHi[4], s.vLo[4])
		s.vHi[5], s.vLo[5] = RoundTo128(s.vHi[5], s.vLo[5])
		s.vHi[6], s.vLo[6] = RoundTo128(s.vHi[6], s.vLo[6])
		s.vHi[7], s.vLo[7] = RoundTo128(s.vHi[7], s.vLo[7])

		for i := 0; i < lenOut; i++ {
			outi := (*[8]uint64)(unsafe.Pointer(&out[i][k]))
			outi[0] = mod.Reduce128(s.vHi[0], s.vLo[0], s.modOut[i])
			outi[1] = mod.Reduce128(s.vHi[1], s.vLo[1], s.modOut[i])
			outi[2] = mod.Reduce128(s.vHi[2], s.vLo[2], s.modOut[i])
			outi[3] = mod.Reduce128(s.vHi[3], s.vLo[3], s.modOut[i])
			outi[4] = mod.Reduce128(s.vHi[4], s.vLo[4], s.modOut[i])
			outi[5] = mod.Reduce128(s.vHi[5], s.vLo[5], s.modOut[i])
			outi[6] = mod.Reduce128(s.vHi[6], s.vLo[6], s.modOut[i])
			outi[7] = mod.Reduce128(s.vHi[7], s.vLo[7], s.modOut[i])

			for j := 0; j < lenIn; j++ {
				outi[0] = mod.Add(outi[0], mod.SMul(s.buff[j][0], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[1] = mod.Add(outi[1], mod.SMul(s.buff[j][1], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[2] = mod.Add(outi[2], mod.SMul(s.buff[j][2], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[3] = mod.Add(outi[3], mod.SMul(s.buff[j][3], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[4] = mod.Add(outi[4], mod.SMul(s.buff[j][4], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[5] = mod.Add(outi[5], mod.SMul(s.buff[j][5], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[6] = mod.Add(outi[6], mod.SMul(s.buff[j][6], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
				outi[7] = mod.Add(outi[7], mod.SMul(s.buff[j][7], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
			}
		}
	}

	for k := M; k < degree; k++ {
		vLo0, vHi0 := uint64(0), uint64(0)
		for i := 0; i < lenIn; i++ {
			s.buff[i][0] = mod.SMul(in[i][k], s.compInv[i], s.compInvS[i], s.modIn[i])
			tmpHi, tmpLo = MultAndFloor(s.buff[i][0], s.decScHi[i], s.decScLo[i])
			vLo0, carry = bits.Add64(vLo0, tmpLo, 0)
			vHi0, _ = bits.Add64(vHi0, tmpHi, carry)
		}
		vHi0, vLo0 = RoundTo128(vHi0, vLo0)

		for i := 0; i < lenOut; i++ {
			out[i][k] = mod.Reduce128(vHi0, vLo0, s.modOut[i])
			for j := 0; j < lenIn; j++ {
				out[i][k] = mod.Add(out[i][k], mod.SMul(s.buff[j][0], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
			}
		}
	}
}

type ScaleEmbedder struct {
	// modIn is the input modulus of the ScaleEmbedder.
	modIn []*mod.Modulus
	// modOut is the output modulus of the ScaleEmbedder.
	modOut []*mod.Modulus

	// compInv is the modular inverse of the compliment of the input modulus limb.
	compInv []uint64
	// compInvS is the Shoup multiplication form of the modular inverse of the compliment of the input modulus limb.
	compInvS []uint64

	// invLo is the low 64 bits of the fixed-point approximation of the inverse of the input modulus limb.
	invLo []uint64
	invHi []uint64

	intSc  [][]uint64
	intScS [][]uint64

	decScLo []uint64
	decScHi []uint64

	intOv   []uint64
	intOvS  []uint64
	decOvLo uint64
	decOvHi uint64

	vLo [8]uint64
	vHi [8]uint64
	v   [8]uint64

	buff [][8]uint64
}

func NewScaleEmbedder(modIn []*mod.Modulus, modOut []*mod.Modulus, scale *big.Rat) *ScaleEmbedder {
	lenIn, lenOut := len(modIn), len(modOut)

	compInv := make([]uint64, lenIn)
	compInvS := make([]uint64, lenIn)

	invLo := make([]uint64, lenIn)
	invHi := make([]uint64, lenIn)

	intSc := make([][]uint64, lenOut)
	intScS := make([][]uint64, lenOut)
	for i := 0; i < lenOut; i++ {
		intSc[i] = make([]uint64, lenIn)
		intScS[i] = make([]uint64, lenIn)
	}

	decScLo := make([]uint64, lenIn)
	decScHi := make([]uint64, lenIn)

	intOv := make([]uint64, lenOut)
	intOvS := make([]uint64, lenOut)

	vLo := [8]uint64{}
	vHi := [8]uint64{}
	v := [8]uint64{}

	buff := make([][8]uint64, lenIn)
	for i := 0; i < lenIn; i++ {
		buff[i] = [8]uint64{}
	}

	// Compute the modular inverse of the input modulus limbs
	for i := 0; i < lenIn; i++ {
		compInv[i] = 1
		for j := 0; j < lenIn; j++ {
			if i != j {
				compInv[i] = mod.Mul(compInv[i], mod.Inv(modIn[j].Value(), modIn[i]), modIn[i])
			}
		}
		compInvS[i] = mod.SForm(compInv[i], modIn[i])
	}

	Qbig := big.NewInt(1)
	tmpInt := big.NewInt(0)
	tmpRat := big.NewRat(0, 1)
	intScBig := big.NewInt(0)

	for i := 0; i < lenIn; i++ {
		Qbig.Mul(Qbig, tmpInt.SetUint64(modIn[i].Value()))

		// Compute a fixed-point approximation of 1/modIn[i]
		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Num().Lsh(tmpRat.Num().SetInt64(1), fixedPrec-64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		invHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)

		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		invLo[i] = tmpInt.Uint64()
	}

	// Compute the fixed-point approximation of Q/Qi * scale,
	// where Q is the input modulus and Qi is the i-th input modulus limb.
	for i := 0; i < lenIn; i++ {
		tmpRat.Num().Set(Qbig)
		tmpRat.Denom().SetUint64(modIn[i].Value())
		tmpRat.Mul(tmpRat, scale)

		intScBig.Div(tmpRat.Num(), tmpRat.Denom())
		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(intScBig, tmpRat.Denom()))

		for j := 0; j < lenOut; j++ {
			tmpInt.Mod(intScBig, tmpInt.SetUint64(modOut[j].Value()))
			intSc[j][i] = tmpInt.Uint64()
			intScS[j][i] = mod.SForm(intSc[j][i], modOut[j])
		}

		tmpRat.Num().Lsh(tmpRat.Num(), fixedPrec-64)
		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		decScHi[i] = tmpInt.Uint64()

		tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
		tmpRat.Num().Lsh(tmpRat.Num(), 64)
		tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
		decScLo[i] = tmpInt.Uint64()
	}

	intOvBig := big.NewInt(0)

	tmpRat.Set(scale)
	tmpRat.Num().Mul(tmpRat.Num(), Qbig)
	intOvBig.Div(tmpRat.Num(), tmpRat.Denom())
	for i := 0; i < lenOut; i++ {
		tmpInt.Mod(intOvBig, tmpInt.SetUint64(modOut[i].Value()))
		intOv[i] = tmpInt.Uint64()
		intOvS[i] = mod.SForm(intOv[i], modOut[i])
	}

	tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(intOvBig, tmpRat.Denom()))
	tmpRat.Num().Lsh(tmpRat.Num(), fixedPrec-64)
	tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
	decOvHi := tmpInt.Uint64()

	tmpRat.Num().Sub(tmpRat.Num(), tmpInt.Mul(tmpInt, tmpRat.Denom()))
	tmpRat.Num().Lsh(tmpRat.Num(), 64)
	tmpInt.Div(tmpRat.Num(), tmpRat.Denom())
	decOvLo := tmpInt.Uint64()

	fmt.Println(decScLo)
	fmt.Println(decScHi)

	return &ScaleEmbedder{
		modIn:    modIn,
		modOut:   modOut,
		compInv:  compInv,
		compInvS: compInvS,
		invLo:    invLo,
		invHi:    invHi,
		intSc:    intSc,
		intScS:   intScS,
		decScLo:  decScLo,
		decScHi:  decScHi,
		intOv:    intOv,
		intOvS:   intOvS,
		decOvLo:  decOvLo,
		decOvHi:  decOvHi,
		vLo:      vLo,
		vHi:      vHi,
		v:        v,
		buff:     buff,
	}
}

func (s *ScaleEmbedder) ScaleEmbedTo(in, out [][]uint64) {
	degree, lenIn, lenOut := len(in[0]), len(s.modIn), len(s.modOut)
	M := (degree >> 3) << 3

	switch {
	case lenIn != len(in):
		panic("ScaleTo: len(in) != len(s.modIn)")
	case lenOut != len(out):
		panic("ScaleTo: len(out) != len(s.modOut)")
	}

	var tmpLo, tmpHi, carry, borrow uint64

	for k := 0; k < M; k += 8 {
		clear(s.vLo[:])
		clear(s.vHi[:])

		for i := 0; i < lenIn; i++ {
			ini := (*[8]uint64)(unsafe.Pointer(&in[i][k]))
			compInv, compInvS, modIn := s.compInv[i], s.compInvS[i], s.modIn[i]
			invLo, invHi := s.invLo[i], s.invHi[i]

			s.buff[i][0] = mod.SMul(ini[0], compInv, compInvS, modIn)
			s.buff[i][1] = mod.SMul(ini[1], compInv, compInvS, modIn)
			s.buff[i][2] = mod.SMul(ini[2], compInv, compInvS, modIn)
			s.buff[i][3] = mod.SMul(ini[3], compInv, compInvS, modIn)
			s.buff[i][4] = mod.SMul(ini[4], compInv, compInvS, modIn)
			s.buff[i][5] = mod.SMul(ini[5], compInv, compInvS, modIn)
			s.buff[i][6] = mod.SMul(ini[6], compInv, compInvS, modIn)
			s.buff[i][7] = mod.SMul(ini[7], compInv, compInvS, modIn)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][0], invHi, invLo)
			s.vLo[0], carry = bits.Add64(s.vLo[0], tmpLo, 0)
			s.vHi[0], _ = bits.Add64(s.vHi[0], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][1], invHi, invLo)
			s.vLo[1], carry = bits.Add64(s.vLo[1], tmpLo, 0)
			s.vHi[1], _ = bits.Add64(s.vHi[1], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][2], invHi, invLo)
			s.vLo[2], carry = bits.Add64(s.vLo[2], tmpLo, 0)
			s.vHi[2], _ = bits.Add64(s.vHi[2], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][3], invHi, invLo)
			s.vLo[3], carry = bits.Add64(s.vLo[3], tmpLo, 0)
			s.vHi[3], _ = bits.Add64(s.vHi[3], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][4], invHi, invLo)
			s.vLo[4], carry = bits.Add64(s.vLo[4], tmpLo, 0)
			s.vHi[4], _ = bits.Add64(s.vHi[4], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][5], invHi, invLo)
			s.vLo[5], carry = bits.Add64(s.vLo[5], tmpLo, 0)
			s.vHi[5], _ = bits.Add64(s.vHi[5], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][6], invHi, invLo)
			s.vLo[6], carry = bits.Add64(s.vLo[6], tmpLo, 0)
			s.vHi[6], _ = bits.Add64(s.vHi[6], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][7], invHi, invLo)
			s.vLo[7], carry = bits.Add64(s.vLo[7], tmpLo, 0)
			s.vHi[7], _ = bits.Add64(s.vHi[7], tmpHi, carry)
		}

		s.v[0] = RoundTo64(s.vHi[0], s.vLo[0])
		s.v[1] = RoundTo64(s.vHi[1], s.vLo[1])
		s.v[2] = RoundTo64(s.vHi[2], s.vLo[2])
		s.v[3] = RoundTo64(s.vHi[3], s.vLo[3])
		s.v[4] = RoundTo64(s.vHi[4], s.vLo[4])
		s.v[5] = RoundTo64(s.vHi[5], s.vLo[5])
		s.v[6] = RoundTo64(s.vHi[6], s.vLo[6])
		s.v[7] = RoundTo64(s.vHi[7], s.vLo[7])

		for i := 0; i < lenOut; i++ {
			outi := (*[8]uint64)(unsafe.Pointer(&out[i][k]))
			intOv, intOvS, modOut := s.intOv[i], s.intOvS[i], s.modOut[i]
			intSc, intScS := s.intSc[i], s.intScS[i]

			outi[0] = mod.SMul(modOut.Value()-s.v[0], intOv, intOvS, modOut)
			outi[1] = mod.SMul(modOut.Value()-s.v[1], intOv, intOvS, modOut)
			outi[2] = mod.SMul(modOut.Value()-s.v[2], intOv, intOvS, modOut)
			outi[3] = mod.SMul(modOut.Value()-s.v[3], intOv, intOvS, modOut)
			outi[4] = mod.SMul(modOut.Value()-s.v[4], intOv, intOvS, modOut)
			outi[5] = mod.SMul(modOut.Value()-s.v[5], intOv, intOvS, modOut)
			outi[6] = mod.SMul(modOut.Value()-s.v[6], intOv, intOvS, modOut)
			outi[7] = mod.SMul(modOut.Value()-s.v[7], intOv, intOvS, modOut)

			for j := 0; j < lenIn; j++ {
				intScj, intScSj := intSc[j], intScS[j]

				outi[0] = mod.Add(outi[0], mod.SMul(s.buff[j][0], intScj, intScSj, modOut), modOut)
				outi[1] = mod.Add(outi[1], mod.SMul(s.buff[j][1], intScj, intScSj, modOut), modOut)
				outi[2] = mod.Add(outi[2], mod.SMul(s.buff[j][2], intScj, intScSj, modOut), modOut)
				outi[3] = mod.Add(outi[3], mod.SMul(s.buff[j][3], intScj, intScSj, modOut), modOut)
				outi[4] = mod.Add(outi[4], mod.SMul(s.buff[j][4], intScj, intScSj, modOut), modOut)
				outi[5] = mod.Add(outi[5], mod.SMul(s.buff[j][5], intScj, intScSj, modOut), modOut)
				outi[6] = mod.Add(outi[6], mod.SMul(s.buff[j][6], intScj, intScSj, modOut), modOut)
				outi[7] = mod.Add(outi[7], mod.SMul(s.buff[j][7], intScj, intScSj, modOut), modOut)
			}
		}

		clear(s.vLo[:])
		clear(s.vHi[:])
		for i := 0; i < lenIn; i++ {
			decScHi, decScLo := s.decScHi[i], s.decScLo[i]

			tmpHi, tmpLo = MultAndFloor(s.buff[i][0], decScHi, decScLo)
			s.vLo[0], carry = bits.Add64(s.vLo[0], tmpLo, 0)
			s.vHi[0], _ = bits.Add64(s.vHi[0], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][1], decScHi, decScLo)
			s.vLo[1], carry = bits.Add64(s.vLo[1], tmpLo, 0)
			s.vHi[1], _ = bits.Add64(s.vHi[1], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][2], decScHi, decScLo)
			s.vLo[2], carry = bits.Add64(s.vLo[2], tmpLo, 0)
			s.vHi[2], _ = bits.Add64(s.vHi[2], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][3], decScHi, decScLo)
			s.vLo[3], carry = bits.Add64(s.vLo[3], tmpLo, 0)
			s.vHi[3], _ = bits.Add64(s.vHi[3], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][4], decScHi, decScLo)
			s.vLo[4], carry = bits.Add64(s.vLo[4], tmpLo, 0)
			s.vHi[4], _ = bits.Add64(s.vHi[4], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][5], decScHi, decScLo)
			s.vLo[5], carry = bits.Add64(s.vLo[5], tmpLo, 0)
			s.vHi[5], _ = bits.Add64(s.vHi[5], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][6], decScHi, decScLo)
			s.vLo[6], carry = bits.Add64(s.vLo[6], tmpLo, 0)
			s.vHi[6], _ = bits.Add64(s.vHi[6], tmpHi, carry)

			tmpHi, tmpLo = MultAndFloor(s.buff[i][7], decScHi, decScLo)
			s.vLo[7], carry = bits.Add64(s.vLo[7], tmpLo, 0)
			s.vHi[7], _ = bits.Add64(s.vHi[7], tmpHi, carry)
		}

		decOvHi, decOvLo := s.decOvHi, s.decOvLo

		tmpHi, tmpLo = MultAndFloor(s.v[0], decOvHi, decOvLo)
		s.vLo[0], borrow = bits.Sub64(s.vLo[0], tmpLo, 0)
		s.vHi[0], _ = bits.Sub64(s.vHi[0], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[1], decOvHi, decOvLo)
		s.vLo[1], borrow = bits.Sub64(s.vLo[1], tmpLo, 0)
		s.vHi[1], _ = bits.Sub64(s.vHi[1], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[2], decOvHi, decOvLo)
		s.vLo[2], borrow = bits.Sub64(s.vLo[2], tmpLo, 0)
		s.vHi[2], _ = bits.Sub64(s.vHi[2], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[3], decOvHi, decOvLo)
		s.vLo[3], borrow = bits.Sub64(s.vLo[3], tmpLo, 0)
		s.vHi[3], _ = bits.Sub64(s.vHi[3], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[4], decOvHi, decOvLo)
		s.vLo[4], borrow = bits.Sub64(s.vLo[4], tmpLo, 0)
		s.vHi[4], _ = bits.Sub64(s.vHi[4], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[5], decOvHi, decOvLo)
		s.vLo[5], borrow = bits.Sub64(s.vLo[5], tmpLo, 0)
		s.vHi[5], _ = bits.Sub64(s.vHi[5], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[6], decOvHi, decOvLo)
		s.vLo[6], borrow = bits.Sub64(s.vLo[6], tmpLo, 0)
		s.vHi[6], _ = bits.Sub64(s.vHi[6], tmpHi, borrow)

		tmpHi, tmpLo = MultAndFloor(s.v[7], decOvHi, decOvLo)
		s.vLo[7], borrow = bits.Sub64(s.vLo[7], tmpLo, 0)
		s.vHi[7], _ = bits.Sub64(s.vHi[7], tmpHi, borrow)

		s.vHi[0], s.vLo[0] = RoundTo128Signed(s.vHi[0], s.vLo[0])
		s.vHi[1], s.vLo[1] = RoundTo128Signed(s.vHi[1], s.vLo[1])
		s.vHi[2], s.vLo[2] = RoundTo128Signed(s.vHi[2], s.vLo[2])
		s.vHi[3], s.vLo[3] = RoundTo128Signed(s.vHi[3], s.vLo[3])
		s.vHi[4], s.vLo[4] = RoundTo128Signed(s.vHi[4], s.vLo[4])
		s.vHi[5], s.vLo[5] = RoundTo128Signed(s.vHi[5], s.vLo[5])
		s.vHi[6], s.vLo[6] = RoundTo128Signed(s.vHi[6], s.vLo[6])
		s.vHi[7], s.vLo[7] = RoundTo128Signed(s.vHi[7], s.vLo[7])

		for i := 0; i < lenOut; i++ {
			outi := (*[8]uint64)(unsafe.Pointer(&out[i][k]))
			modOut := s.modOut[i]

			if s.vHi[0]&0x8000000000000000 != 0 {
				outi[0] = mod.Sub(outi[0], mod.Reduce128(-s.vHi[0], -s.vLo[0], modOut), modOut)
			} else {
				outi[0] = mod.Add(outi[0], mod.Reduce128(s.vHi[0], s.vLo[0], modOut), modOut)
			}

			if s.vHi[1]&0x8000000000000000 != 0 {
				outi[1] = mod.Sub(outi[1], mod.Reduce128(-s.vHi[1], -s.vLo[1], modOut), modOut)
			} else {
				outi[1] = mod.Add(outi[1], mod.Reduce128(s.vHi[1], s.vLo[1], modOut), modOut)
			}

			if s.vHi[2]&0x8000000000000000 != 0 {
				outi[2] = mod.Sub(outi[2], mod.Reduce128(-s.vHi[2], -s.vLo[2], modOut), modOut)
			} else {
				outi[2] = mod.Add(outi[2], mod.Reduce128(s.vHi[2], s.vLo[2], modOut), modOut)
			}

			if s.vHi[3]&0x8000000000000000 != 0 {
				outi[3] = mod.Sub(outi[3], mod.Reduce128(-s.vHi[3], -s.vLo[3], modOut), modOut)
			} else {
				outi[3] = mod.Add(outi[3], mod.Reduce128(s.vHi[3], s.vLo[3], modOut), modOut)
			}

			if s.vHi[4]&0x8000000000000000 != 0 {
				outi[4] = mod.Sub(outi[4], mod.Reduce128(-s.vHi[4], -s.vLo[4], modOut), modOut)
			} else {
				outi[4] = mod.Add(outi[4], mod.Reduce128(s.vHi[4], s.vLo[4], modOut), modOut)
			}

			if s.vHi[5]&0x8000000000000000 != 0 {
				outi[5] = mod.Sub(outi[5], mod.Reduce128(-s.vHi[5], -s.vLo[5], modOut), modOut)
			} else {
				outi[5] = mod.Add(outi[5], mod.Reduce128(s.vHi[5], s.vLo[5], modOut), modOut)
			}

			if s.vHi[6]&0x8000000000000000 != 0 {
				outi[6] = mod.Sub(outi[6], mod.Reduce128(-s.vHi[6], -s.vLo[6], modOut), modOut)
			} else {
				outi[6] = mod.Add(outi[6], mod.Reduce128(s.vHi[6], s.vLo[6], modOut), modOut)
			}

			if s.vHi[7]&0x8000000000000000 != 0 {
				outi[7] = mod.Sub(outi[7], mod.Reduce128(-s.vHi[7], -s.vLo[7], modOut), modOut)
			} else {
				outi[7] = mod.Add(outi[7], mod.Reduce128(s.vHi[7], s.vLo[7], modOut), modOut)
			}
		}
	}

	for k := M; k < degree; k++ {
		vHi0, vLo0 := uint64(0), uint64(0)
		for i := 0; i < lenIn; i++ {
			s.buff[i][0] = mod.SMul(in[i][k], s.compInv[i], s.compInvS[i], s.modIn[i])

			tmpHi, tmpLo = MultAndFloor(s.buff[i][0], s.invHi[i], s.invLo[i])
			vLo0, carry = bits.Add64(vLo0, tmpLo, 0)
			vHi0, _ = bits.Add64(vHi0, tmpHi, carry)
		}

		v := RoundTo64(vHi0, vLo0)

		for i := 0; i < lenOut; i++ {
			out[i][k] = mod.SMul(s.modOut[i].Value()-v, s.intOv[i], s.intOvS[i], s.modOut[i])

			for j := 0; j < lenIn; j++ {
				out[i][k] = mod.Add(out[i][k], mod.SMul(s.buff[j][0], s.intSc[i][j], s.intScS[i][j], s.modOut[i]), s.modOut[i])
			}
		}

		vHi0, vLo0 = uint64(0), uint64(0)
		for i := 0; i < lenIn; i++ {
			tmpHi, tmpLo = MultAndFloor(s.buff[i][0], s.decScHi[i], s.decScLo[i])
			vLo0, carry = bits.Add64(vLo0, tmpLo, 0)
			vHi0, _ = bits.Add64(vHi0, tmpHi, carry)
		}

		tmpHi, tmpLo = MultAndFloor(v, s.decOvHi, s.decOvLo)
		vLo0, borrow = bits.Sub64(vLo0, tmpLo, 0)
		vHi0, _ = bits.Sub64(vHi0, tmpHi, borrow)

		vHi0, vLo0 = RoundTo128Signed(vHi0, vLo0)

		for i := 0; i < lenOut; i++ {
			if vHi0&0x8000000000000000 != 0 {
				out[i][k] = mod.Sub(out[i][k], mod.Reduce128(-vHi0, -vLo0, s.modOut[i]), s.modOut[i])
			} else {
				out[i][k] = mod.Add(out[i][k], mod.Reduce128(vHi0, vLo0, s.modOut[i]), s.modOut[i])
			}
		}
	}
}
