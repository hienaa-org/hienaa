package crt_test

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

type BigEmbedder struct {
	modIn      []*big.Int
	modInProd  *big.Int
	modOut     []*big.Int
	modOutProd *big.Int
	crtIn      []*big.Int
}

func NewBigEmbedder(modOut, modIn []*num.Modulus) *BigEmbedder {
	modInBig := make([]*big.Int, len(modIn))
	modInProd := big.NewInt(1)
	for i := range modIn {
		modInBig[i] = new(big.Int).SetUint64(modIn[i].Value())
		modInProd.Mul(modInProd, modInBig[i])
	}

	modOutBig := make([]*big.Int, len(modOut))
	modOutProd := big.NewInt(1)
	for i := range modOut {
		modOutBig[i] = new(big.Int).SetUint64(modOut[i].Value())
		modOutProd.Mul(modOutProd, modOutBig[i])
	}

	crtIn := make([]*big.Int, len(modIn))
	for i := range crtIn {
		qStar := new(big.Int).Div(modInProd, modInBig[i])
		qInv := new(big.Int).ModInverse(qStar, modInBig[i])
		crtIn[i] = new(big.Int).Mul(qStar, qInv)
		crtIn[i].Mod(crtIn[i], modInProd)
	}

	return &BigEmbedder{
		modIn:      modInBig,
		modInProd:  modInProd,
		modOut:     modOutBig,
		modOutProd: modOutProd,
		crtIn:      crtIn,
	}
}

func (e *BigEmbedder) Embed(v [][]uint64) [][]uint64 {
	lenIn := len(v[0])

	vOut := make([][]uint64, len(e.modOut))
	for i := range vOut {
		vOut[i] = make([]uint64, lenIn)
	}

	for i := 0; i < lenIn; i++ {
		xBig := new(big.Int)
		for j := range v {
			xBig.Add(xBig, new(big.Int).Mul(new(big.Int).SetUint64(v[j][i]), e.crtIn[j]))
		}

		xBig.Mod(xBig, e.modInProd)

		if xBig.Cmp(new(big.Int).Rsh(e.modInProd, 1)) == 1 {
			xBig.Sub(xBig, e.modInProd)
		}

		xBig.Mod(xBig, e.modOutProd)

		for j := range vOut {
			vOut[j][i] = new(big.Int).Mod(xBig, e.modOut[j]).Uint64()
		}
	}

	return vOut
}

func (e *BigEmbedder) Scale(v [][]uint64) [][]uint64 {
	scale := big.NewRat(1, 1)
	scale.Num().Set(e.modOutProd)
	scale.Denom().Set(e.modInProd)
	return e.ScaleEmbed(v, scale)
}

func (e *BigEmbedder) ScaleEmbed(v [][]uint64, scale *big.Rat) [][]uint64 {
	lenIn := len(v[0])

	vOut := make([][]uint64, len(e.modOut))
	for i := range vOut {
		vOut[i] = make([]uint64, lenIn)
	}

	for i := 0; i < lenIn; i++ {
		xBig := new(big.Int)
		for j := range v {
			xBig.Add(xBig, new(big.Int).Mul(new(big.Int).SetUint64(v[j][i]), e.crtIn[j]))
		}

		xBig.Mod(xBig, e.modInProd)

		if xBig.Cmp(new(big.Int).Rsh(e.modInProd, 1)) == 1 {
			xBig.Sub(xBig, e.modInProd)
		}

		xBig.Mul(xBig, scale.Num())
		xBigRem := new(big.Int).Mod(xBig, scale.Denom())
		if xBigRem.Cmp(new(big.Int).Rsh(scale.Denom(), 1)) == 1 {
			xBigRem.Sub(xBigRem, scale.Denom())
		}
		xBig.Sub(xBig, xBigRem)
		xBig.Div(xBig, scale.Denom())
		xBig.Mod(xBig, e.modOutProd)

		for j := range vOut {
			vOut[j][i] = new(big.Int).Mod(xBig, e.modOut[j]).Uint64()
		}
	}

	return vOut
}

func genModInOut(modInLen, modOutLen int) (modIn, modOut []*num.Modulus) {
	p := uint64(1)<<20 + 1

	modIn = make([]*num.Modulus, modInLen)
	for i := range modIn {
		p = num.MustNextPrime(p, 2)
		modIn[i] = num.NewModulus(p)
	}

	modOut = make([]*num.Modulus, modOutLen)
	for i := range modOut {
		p = num.MustNextPrime(p, 2)
		modOut[i] = num.NewModulus(p)
	}

	return
}

func TestEmbedder(t *testing.T) {
	modInLen := int(rSrc.SampleN(20))
	modOutLen := int(rSrc.SampleN(20))
	vLen := int(rSrc.SampleN(1 << 5))

	modIn, modOut := genModInOut(modInLen, modOutLen)

	embBig := NewBigEmbedder(modOut, modIn)
	emb := crt.NewEmbedder(modOut, modIn)

	v := randPoly(vLen, modIn).Coeffs

	vOut := emb.EmbedVec(v)
	vOutRef := embBig.Embed(v)

	assert.Equal(t, vOutRef, vOut)
}

func TestScaler(t *testing.T) {
	modInLen := int(rSrc.SampleN(20))
	modOutLen := int(rSrc.SampleN(20))
	vLen := int(rSrc.SampleN(1 << 10))

	modIn, modOut := genModInOut(modInLen, modOutLen)

	embBig := NewBigEmbedder(modOut, modIn)
	sc := crt.NewScaler(modOut, modIn)

	v := randPoly(vLen, modIn).Coeffs

	vOut := sc.ScaleVec(v)
	vOutRef := embBig.Scale(v)

	assert.Equal(t, vOutRef, vOut)
}

func TestScaleEmbedder(t *testing.T) {
	modInLen := int(rSrc.SampleN(20))
	modOutLen := int(rSrc.SampleN(20))
	vLen := int(rSrc.SampleN(1 << 10))

	modIn, modOut := genModInOut(modInLen, modOutLen)

	embBig := NewBigEmbedder(modOut, modIn)

	rSrc := rand.New(rand.NewSource(0))
	scale := big.NewRat(1, 1)
	scale.Num().Rand(rSrc, embBig.modOutProd)
	scale.Denom().Rand(rSrc, embBig.modInProd)

	scEmb := crt.NewScaleEmbedder(modOut, modIn, scale)

	v := randPoly(vLen, modIn).Coeffs

	vOut := scEmb.ScaleEmbedVec(v)
	vOutRef := embBig.ScaleEmbed(v, scale)

	assert.Equal(t, vOutRef, vOut)
}
