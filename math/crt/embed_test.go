package crt_test

import (
	"math/big"
	"math/rand/v2"
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
	primes := make([]uint64, modInLen+modOutLen)
	start := uint64(num.MaxModulus) | 1
	for i := range primes {
		primes[i] = num.PrevPrime(start-2, 2)
		start = primes[i]
	}

	rSrc := rand.New(rand.NewChaCha8([32]byte{}))
	rSrc.Shuffle(len(primes), func(i, j int) {
		primes[i], primes[j] = primes[j], primes[i]
	})

	modIn = make([]*num.Modulus, modInLen)
	for i, q := range primes[:modInLen] {
		modIn[i] = num.NewModulus(q)
	}
	modOut = make([]*num.Modulus, modOutLen)
	for i, q := range primes[modInLen:] {
		modOut[i] = num.NewModulus(q)
	}

	return
}

func TestEmbedder(t *testing.T) {
	modInLen := rSrc.SampleN(1, 20)
	modOutLen := rSrc.SampleN(1, 20)
	vLen := rSrc.SampleN(1, 1<<10)

	modIn, modOut := genModInOut(modInLen, modOutLen)

	embBig := NewBigEmbedder(modOut, modIn)
	emb := crt.NewVecEmbedder(modOut, modIn)

	v := randPoly(vLen, modIn).Coeffs

	vOut := emb.Embed(v)
	vOutRef := embBig.Embed(v)

	assert.Equal(t, vOutRef, vOut)
}

func TestScaler(t *testing.T) {
	modInLen := rSrc.SampleN(3, 20)
	modOutLen := rSrc.SampleN(3, 20)
	vLen := rSrc.SampleN(1, 1<<10)

	modIn, modOut := genModInOut(modInLen, modOutLen)

	embBig := NewBigEmbedder(modOut, modIn)
	sc := crt.NewVecScaler(modOut, modIn)

	v := randPoly(vLen, modIn).Coeffs

	vOut := sc.Scale(v)
	vOutRef := embBig.Scale(v)

	assert.Equal(t, vOutRef, vOut)

	t.Run("modIn|modOut", func(t *testing.T) {
		modIn = modOut[:len(modOut)>>1]

		embBig := NewBigEmbedder(modOut, modIn)
		sc := crt.NewVecScaler(modOut, modIn)

		v := randPoly(vLen, modIn).Coeffs

		vOut := sc.Scale(v)
		vOutRef := embBig.Scale(v)

		assert.Equal(t, vOutRef, vOut)
	})

	t.Run("modOut|modIn", func(t *testing.T) {
		modOut = modIn[:len(modIn)>>1]

		embBig := NewBigEmbedder(modOut, modIn)
		sc := crt.NewVecScaler(modOut, modIn)

		v := randPoly(vLen, modIn).Coeffs

		vOut := sc.Scale(v)
		vOutRef := embBig.Scale(v)

		assert.Equal(t, vOutRef, vOut)
	})
}
