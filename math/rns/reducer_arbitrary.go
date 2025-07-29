package rns

import (
	"math"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

type reducerAnyModulus struct {
	// modulus is the modulus.
	modulus *mod.Modulus
	// ambModulus is the modulus for the arbitrary modulus reducer.
	ambModulus []*mod.Modulus
	// embedder is the embedder for the arbitrary modulus reducer.
	embedder *Embedder

	// deg is the degree of the modulo polynomial.
	deg int
	// maxDeg is the maximum degree of the input polynomial.
	maxDeg int
	// diffDeg is the degree of the difference between the cyclotomic polynomial and the intermediate reducing polynomial Q_sp.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the cyclotomic polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT []singleTransformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT []singleTransformer

	// modPoly is the polynomial we target to reduce to.
	modPoly [][]uint64
	// quotientPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is ⌊X^d/modPoly⌋ modulo the modulus, where d is the maximum degree of the input polynomial.
	quotientPoly [][]uint64

	// buf is the polynomial buffer for the reducer.
	buf reducerNTTBuffer
}

func NewReducerAnyModulus(mPoly []uint64, maxDeg int, modulus *mod.Modulus) *reducerAnyModulus {
	deg := len(mPoly) - 1

	var ambModulus []*mod.Modulus
	var embedder *Embedder
	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT []singleTransformer
	var modPoly, quotientPoly [][]uint64
	var buf reducerNTTBuffer

	degNext = num.NextProdPower(uint64(deg), []uint64{2})
	diffDeg = uint64(maxDeg - deg)
	diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

	lenAmbMod := int(math.Ceil((2.0*math.Log2(float64(modulus.Value())) + math.Log2(float64(max(degNext, diffDegNext)))) / 62.0))
	ambModulus = FindPrevNTTPrimes(RingParameters{0, int(2 * max(degNext, diffDegNext)), Cyclic}, 61, lenAmbMod)
	embedder = NewEmbedder(ambModulus, []*mod.Modulus{modulus})

	degNextParams := RingParameters{0, int(degNext), Cyclic}
	degNextNTT = make([]singleTransformer, lenAmbMod)
	for i := range degNextNTT {
		degNextNTT[i] = newCyclicPow235Transformer(degNextParams, ambModulus[i])
	}

	diffDegNextParams := RingParameters{0, int(diffDegNext), Cyclic}
	diffDegNextNTT = make([]singleTransformer, lenAmbMod)
	for i := range diffDegNextNTT {
		diffDegNextNTT[i] = newCyclicPow235Transformer(diffDegNextParams, ambModulus[i])
	}

	modPoly = make([][]uint64, lenAmbMod)
	for i := range modPoly {
		modPoly[i] = make([]uint64, degNext)
		copy(modPoly[i][:deg+1], mPoly)
	}

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	qPoly := quotientPolynomialMod(dividend, mPoly[:deg+1], modulus)
	quotientPoly = make([][]uint64, lenAmbMod)
	for i := range quotientPoly {
		quotientPoly[i] = make([]uint64, diffDegNext)
		copy(quotientPoly[i][:maxDeg-deg+1], qPoly)
	}

	for i := 0; i < lenAmbMod; i++ {
		degNextNTT[i].nttInPlace(modPoly[i])
		diffDegNextNTT[i].nttInPlace(quotientPoly[i])
	}

	buf = newReducerNTTBuffer(lenAmbMod, int(maxDeg+1), int(diffDegNext), int(degNext))

	return &reducerAnyModulus{
		modulus:    modulus,
		ambModulus: ambModulus,
		embedder:   embedder,

		deg:         int(deg),
		maxDeg:      int(maxDeg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly:      modPoly,
		quotientPoly: quotientPoly,

		buf: buf,
	}
}

func (r *reducerAnyModulus) ReduceTo(pOut, pIn []uint64) {
	copy(r.buf.pIn[0], pIn)

	// Compute pQuo = ⌊pIn/X^deg⌋
	for i := 0; i < len(r.buf.pQuo); i++ {
		clear(r.buf.pQuo[i])
		for j := 0; j <= r.diffDeg; j++ {
			r.buf.pQuo[i][j] = r.buf.pIn[0][r.deg+j]
		}
	}

	// Compute pQuo = pQuo × ⌊X^(deg+diffDeg)/\Phi_m(X)⌋
	for i := 0; i < len(r.buf.pQuo); i++ {
		r.diffDegNextNTT[i].nttInPlace(r.buf.pQuo[i])
		mod.MMulLazyVecTo(r.buf.pQuo[i], r.buf.pQuo[i], r.quotientPoly[i], r.ambModulus[i])
		r.diffDegNextNTT[i].invNTTInPlace(r.buf.pQuo[i])
	}
	r.embedder.EmbedVecTo(r.buf.pQuo[0:1], r.buf.pQuo)

	// Compute pRem = ⌊pQuo/X^diffDeg⌋ % (X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.diffDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pQuo[0][r.diffDeg+j] = mod.Add(r.buf.pQuo[0][r.diffDeg+j], r.buf.pQuo[0][r.diffDeg+i*r.degNext+j], r.modulus)
			r.buf.pQuo[0][r.diffDeg+i*r.degNext+j] = 0
		}
	}

	for i := 0; i < len(r.buf.pRem); i++ {
		clear(r.buf.pRem[i])
		for j := 0; j < min(r.diffDeg, r.degNext); j++ {
			r.buf.pRem[i][j] = r.buf.pQuo[0][r.diffDeg+j]
		}
	}

	// Compute pRem = pRem × quotientPoly (mod X^degNext - 1)
	for i := 0; i < len(r.buf.pRem); i++ {
		r.degNextNTT[i].nttInPlace(r.buf.pRem[i])
		mod.MMulLazyVecTo(r.buf.pRem[i], r.buf.pRem[i], r.modPoly[i], r.ambModulus[i])
		r.degNextNTT[i].invNTTInPlace(r.buf.pRem[i])
	}
	r.embedder.EmbedVecTo(r.buf.pRem[0:1], r.buf.pRem)

	// Compute pIn = pIn (mod X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.maxDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.maxDeg {
				break
			}
			r.buf.pIn[0][j] = mod.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.modulus)
			r.buf.pIn[0][i*r.degNext+j] = 0
		}
	}

	// Compute pOut = pIn - pRem
	for i := 0; i < r.deg; i++ {
		pOut[i] = mod.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.modulus)
	}
}

type reducerNTTModulus struct {
	// modulus is the modulus.
	modulus *mod.Modulus

	// deg is the degree of the modulo polynomial.
	deg int
	// maxDeg is the maximum degree of the input polynomial.
	maxDeg int
	// diffDeg is the degree of the difference between the modulo polynomial and the maximum degree of the input polynomial.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the modulo polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT singleTransformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT singleTransformer

	// modPoly is the polynomial we target to reduce to.
	modPoly []uint64
	// quotientPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is ⌊X^d/modPoly⌋ modulo the modulus, where d is the maximum degree of the input polynomial.
	quotientPoly []uint64

	// buf is the polynomial buffer for the reducer.
	buf reducerNTTBuffer
}

func NewReducerNTTModulus(modPoly []uint64, maxDeg int, modulus *mod.Modulus) *reducerNTTModulus {
	deg := len(modPoly) - 1

	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT singleTransformer
	var quotientPoly []uint64
	var buf reducerNTTBuffer

	degNext = num.NextProdPower(uint64(deg), []uint64{2})
	diffDeg = uint64(maxDeg - deg)
	diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

	degNextParams := RingParameters{0, int(degNext), Cyclic}
	degNextNTT = newCyclicPow235Transformer(degNextParams, modulus)

	diffDegNextParams := RingParameters{0, int(diffDegNext), Cyclic}
	diffDegNextNTT = newCyclicPow235Transformer(diffDegNextParams, modulus)

	modPoly = append(modPoly, make([]uint64, int(degNext)-deg-1)...)

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	quotientPoly = quotientPolynomialMod(dividend, modPoly[:deg+1], modulus)
	quotientPoly = append(quotientPoly, make([]uint64, int(diffDegNext)-maxDeg+deg-1)...)

	degNextNTT.nttInPlace(modPoly)
	diffDegNextNTT.nttInPlace(quotientPoly)

	buf = newReducerNTTBuffer(1, int(maxDeg+1), int(diffDegNext), int(degNext))

	return &reducerNTTModulus{
		modulus: modulus,

		deg:         int(deg),
		maxDeg:      int(maxDeg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly:      modPoly,
		quotientPoly: quotientPoly,

		buf: buf,
	}
}

func (r *reducerNTTModulus) ReduceTo(pOut, pIn []uint64) {
	copy(r.buf.pIn[0], pIn)

	// Compute pQuo = ⌊pIn/X^deg⌋
	clear(r.buf.pQuo[0])
	for i := 0; i <= r.diffDeg; i++ {
		r.buf.pQuo[0][i] = r.buf.pIn[0][r.deg+i]
	}

	// Compute pQuo = pQuo × ⌊X^(deg+diffDeg)/modPoly⌋
	r.diffDegNextNTT.nttInPlace(r.buf.pQuo[0])
	mod.MMulLazyVecTo(r.buf.pQuo[0], r.buf.pQuo[0], r.quotientPoly, r.modulus)
	r.diffDegNextNTT.invNTTInPlace(r.buf.pQuo[0])

	// Compute pRem = ⌊pQuo/X^diffDeg⌋ % (X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.diffDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pQuo[0][r.diffDeg+j] = mod.Add(r.buf.pQuo[0][r.diffDeg+j], r.buf.pQuo[0][r.diffDeg+i*r.degNext+j], r.modulus)
			r.buf.pQuo[0][r.diffDeg+i*r.degNext+j] = 0
		}
	}

	clear(r.buf.pRem[0])
	for i := 0; i < min(r.diffDeg, r.degNext); i++ {
		r.buf.pRem[0][i] = r.buf.pQuo[0][r.diffDeg+i]
	}

	// Compute pRem = pRem × modPoly (mod X^degNext - 1)
	r.degNextNTT.nttInPlace(r.buf.pRem[0])
	mod.MMulLazyVecTo(r.buf.pRem[0], r.buf.pRem[0], r.modPoly, r.modulus)
	r.degNextNTT.invNTTInPlace(r.buf.pRem[0])

	// Compute pIn = pIn (mod X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.maxDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.maxDeg {
				break
			}
			r.buf.pIn[0][j] = mod.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.modulus)
			r.buf.pIn[0][i*r.degNext+j] = 0
		}
	}

	// Compute pOut = pIn - pRem
	for i := 0; i < r.deg; i++ {
		pOut[i] = mod.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.modulus)
	}
}
