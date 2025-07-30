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
	// quoPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is floor(X^d/modPoly) modulo the modulus, where d is the maximum degree of the input polynomial.
	quoPoly [][]uint64

	buf reducerBuffer
}

// newReducerAnyModulus creates a new [reducerAnyModulus].
func newReducerAnyModulus(maxDeg int, modPoly []uint64, modulus *mod.Modulus) *reducerAnyModulus {
	deg := len(modPoly) - 1

	degNext := num.NextProdPower(uint64(deg), []uint64{2})
	diffDeg := uint64(maxDeg - deg)
	diffDegNext := num.NextProdPower(2*diffDeg+1, []uint64{2})

	lenAmbMod := int(math.Ceil((2.0*math.Log2(float64(modulus.Value())) + math.Log2(float64(max(degNext, diffDegNext)))) / mod.MaxModulusBits))
	ambModulus := FindPrevNTTPrimes(NewCyclicParameters(int(2*max(degNext, diffDegNext))), mod.MaxModulusBits, lenAmbMod)
	embedder := NewEmbedder(ambModulus, []*mod.Modulus{modulus})

	degNextParams := NewCyclicParameters(int(degNext))
	degNextNTT := make([]singleTransformer, lenAmbMod)
	for i := range degNextNTT {
		degNextNTT[i] = newCyclicPow235Transformer(degNextParams, ambModulus[i])
	}

	diffDegNextParams := NewCyclicParameters(int(diffDegNext))
	diffDegNextNTT := make([]singleTransformer, lenAmbMod)
	for i := range diffDegNextNTT {
		diffDegNextNTT[i] = newCyclicPow235Transformer(diffDegNextParams, ambModulus[i])
	}

	ambModPoly := make([][]uint64, lenAmbMod)
	for i := range ambModPoly {
		ambModPoly[i] = make([]uint64, degNext)
		copy(ambModPoly[i][:deg+1], modPoly)
	}

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	quoPoly := quotientPolynomial(dividend, modPoly[:deg+1], modulus)
	ambQuoPoly := make([][]uint64, lenAmbMod)
	for i := range ambQuoPoly {
		ambQuoPoly[i] = make([]uint64, diffDegNext)
		copy(ambQuoPoly[i][:maxDeg-deg+1], quoPoly)
	}

	for i := 0; i < lenAmbMod; i++ {
		degNextNTT[i].nttInPlace(ambModPoly[i])
		diffDegNextNTT[i].nttInPlace(ambQuoPoly[i])
	}

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

		modPoly: ambModPoly,
		quoPoly: ambQuoPoly,

		buf: newReducerBuffer(lenAmbMod, int(maxDeg+1), int(diffDegNext), int(degNext)),
	}
}

func (r *reducerAnyModulus) reduceTo(pOut, p []uint64) {
	copy(r.buf.pIn[0], p)

	// Compute pQuo = floor(pIn/X^deg)
	for i := 0; i < len(r.buf.pQuo); i++ {
		clear(r.buf.pQuo[i])
		for j := 0; j <= r.diffDeg; j++ {
			r.buf.pQuo[i][j] = r.buf.pIn[0][r.deg+j]
		}
	}

	// Compute pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
	for i := 0; i < len(r.buf.pQuo); i++ {
		r.diffDegNextNTT[i].nttInPlace(r.buf.pQuo[i])
		mod.MMulLazyVecTo(r.buf.pQuo[i], r.buf.pQuo[i], r.quoPoly[i], r.ambModulus[i])
		r.diffDegNextNTT[i].invNTTInPlace(r.buf.pQuo[i])
	}
	r.embedder.EmbedVecTo(r.buf.pQuo[0:1], r.buf.pQuo)

	// Compute pRem = floor(pQuo/X^diffDeg) % (X^degNext - 1)
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

	// Compute pRem = pRem * quoPoly (mod X^degNext - 1)
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

func (r *reducerAnyModulus) safeCopy() reducer {
	diffDegNextNTT := make([]singleTransformer, len(r.diffDegNextNTT))
	for i := range diffDegNextNTT {
		diffDegNextNTT[i] = r.diffDegNextNTT[i].safeCopy()
	}

	degNextNTT := make([]singleTransformer, len(r.degNextNTT))
	for i := range degNextNTT {
		degNextNTT[i] = r.degNextNTT[i].safeCopy()
	}

	return &reducerAnyModulus{
		modulus:    r.modulus,
		ambModulus: r.ambModulus,
		embedder:   r.embedder.SafeCopy(),

		deg:         r.deg,
		maxDeg:      r.maxDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly: r.modPoly,
		quoPoly: r.quoPoly,

		buf: newReducerBuffer(len(r.ambModulus), r.maxDeg+1, r.diffDegNext, r.degNext),
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
	// quoPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is floor(X^d/modPoly) modulo the modulus, where d is the maximum degree of the input polynomial.
	quoPoly []uint64

	buf reducerBuffer
}

// newReducerNTTModulus creates a new [reducerNTTModulus].
func newReducerNTTModulus(maxDeg int, modPoly []uint64, modulus *mod.Modulus) *reducerNTTModulus {
	deg := len(modPoly) - 1

	degNext := num.NextProdPower(uint64(deg), []uint64{2})
	diffDeg := uint64(maxDeg - deg)
	diffDegNext := num.NextProdPower(2*diffDeg+1, []uint64{2})

	degNextParams := NewCyclicParameters(int(degNext))
	degNextNTT := newCyclicPow235Transformer(degNextParams, modulus)

	diffDegNextParams := NewCyclicParameters(int(diffDegNext))
	diffDegNextNTT := newCyclicPow235Transformer(diffDegNextParams, modulus)

	modPolyExtended := make([]uint64, degNext)
	copy(modPolyExtended, modPoly)

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	quoPoly := quotientPolynomial(dividend, modPolyExtended[:deg+1], modulus)
	quoPoly = append(quoPoly, make([]uint64, int(diffDegNext)-maxDeg+deg-1)...)

	degNextNTT.nttInPlace(modPolyExtended)
	diffDegNextNTT.nttInPlace(quoPoly)

	return &reducerNTTModulus{
		modulus: modulus,

		deg:         int(deg),
		maxDeg:      int(maxDeg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly: modPolyExtended,
		quoPoly: quoPoly,

		buf: newReducerBuffer(1, int(maxDeg+1), int(diffDegNext), int(degNext)),
	}
}

func (r *reducerNTTModulus) reduceTo(pOut, p []uint64) {
	copy(r.buf.pIn[0], p)

	// Compute pQuo = floor(pIn/X^deg)
	clear(r.buf.pQuo[0])
	for i := 0; i <= r.diffDeg; i++ {
		r.buf.pQuo[0][i] = r.buf.pIn[0][r.deg+i]
	}

	// Compute pQuo = pQuo * floor(X^(deg+diffDeg)/modPoly)
	r.diffDegNextNTT.nttInPlace(r.buf.pQuo[0])
	mod.MMulLazyVecTo(r.buf.pQuo[0], r.buf.pQuo[0], r.quoPoly, r.modulus)
	r.diffDegNextNTT.invNTTInPlace(r.buf.pQuo[0])

	// Compute pRem = floor(pQuo/X^diffDeg) % (X^degNext - 1)
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

	// Compute pRem = pRem * modPoly (mod X^degNext - 1)
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

func (r *reducerNTTModulus) safeCopy() reducer {
	return &reducerNTTModulus{
		modulus: r.modulus,

		deg:         r.deg,
		maxDeg:      r.maxDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: r.diffDegNextNTT.safeCopy(),
		degNextNTT:     r.degNextNTT.safeCopy(),

		modPoly: r.modPoly,
		quoPoly: r.quoPoly,

		buf: newReducerBuffer(1, r.maxDeg+1, r.diffDegNext, r.degNext),
	}
}
