package crt

import (
	"math"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type reducerAnyModulus struct {
	// mod is the mod.
	mod *num.Modulus
	// ambMod is the modulus for the arbitrary modulus reducer.
	ambMod []*num.Modulus
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
	diffDegNextNTT []dft.Transformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT []dft.Transformer

	// modPoly is the polynomial we target to reduce to.
	modPoly [][]uint64
	// quoPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is floor(X^d/modPoly) modulo the modulus, where d is the maximum degree of the input polynomial.
	quoPoly [][]uint64

	buf reducerBuffer
}

// newReducerAnyModulus creates a new [reducerAnyModulus].
func newReducerAnyModulus(maxDeg int, mod *num.Modulus, modPoly []int64) *reducerAnyModulus {
	modPolyReduced := make([]uint64, len(modPoly))
	for i := range modPolyReduced {
		modPolyReduced[i] = dftops.ReduceInt(modPoly[i]%int64(mod.Value()), mod)
	}

	deg := len(modPoly) - 1

	degNext := num.NextProdPower(deg, []int{2})
	diffDeg := maxDeg - deg
	diffDegNext := num.NextProdPower(2*diffDeg+1, []int{2})

	maxBits := 2*num.Log2(mod.Value()) + num.Log2(max(degNext, diffDegNext))
	lenAmbMod := int(math.Ceil(maxBits / num.MaxModulusBits))
	ambMod := dft.FindPrevNTTPrimes(dft.NewCyclicParameters(2*max(degNext, diffDegNext)), num.MaxModulusBits, lenAmbMod)
	embedder := NewEmbedder([]*num.Modulus{mod}, ambMod)

	degNextParams := dft.NewCyclicParameters(degNext)
	degNextNTT := make([]dft.Transformer, lenAmbMod)
	for i := range degNextNTT {
		degNextNTT[i] = dft.NewTransformer(degNextParams, ambMod[i])
	}

	diffDegNextParams := dft.NewCyclicParameters(diffDegNext)
	diffDegNextNTT := make([]dft.Transformer, lenAmbMod)
	for i := range diffDegNextNTT {
		diffDegNextNTT[i] = dft.NewTransformer(diffDegNextParams, ambMod[i])
	}

	ambModPoly := make([][]uint64, lenAmbMod)
	for i := range ambModPoly {
		ambModPoly[i] = make([]uint64, degNext)
		copy(ambModPoly[i][:deg+1], modPolyReduced)
	}

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	quoPoly := dftops.Quotient(dividend, modPolyReduced[:deg+1], mod)
	ambQuoPoly := make([][]uint64, lenAmbMod)
	for i := range ambQuoPoly {
		ambQuoPoly[i] = make([]uint64, diffDegNext)
		copy(ambQuoPoly[i][:maxDeg-deg+1], quoPoly)
	}

	for i := 0; i < lenAmbMod; i++ {
		degNextNTT[i].ForwardInPlace(ambModPoly[i])
		diffDegNextNTT[i].ForwardInPlace(ambQuoPoly[i])
	}

	return &reducerAnyModulus{
		mod:      mod,
		ambMod:   ambMod,
		embedder: embedder,

		deg:         deg,
		maxDeg:      maxDeg,
		diffDeg:     diffDeg,
		diffDegNext: diffDegNext,
		degNext:     degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly: ambModPoly,
		quoPoly: ambQuoPoly,

		buf: newReducerBuffer(lenAmbMod, maxDeg+1, diffDegNext, degNext),
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
		r.diffDegNextNTT[i].ForwardInPlace(r.buf.pQuo[i])
		vec.MMulTo(r.buf.pQuo[i], r.buf.pQuo[i], r.quoPoly[i], r.ambMod[i])
		r.diffDegNextNTT[i].InverseInPlace(r.buf.pQuo[i])
	}
	r.embedder.EmbedVecTo(r.buf.pQuo[0:1], r.buf.pQuo)

	// Compute pRem = floor(pQuo/X^diffDeg) % (X^degNext - 1)
	for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pQuo[0][r.diffDeg+j] = num.Add(r.buf.pQuo[0][r.diffDeg+j], r.buf.pQuo[0][r.diffDeg+i*r.degNext+j], r.mod)
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
		r.degNextNTT[i].ForwardInPlace(r.buf.pRem[i])
		vec.MMulTo(r.buf.pRem[i], r.buf.pRem[i], r.modPoly[i], r.ambMod[i])
		r.degNextNTT[i].InverseInPlace(r.buf.pRem[i])
	}
	r.embedder.EmbedVecTo(r.buf.pRem[0:1], r.buf.pRem)

	// Compute pIn = pIn (mod X^degNext - 1)
	for i := 1; i <= num.DivCeil(r.maxDeg, r.degNext); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.maxDeg {
				break
			}
			r.buf.pIn[0][j] = num.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.mod)
			r.buf.pIn[0][i*r.degNext+j] = 0
		}
	}

	// Compute pOut = pIn - pRem
	for i := 0; i < r.deg; i++ {
		pOut[i] = num.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.mod)
	}
}

func (r *reducerAnyModulus) safeCopy() reducer {
	diffDegNextNTT := make([]dft.Transformer, len(r.diffDegNextNTT))
	for i := range diffDegNextNTT {
		diffDegNextNTT[i] = r.diffDegNextNTT[i].SafeCopy()
	}

	degNextNTT := make([]dft.Transformer, len(r.degNextNTT))
	for i := range degNextNTT {
		degNextNTT[i] = r.degNextNTT[i].SafeCopy()
	}

	return &reducerAnyModulus{
		mod:      r.mod,
		ambMod:   r.ambMod,
		embedder: r.embedder.SafeCopy(),

		deg:         r.deg,
		maxDeg:      r.maxDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly: r.modPoly,
		quoPoly: r.quoPoly,

		buf: newReducerBuffer(len(r.ambMod), r.maxDeg+1, r.diffDegNext, r.degNext),
	}
}

type reducerNTTModulus struct {
	mod *num.Modulus

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
	diffDegNextNTT dft.Transformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT dft.Transformer

	// modPoly is the polynomial we target to reduce to.
	modPoly []uint64
	// quoPoly is rounding of a monomial over the modPoly modulo the modulus.
	// Precisely, it is floor(X^d/modPoly) modulo the modulus, where d is the maximum degree of the input polynomial.
	quoPoly []uint64

	buf reducerBuffer
}

// newReducerNTTModulus creates a new [reducerNTTModulus].
func newReducerNTTModulus(maxDeg int, mod *num.Modulus, modPoly []int64) *reducerNTTModulus {
	modPolyReduced := make([]uint64, len(modPoly))
	for i := range modPolyReduced {
		modPolyReduced[i] = dftops.ReduceInt(modPoly[i]%int64(mod.Value()), mod)
	}

	deg := len(modPoly) - 1

	degNext := num.NextProdPower(deg, []int{2})
	diffDeg := maxDeg - deg
	diffDegNext := num.NextProdPower(2*diffDeg+1, []int{2})

	degNextParams := dft.NewCyclicParameters(int(degNext))
	degNextNTT := dft.NewTransformer(degNextParams, mod)

	diffDegNextParams := dft.NewCyclicParameters(int(diffDegNext))
	diffDegNextNTT := dft.NewTransformer(diffDegNextParams, mod)

	modPolyExtended := make([]uint64, degNext)
	copy(modPolyExtended, modPolyReduced)

	dividend := make([]uint64, maxDeg+1)
	dividend[maxDeg] = 1
	quoPoly := dftops.Quotient(dividend, modPolyExtended[:deg+1], mod)
	quoPoly = append(quoPoly, make([]uint64, int(diffDegNext)-maxDeg+deg-1)...)

	degNextNTT.ForwardInPlace(modPolyExtended)
	diffDegNextNTT.ForwardInPlace(quoPoly)

	return &reducerNTTModulus{
		mod: mod,

		deg:         deg,
		maxDeg:      maxDeg,
		diffDeg:     diffDeg,
		diffDegNext: diffDegNext,
		degNext:     degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		modPoly: modPolyExtended,
		quoPoly: quoPoly,

		buf: newReducerBuffer(1, maxDeg+1, diffDegNext, degNext),
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
	r.diffDegNextNTT.ForwardInPlace(r.buf.pQuo[0])
	vec.MMulTo(r.buf.pQuo[0], r.buf.pQuo[0], r.quoPoly, r.mod)
	r.diffDegNextNTT.InverseInPlace(r.buf.pQuo[0])

	// Compute pRem = floor(pQuo/X^diffDeg) % (X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.diffDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pQuo[0][r.diffDeg+j] = num.Add(r.buf.pQuo[0][r.diffDeg+j], r.buf.pQuo[0][r.diffDeg+i*r.degNext+j], r.mod)
			r.buf.pQuo[0][r.diffDeg+i*r.degNext+j] = 0
		}
	}

	clear(r.buf.pRem[0])
	for i := 0; i < min(r.diffDeg, r.degNext); i++ {
		r.buf.pRem[0][i] = r.buf.pQuo[0][r.diffDeg+i]
	}

	// Compute pRem = pRem * modPoly (mod X^degNext - 1)
	r.degNextNTT.ForwardInPlace(r.buf.pRem[0])
	vec.MMulTo(r.buf.pRem[0], r.buf.pRem[0], r.modPoly, r.mod)
	r.degNextNTT.InverseInPlace(r.buf.pRem[0])

	// Compute pIn = pIn (mod X^degNext - 1)
	for i := 1; i <= int(math.Ceil(float64(r.maxDeg)/float64(r.degNext))); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.maxDeg {
				break
			}
			r.buf.pIn[0][j] = num.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.mod)
			r.buf.pIn[0][i*r.degNext+j] = 0
		}
	}

	// Compute pOut = pIn - pRem
	for i := 0; i < r.deg; i++ {
		pOut[i] = num.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.mod)
	}
}

func (r *reducerNTTModulus) safeCopy() reducer {
	return &reducerNTTModulus{
		mod: r.mod,

		deg:         r.deg,
		maxDeg:      r.maxDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: r.diffDegNextNTT.SafeCopy(),
		degNextNTT:     r.degNextNTT.SafeCopy(),

		modPoly: r.modPoly,
		quoPoly: r.quoPoly,

		buf: newReducerBuffer(1, r.maxDeg+1, r.diffDegNext, r.degNext),
	}
}
