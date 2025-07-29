package rns

import (
	"math"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
)

type cyclotomicReducerAnyModulus struct {
	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isPrimePow is true if the cyclotomic polynomial is a prime power.
	isPrimePow bool

	// params is the ring parameters.
	params RingParameters
	// modulus is the modulus for the cyclotomic reducer.
	modulus *mod.Modulus
	// ambModulus is the modulus for the arbitrary modulus reducer.
	ambModulus []*mod.Modulus
	// embedder is the embedder for the arbitrary modulus reducer.
	embedder *Embedder

	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
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

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly [][]uint64
	// quotientPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is ⌊X^d_qs/\Phi_m(X)⌋ modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quotientPoly [][]uint64

	// buf is the polynomial buffer for the reducer.
	buf reducerNTTBuffer
}

func NewCyclotomicReducerAnyModulus(params RingParameters, modulus *mod.Modulus) *cyclotomicReducerAnyModulus {
	cycloDeg, deg := uint64(params.CycloDegree()), uint64(params.Degree())

	factors := num.Factor(cycloDeg)
	leastFactor := cycloDeg
	for key := range factors {
		if key < leastFactor {
			leastFactor = key
		}
	}
	redDeg := cycloDeg - cycloDeg/leastFactor

	var isPrimePower bool
	var ambModulus []*mod.Modulus
	var embedder *Embedder
	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT []singleTransformer
	var cycloPoly, quotientPoly [][]uint64
	var buf reducerNTTBuffer

	if redDeg == deg {
		isPrimePower = true
		diffDeg, diffDegNext, degNext = 0, 0, 0
		diffDegNextNTT = []singleTransformer{}
		ambModulus = make([]*mod.Modulus, 0)
		embedder = nil
		degNextNTT = []singleTransformer{}
		cycloPoly, quotientPoly = [][]uint64{}, [][]uint64{}
		buf = newReducerNTTBuffer(0, 0, 0, 0)
	} else {
		degNext = num.NextProdPower(uint64(deg), []uint64{2})
		diffDeg = redDeg - deg
		diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

		lenAmbMod := int(math.Ceil((2.0*math.Log2(float64(modulus.Value())) + math.Log2(float64(max(degNext, diffDegNext)))) / 62.0))
		ambModulus = FindPrevNTTPrimes(params, 61, lenAmbMod)
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

		cycloPoly = make([][]uint64, lenAmbMod)
		cycloPolyInt := computeCyclotomicPolynomial(cycloDeg)
		for i := range cycloPoly {
			cycloPoly[i] = make([]uint64, degNext)
			for j := range cycloPolyInt {
				cycloPoly[i][j] = reduceInt(cycloPolyInt[j], ambModulus[i])
			}
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quotientPoly = make([][]uint64, lenAmbMod)
		for i := range quotientPoly {
			quotientPoly[i] = append(quotientPolynomialMod(dividend, cycloPoly[i][:deg+1], ambModulus[i]), make([]uint64, int(diffDegNext)-len(quotientPoly[i]))...)
		}

		for i := 0; i < lenAmbMod; i++ {
			degNextNTT[i].nttInPlace(cycloPoly[i])
			diffDegNextNTT[i].nttInPlace(quotientPoly[i])
		}

		buf = newReducerNTTBuffer(lenAmbMod, int(cycloDeg), int(diffDegNext), int(degNext))
	}

	return &cyclotomicReducerAnyModulus{
		leastFac:   int(leastFactor),
		isPrimePow: isPrimePower,

		params:     params,
		modulus:    modulus,
		ambModulus: ambModulus,
		embedder:   embedder,

		redDeg:      int(redDeg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly:    cycloPoly,
		quotientPoly: quotientPoly,

		buf: buf,
	}
}

func (r *cyclotomicReducerAnyModulus) ReduceTo(pOut, pIn []uint64) {
	copy(r.buf.pIn[0], pIn)

	leastFac := r.leastFac
	cycloDeg := r.params.CycloDegree()
	deg := r.params.Degree()
	skip := cycloDeg / leastFac

	for j := 0; j < skip; j++ {
		for i := 0; i < leastFac-1; i++ {
			r.buf.pIn[0][i*skip+j] = mod.Sub(r.buf.pIn[0][i*skip+j], r.buf.pIn[0][cycloDeg-skip+j], r.modulus)
		}
		r.buf.pIn[0][cycloDeg-skip+j] = 0
	}

	if !r.isPrimePow {
		// Compute pQuo = ⌊pIn/X^deg⌋
		for i := 0; i < len(r.buf.pQuo); i++ {
			clear(r.buf.pQuo[i])
			for j := 0; j < r.diffDeg; j++ {
				r.buf.pQuo[i][j] = r.buf.pIn[0][deg+j]
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
			mod.MMulLazyVecTo(r.buf.pRem[i], r.buf.pRem[i], r.cycloPoly[i], r.ambModulus[i])
			r.degNextNTT[i].invNTTInPlace(r.buf.pRem[i])
		}
		r.embedder.EmbedVecTo(r.buf.pRem[0:1], r.buf.pRem)

		// Compute pIn = pIn (mod X^degNext - 1)
		for i := 1; i <= int(math.Ceil(float64(r.redDeg)/float64(r.degNext))); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[0][j] = mod.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.modulus)
				r.buf.pIn[0][i*r.degNext+j] = 0
			}
		}

		// Compute pOut = pIn - pRem
		for i := 0; i < r.params.Degree(); i++ {
			pOut[i] = mod.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.modulus)
		}
	} else {
		copy(pOut, r.buf.pIn[0][:r.params.Degree()])
	}
}

// cyclotomicReducerNTTModulus reduces a polynomial modulo a cyclotomic polynomial.
// In other words, it computes a(X) mod \Phi_m(X), where a(X) is at most degree m-1.
// It uses Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
// In this implementation, Q_sp = (Xᵐ-1)/(X^(m/p)-1) for the smallest prime factor p.
type cyclotomicReducerNTTModulus struct {
	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isPrimePow is true if the cyclotomic polynomial is a prime power.
	isPrimePow bool

	// params is the ring parameters.
	params RingParameters
	// modulus is the modulus.
	modulus *mod.Modulus

	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
	// diffDeg is the degree of the difference between the cyclotomic polynomial and the intermediate reducing polynomial Q_sp.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the cyclotomic polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT singleTransformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT singleTransformer

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly []uint64
	// quotientPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is ⌊X^d_qs/\Phi_m(X)⌋ modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quotientPoly []uint64

	// buf is the polynomial buffer for the reducer.
	buf reducerNTTBuffer
}

func NewCyclotomicReducerNTTModulus(params RingParameters, modulus *mod.Modulus) *cyclotomicReducerNTTModulus {
	cycloDeg, deg := uint64(params.CycloDegree()), uint64(params.Degree())

	factors := num.Factor(cycloDeg)
	leastFactor := cycloDeg
	for key := range factors {
		if key < leastFactor {
			leastFactor = key
		}
	}
	redDeg := cycloDeg - cycloDeg/leastFactor

	var isPrimePower bool
	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT singleTransformer
	var cycloPoly, quotientPoly []uint64
	var buf reducerNTTBuffer

	if redDeg == deg {
		isPrimePower = true
		diffDeg, diffDegNext, degNext = 0, 0, 0
		diffDegNextNTT, degNextNTT = newTrivialTransformer(modulus), newTrivialTransformer(modulus)
		cycloPoly, quotientPoly = []uint64{}, []uint64{}
		buf = newReducerNTTBuffer(0, 0, 0, 0)
	} else {
		degNext = num.NextProdPower(uint64(deg), []uint64{2})
		diffDeg = redDeg - deg
		diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

		degNextParams := RingParameters{0, int(degNext), Cyclic}
		degNextNTT = newCyclicPow235Transformer(degNextParams, modulus)

		diffDegNextParams := RingParameters{0, int(diffDegNext), Cyclic}
		diffDegNextNTT = newCyclicPow235Transformer(diffDegNextParams, modulus)

		cycloPoly = make([]uint64, degNext)
		cycloPolyInt := computeCyclotomicPolynomial(cycloDeg)
		for i := range cycloPolyInt {
			cycloPoly[i] = reduceInt(cycloPolyInt[i], modulus)
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quotientPoly = append(quotientPolynomialMod(dividend, cycloPoly[:deg+1], modulus), make([]uint64, int(diffDegNext)-len(quotientPoly))...)

		degNextNTT.nttInPlace(cycloPoly)
		diffDegNextNTT.nttInPlace(quotientPoly)

		buf = newReducerNTTBuffer(1, int(cycloDeg), int(diffDegNext), int(degNext))
	}

	return &cyclotomicReducerNTTModulus{
		leastFac:   int(leastFactor),
		isPrimePow: isPrimePower,

		params:  params,
		modulus: modulus,

		redDeg:      int(redDeg),
		diffDeg:     int(diffDeg),
		diffDegNext: int(diffDegNext),
		degNext:     int(degNext),

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly:    cycloPoly,
		quotientPoly: quotientPoly,

		buf: buf,
	}
}

func (r *cyclotomicReducerNTTModulus) ReduceTo(pOut, pIn []uint64) {
	copy(r.buf.pIn[0], pIn)

	leastFac := r.leastFac
	cycloDeg := r.params.CycloDegree()
	deg := r.params.Degree()
	skip := cycloDeg / leastFac

	for j := 0; j < skip; j++ {
		for i := 0; i < leastFac-1; i++ {
			r.buf.pIn[0][i*skip+j] = mod.Sub(r.buf.pIn[0][i*skip+j], r.buf.pIn[0][cycloDeg-skip+j], r.modulus)
		}
		r.buf.pIn[0][cycloDeg-skip+j] = 0
	}

	if !r.isPrimePow {
		// Compute pQuo = ⌊pIn/X^deg⌋
		clear(r.buf.pQuo[0])
		for i := 0; i < r.diffDeg; i++ {
			r.buf.pQuo[0][i] = r.buf.pIn[0][deg+i]
		}

		// Compute pQuo = pQuo × ⌊X^(deg+diffDeg)/\Phi_m(X)⌋
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

		// Compute pRem = pRem × quotientPoly (mod X^degNext - 1)
		r.degNextNTT.nttInPlace(r.buf.pRem[0])
		mod.MMulLazyVecTo(r.buf.pRem[0], r.buf.pRem[0], r.cycloPoly, r.modulus)
		r.degNextNTT.invNTTInPlace(r.buf.pRem[0])

		// Compute pIn = pIn (mod X^degNext - 1)
		for i := 1; i <= int(math.Ceil(float64(r.redDeg)/float64(r.degNext))); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[0][j] = mod.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.modulus)
				r.buf.pIn[0][i*r.degNext+j] = 0
			}
		}

		// Compute pOut = pIn - pRem
		for i := 0; i < r.params.Degree(); i++ {
			pOut[i] = mod.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.modulus)
		}
	} else {
		copy(pOut, r.buf.pIn[0][:r.params.Degree()])
	}
}
