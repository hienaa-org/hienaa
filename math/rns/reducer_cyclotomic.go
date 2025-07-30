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
	// quoPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/\Phi_m(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quoPoly [][]uint64

	buf reducerBuffer
}

// newCyclotomicReducerAnyModulus creates a new [cyclotomicReducerAnyModulus].
func newCyclotomicReducerAnyModulus(params RingParameters, modulus *mod.Modulus) *cyclotomicReducerAnyModulus {
	cycloDeg, deg := uint64(params.CycloDegree()), uint64(params.Degree())

	factors := num.Factor(cycloDeg)
	leastFactor := cycloDeg
	for key := range factors {
		if key < leastFactor {
			leastFactor = key
		}
	}
	redDeg := cycloDeg - cycloDeg/leastFactor

	isPrimePower := redDeg == deg

	var ambModulus []*mod.Modulus
	var embedder *Embedder
	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT []singleTransformer
	var cycloPoly, quoPoly [][]uint64
	var buf reducerBuffer

	if !isPrimePower {
		degNext = num.NextProdPower(uint64(deg), []uint64{2})
		diffDeg = redDeg - deg
		diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

		lenAmbMod := int(math.Ceil((2*math.Log2(float64(modulus.Value())) + math.Log2(float64(max(degNext, diffDegNext)))) / mod.MaxModulusBits))
		ambModulus = FindPrevNTTPrimes(params, 61, lenAmbMod)
		embedder = NewEmbedder(ambModulus, []*mod.Modulus{modulus})

		degNextParams := NewCyclicParameters(int(degNext))
		degNextNTT = make([]singleTransformer, lenAmbMod)
		for i := range degNextNTT {
			degNextNTT[i] = newCyclicPow235Transformer(degNextParams, ambModulus[i])
		}

		diffDegNextParams := NewCyclicParameters(int(diffDegNext))
		diffDegNextNTT = make([]singleTransformer, lenAmbMod)
		for i := range diffDegNextNTT {
			diffDegNextNTT[i] = newCyclicPow235Transformer(diffDegNextParams, ambModulus[i])
		}

		cycloPoly = make([][]uint64, lenAmbMod)
		cycloPolySigned := cyclotomicPolynomial(cycloDeg)
		for i := range cycloPoly {
			cycloPoly[i] = make([]uint64, degNext)
			for j := range cycloPolySigned {
				cycloPoly[i][j] = reduceInt(cycloPolySigned[j], ambModulus[i])
			}
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quoPoly = make([][]uint64, lenAmbMod)
		for i := range quoPoly {
			quoPoly[i] = append(quotientPolynomial(dividend, cycloPoly[i][:deg+1], ambModulus[i]), make([]uint64, int(diffDegNext)-len(quoPoly[i]))...)
		}

		for i := 0; i < lenAmbMod; i++ {
			degNextNTT[i].nttInPlace(cycloPoly[i])
			diffDegNextNTT[i].nttInPlace(quoPoly[i])
		}

		buf = newReducerBuffer(lenAmbMod, int(cycloDeg), int(diffDegNext), int(degNext))
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

		cycloPoly: cycloPoly,
		quoPoly:   quoPoly,

		buf: buf,
	}
}

func (r *cyclotomicReducerAnyModulus) reduceTo(pOut, p []uint64) {
	copy(r.buf.pIn[0], p)

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
		// Compute pQuo = floor(pIn/X^deg)
		for i := 0; i < len(r.buf.pQuo); i++ {
			clear(r.buf.pQuo[i])
			for j := 0; j < r.diffDeg; j++ {
				r.buf.pQuo[i][j] = r.buf.pIn[0][deg+j]
			}
		}

		// Compute pQuo = pQuo × floor(X^(deg+diffDeg)/\Phi_m(X))
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

func (r *cyclotomicReducerAnyModulus) safeCopy() reducer {
	diffDegNextNTT := make([]singleTransformer, len(r.diffDegNextNTT))
	for i := range diffDegNextNTT {
		diffDegNextNTT[i] = r.diffDegNextNTT[i].safeCopy()
	}

	degNextNTT := make([]singleTransformer, len(r.degNextNTT))
	for i := range degNextNTT {
		degNextNTT[i] = r.degNextNTT[i].safeCopy()
	}

	var embedder *Embedder
	var buf reducerBuffer
	if !r.isPrimePow {
		embedder = r.embedder.SafeCopy()
		buf = newReducerBuffer(len(r.ambModulus), r.params.cycloDegree, r.diffDegNext, r.degNext)
	}

	return &cyclotomicReducerAnyModulus{
		leastFac:   r.leastFac,
		isPrimePow: r.isPrimePow,

		params:     r.params,
		modulus:    r.modulus,
		ambModulus: r.ambModulus,
		embedder:   embedder,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly: r.cycloPoly,
		quoPoly:   r.quoPoly,

		buf: buf,
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
	// quoPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/\Phi_m(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quoPoly []uint64

	buf reducerBuffer
}

// newCyclotomicReducerNTTModulus creates a new [cyclotomicReducerNTTModulus].
func newCyclotomicReducerNTTModulus(params RingParameters, modulus *mod.Modulus) *cyclotomicReducerNTTModulus {
	cycloDeg, deg := uint64(params.CycloDegree()), uint64(params.Degree())

	factors := num.Factor(cycloDeg)
	leastFactor := cycloDeg
	for key := range factors {
		if key < leastFactor {
			leastFactor = key
		}
	}
	redDeg := cycloDeg - cycloDeg/leastFactor

	isPrimePower := redDeg == deg

	var diffDeg, diffDegNext, degNext uint64
	var diffDegNextNTT, degNextNTT singleTransformer
	var cycloPoly, quoPoly []uint64
	var buf reducerBuffer

	if !isPrimePower {
		degNext = num.NextProdPower(uint64(deg), []uint64{2})
		diffDeg = redDeg - deg
		diffDegNext = num.NextProdPower(2*diffDeg+1, []uint64{2})

		degNextParams := RingParameters{0, int(degNext), Cyclic}
		degNextNTT = newCyclicPow235Transformer(degNextParams, modulus)

		diffDegNextParams := RingParameters{0, int(diffDegNext), Cyclic}
		diffDegNextNTT = newCyclicPow235Transformer(diffDegNextParams, modulus)

		cycloPoly = make([]uint64, degNext)
		cycloPolySigned := cyclotomicPolynomial(cycloDeg)
		for i := range cycloPolySigned {
			cycloPoly[i] = reduceInt(cycloPolySigned[i], modulus)
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quoPoly = append(quotientPolynomial(dividend, cycloPoly[:deg+1], modulus), make([]uint64, int(diffDegNext)-len(quoPoly))...)

		degNextNTT.nttInPlace(cycloPoly)
		diffDegNextNTT.nttInPlace(quoPoly)

		buf = newReducerBuffer(1, int(cycloDeg), int(diffDegNext), int(degNext))
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

		cycloPoly: cycloPoly,
		quoPoly:   quoPoly,

		buf: buf,
	}
}

func (r *cyclotomicReducerNTTModulus) reduceTo(pOut, p []uint64) {
	copy(r.buf.pIn[0], p)

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
		// Compute pQuo = floor(pIn/X^deg)
		clear(r.buf.pQuo[0])
		for i := 0; i < r.diffDeg; i++ {
			r.buf.pQuo[0][i] = r.buf.pIn[0][deg+i]
		}

		// Compute pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
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

		// Compute pRem = pRem * quoPoly (mod X^degNext - 1)
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

func (r *cyclotomicReducerNTTModulus) safeCopy() reducer {
	var diffDegNextNTT, degNextNTT singleTransformer
	var buf reducerBuffer
	if !r.isPrimePow {
		diffDegNextNTT = r.diffDegNextNTT.safeCopy()
		degNextNTT = r.degNextNTT.safeCopy()
		buf = newReducerBuffer(1, r.params.cycloDegree, r.diffDegNext, r.degNext)
	}

	return &cyclotomicReducerNTTModulus{
		leastFac:   r.leastFac,
		isPrimePow: r.isPrimePow,

		params:  r.params,
		modulus: r.modulus,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly: r.cycloPoly,
		quoPoly:   r.quoPoly,

		buf: buf,
	}
}
