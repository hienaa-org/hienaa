package crt

import (
	"math"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/internal/dftops"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type cyclotomicReducerNTTModulus struct {
	params dft.RingParameters
	mod    *num.Modulus

	reducer *dftops.CyclotomicReducerNTTModulus
}

func newCyclotomicReducerNTTModulus(params dft.RingParameters, mod *num.Modulus) *cyclotomicReducerNTTModulus {
	return &cyclotomicReducerNTTModulus{
		params: params,
		mod:    mod,

		reducer: dftops.NewCyclotomicReducerNTTModulus(params.CycloOrder(), params.Rank(), mod),
	}
}

func (r *cyclotomicReducerNTTModulus) reduceTo(pOut, p []uint64) {
	r.reducer.ReduceTo(pOut, p)
}

func (r *cyclotomicReducerNTTModulus) safeCopy() reducer {
	return &cyclotomicReducerNTTModulus{
		params: r.params,
		mod:    r.mod,

		reducer: r.reducer.SafeCopy(),
	}
}

type cyclotomicReducerAnyModulus struct {
	params dft.RingParameters
	mod    *num.Modulus

	// ambMod is the modulus for the arbitrary modulus reducer.
	ambMod []*num.Modulus
	// embedder is the embedder for the arbitrary modulus reducer.
	embedder *Embedder

	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isPrimePow is true if the cyclotomic polynomial is a prime power.
	isPrimePow bool

	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
	// diffDeg is the degree of the difference between the cyclotomic polynomial and the intermediate reducing polynomial Q_sp.
	diffDeg int
	// diffDegNext is the smallest power of 2 that is greater than 2*diffDeg+1.
	diffDegNext int
	// degNext is the smallest power of 2 that is greater than the degree of the cyclotomic polynomial.
	degNext int

	// diffDegNextNTT is the NTT transformer for degree diffDegNext.
	diffDegNextNTT []*dftops.CyclicPow2Transformer
	// degNextNTT is the NTT transformer for degree degNext.
	degNextNTT []*dftops.CyclicPow2Transformer

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly [][]uint64
	// quoPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/\Phi_m(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	quoPoly [][]uint64

	buf reducerBuffer
}

// newCyclotomicReducerAnyModulus creates a new [cyclotomicReducerAnyModulus].
func newCyclotomicReducerAnyModulus(params dft.RingParameters, mod *num.Modulus) *cyclotomicReducerAnyModulus {
	cycloOrd, rank := params.CycloOrder(), params.Rank()

	primes, _ := num.Factor(cycloOrd)
	leastFactor := cycloOrd
	for _, p := range primes {
		if p < leastFactor {
			leastFactor = p
		}
	}
	redDeg := cycloOrd - cycloOrd/leastFactor

	isPrimePower := redDeg == rank

	var ambMod []*num.Modulus
	var embedder *Embedder
	var diffDeg, diffDegNext, degNext int
	var diffDegNextNTT, degNextNTT []*dftops.CyclicPow2Transformer
	var cycloPoly, quoPoly [][]uint64
	var buf reducerBuffer

	if !isPrimePower {
		degNext = num.NextProdPower(rank, []int{2})
		diffDeg = redDeg - rank
		diffDegNext = num.NextProdPower(2*diffDeg+1, []int{2})

		maxBits := 2*num.Log2(mod.Value()) + num.Log2(max(degNext, diffDegNext))
		lenAmbMod := int(math.Ceil(maxBits / num.MaxModulusBits))
		ambMod = dft.FindPrevNTTPrimes(params, num.MaxModulusBits, lenAmbMod)
		embedder = NewEmbedder([]*num.Modulus{mod}, ambMod)

		degNextNTT = make([]*dftops.CyclicPow2Transformer, lenAmbMod)
		for i := range degNextNTT {
			degNextNTT[i] = dftops.NewCyclicPow2Transformer(degNext, ambMod[i])
		}

		diffDegNextNTT = make([]*dftops.CyclicPow2Transformer, lenAmbMod)
		for i := range diffDegNextNTT {
			diffDegNextNTT[i] = dftops.NewCyclicPow2Transformer(diffDegNext, ambMod[i])
		}

		cycloPoly = make([][]uint64, lenAmbMod)
		cycloPolySigned := dftops.CyclotomicPolynomial(params.CycloOrder())
		for i := range cycloPoly {
			cycloPoly[i] = make([]uint64, degNext)
			for j := range cycloPolySigned {
				cycloPoly[i][j] = dftops.ReduceInt(cycloPolySigned[j], ambMod[i])
			}
		}

		dividend := make([]uint64, redDeg+1)
		dividend[redDeg] = 1
		quoPoly = make([][]uint64, lenAmbMod)
		for i := range quoPoly {
			quoPoly[i] = dftops.Quotient(dividend, cycloPoly[i][:rank+1], ambMod[i])
			quoPoly[i] = append(quoPoly[i], make([]uint64, diffDegNext-len(quoPoly[i]))...)
		}

		for i := 0; i < lenAmbMod; i++ {
			degNextNTT[i].ForwardInPlace(cycloPoly[i])
			diffDegNextNTT[i].ForwardInPlace(quoPoly[i])
		}

		buf = newReducerBuffer(lenAmbMod, cycloOrd, diffDegNext, degNext)
	} else {
		buf = newReducerBuffer(1, 0, 0, 0)
	}

	return &cyclotomicReducerAnyModulus{
		leastFac:   int(leastFactor),
		isPrimePow: isPrimePower,

		params:   params,
		mod:      mod,
		ambMod:   ambMod,
		embedder: embedder,

		redDeg:      redDeg,
		diffDeg:     diffDeg,
		diffDegNext: diffDegNext,
		degNext:     degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		cycloPoly: cycloPoly,
		quoPoly:   quoPoly,

		buf: buf,
	}
}

func (r *cyclotomicReducerAnyModulus) reduceTo(pOut, p []uint64) {
	copy(r.buf.pIn[0], p)

	cycloOrd := r.params.CycloOrder()
	rank := r.params.Rank()
	skip := cycloOrd / r.leastFac

	for j := 0; j < skip; j++ {
		for i := 0; i < r.leastFac-1; i++ {
			r.buf.pIn[0][i*skip+j] = num.Sub(r.buf.pIn[0][i*skip+j], r.buf.pIn[0][cycloOrd-skip+j], r.mod)
		}
		r.buf.pIn[0][cycloOrd-skip+j] = 0
	}

	if !r.isPrimePow {
		// Compute pQuo = floor(pIn/X^deg)
		for i := 0; i < len(r.buf.pQuo); i++ {
			clear(r.buf.pQuo[i])
			for j := 0; j < r.diffDeg; j++ {
				r.buf.pQuo[i][j] = r.buf.pIn[0][rank+j]
			}
		}

		// Compute pQuo = pQuo × floor(X^(deg+diffDeg)/\Phi_m(X))
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
			vec.MMulTo(r.buf.pRem[i], r.buf.pRem[i], r.cycloPoly[i], r.ambMod[i])
			r.degNextNTT[i].InverseInPlace(r.buf.pRem[i])
		}
		r.embedder.EmbedVecTo(r.buf.pRem[0:1], r.buf.pRem)

		// Compute pIn = pIn (mod X^degNext - 1)
		for i := 1; i <= num.DivCeil(r.redDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[0][j] = num.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.mod)
				r.buf.pIn[0][i*r.degNext+j] = 0
			}
		}

		// Compute pOut = pIn - pRem
		for i := 0; i < r.params.Rank(); i++ {
			pOut[i] = num.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.mod)
		}
	} else {
		copy(pOut, r.buf.pIn[0][:rank])
	}
}

func (r *cyclotomicReducerAnyModulus) safeCopy() reducer {
	var embedder *Embedder
	var buf reducerBuffer
	if !r.isPrimePow {
		embedder = r.embedder.SafeCopy()
		buf = newReducerBuffer(len(r.ambMod), r.params.CycloOrder(), r.diffDegNext, r.degNext)
	} else {
		buf = newReducerBuffer(1, 0, 0, 0)
	}

	return &cyclotomicReducerAnyModulus{
		leastFac:   r.leastFac,
		isPrimePow: r.isPrimePow,

		params:   r.params,
		mod:      r.mod,
		ambMod:   r.ambMod,
		embedder: embedder,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: r.diffDegNextNTT,
		degNextNTT:     r.degNextNTT,

		cycloPoly: r.cycloPoly,
		quoPoly:   r.quoPoly,

		buf: buf,
	}
}
