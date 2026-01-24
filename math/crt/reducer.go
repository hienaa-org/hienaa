package crt

import (
	"math"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// reducerBuffer is a buffer for [Reducer].
type reducerBuffer struct {
	// pIn is a buffer for the input polynomial.
	pIn [][]uint64
	// pQuo is a buffer for the quotient polynomial.
	pQuo [][]uint64
	// pRem is a buffer for the remainder polynomial.
	pRem [][]uint64
}

// newReducerBuffer creates a new [reducerBuffer].
func newReducerBuffer(lenAmbMod, in, quo, rem int) reducerBuffer {
	pIn := make([][]uint64, lenAmbMod)
	pQuo := make([][]uint64, lenAmbMod)
	pRem := make([][]uint64, lenAmbMod)

	for i := range pIn {
		pIn[i] = make([]uint64, in)
		pQuo[i] = make([]uint64, quo)
		pRem[i] = make([]uint64, rem)
	}

	return reducerBuffer{
		pIn:  pIn,
		pQuo: pQuo,
		pRem: pRem,
	}
}

// LongDivReducer reduces a polynomial modulo modulus polynomial using long division.
// In cases where the modulus polynomial is small or sparse, this might be more efficient than [Reducer].
type LongDivReducer struct {
	params dft.RingParameters
	mod    []*num.Modulus

	// maxRank is the maximum rank of the input polynomial.
	maxRank int

	// modPoly is the polynomial we target to reduce to.
	modPoly [][]uint64

	buf reducerBuffer
}

// NewLongDivReducer creates a new [LongDivReducer].
func NewLongDivReducer(maxRank int, mod []*num.Modulus, modPoly []int64) *LongDivReducer {
	if maxRank < len(modPoly)-1 {
		panic("maxRank must be greater than or equal to modPoly degree")
	} else if modPoly[len(modPoly)-1] != 1 {
		panic("modPoly must be monic")
	}

	modPolyRed := make([][]uint64, len(mod))
	for i := range mod {
		modPolyRed[i] = vec.Reduce(modPoly, mod[i])
	}

	return &LongDivReducer{
		params: dft.NewOtherParameters(modPoly),
		mod:    mod,

		maxRank: maxRank,
		modPoly: modPolyRed,

		buf: newReducerBuffer(1, maxRank, maxRank-len(modPoly)+1, maxRank),
	}
}

// quoRemTo computes quotient and remainder of p to pQuo, pRem with the idx-th modulus.
func (r *LongDivReducer) quoRemTo(pQuo, pRem, p []uint64, idx int) {
	clear(pQuo)
	clear(r.buf.pIn[0])
	copy(r.buf.pIn[0], p)

	for i := 0; i <= len(p)-len(r.modPoly[idx]); i++ {
		if r.buf.pIn[0][len(r.buf.pIn[0])-i-1] != 0 {
			pQuo[len(pQuo)-i-1] = r.buf.pIn[0][len(r.buf.pIn[0])-i-1]
			vec.ScalarMulSubTo(
				r.buf.pIn[0][len(r.buf.pIn[0])-i-len(r.modPoly[idx]):len(r.buf.pIn[0])-i],
				r.modPoly[idx],
				pQuo[len(pQuo)-i-1],
				r.mod[idx],
			)
		}
	}

	copy(pRem, r.buf.pIn[0][:r.params.Rank()])
}

// Reduce reduces p.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) Reduce(p *Poly) *Poly {
	pOut := NewPoly(r.params.Rank(), p.ModLen())
	r.ReduceTo(pOut, p)
	return pOut
}

// ReduceTo reduces p to pOut.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) ReduceTo(pOut, p *Poly) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.maxRank {
		panic("rank must be less than or equal to maxRank")
	} else if pOut.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	for i := range r.mod {
		r.quoRemTo(r.buf.pQuo[0], pOut.Coeffs[i], p.Coeffs[i], i)
	}
}

// Quotient returns p / modPoly.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) Quotient(p *Poly) *Poly {
	pOut := NewPoly(p.Rank()-r.params.Rank(), p.ModLen())
	r.QuotientTo(pOut, p)
	return pOut
}

// QuotientTo computes pOut = p / modPoly.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) QuotientTo(pOut, p *Poly) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.maxRank {
		panic("rank must be less than or equal to maxRank")
	} else if pOut.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	for i := range r.mod {
		r.quoRemTo(pOut.Coeffs[i], r.buf.pRem[0], p.Coeffs[i], i)
	}
}

// QuoRem returns p / modPoly and p % modPoly.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) QuoRem(p *Poly) (pQuo, pRem *Poly) {
	pQuo = NewPoly(p.Rank()-r.params.Rank(), p.ModLen())
	pRem = NewPoly(r.params.Rank(), p.ModLen())
	r.QuoRemTo(pQuo, pRem, p)
	return
}

// QuoRemTo computes pQuo = p / modPoly and pRem = p % modPoly.
//
// Panics when p, pQuo or pRem is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) QuoRemTo(pQuo, pRem, p *Poly) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.maxRank {
		panic("rank must be less than or equal to maxRank")
	} else if pQuo.ModLen() != len(r.mod) || pRem.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	for i := range r.mod {
		r.quoRemTo(pQuo.Coeffs[i], pRem.Coeffs[i], p.Coeffs[i], i)
	}
}

// Modulus returns the modulus.
func (r *LongDivReducer) Modulus() []*num.Modulus {
	return r.mod
}

// SubReducer returns a reducer for modulus of given indices.
func (r *LongDivReducer) SubReducer(idx ...int) *LongDivReducer {
	modCopy := make([]*num.Modulus, len(idx))
	modPolyCopy := make([][]uint64, len(idx))
	for i := range idx {
		modCopy[i] = r.mod[idx[i]]
		modPolyCopy[i] = r.modPoly[idx[i]]
	}

	return &LongDivReducer{
		params: r.params,
		mod:    modCopy,

		maxRank: r.maxRank,

		modPoly: modPolyCopy,

		buf: newReducerBuffer(1, r.maxRank, r.maxRank-len(modPolyCopy[0])+1, r.maxRank),
	}
}

// SafeCopy returns a thread-safe copy.
func (r *LongDivReducer) SafeCopy() *LongDivReducer {
	return &LongDivReducer{
		params: r.params,
		mod:    r.mod,

		maxRank: r.maxRank,

		modPoly: r.modPoly,

		buf: newReducerBuffer(1, r.maxRank, r.maxRank-len(r.modPoly[0])+1, r.maxRank),
	}
}

// CyclotomicReducer is an optimized [Reducer] for cyclotomic polynomial.
// In other words, it computes a(X) mod \Phi_m(X), where a(X) is at most degree m-1.
// It uses Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
// In this implementation, Q_sp = (X^m-1)/(X^(m/p)-1) for the smallest prime factor p.
type CyclotomicReducer struct {
	params dft.RingParameters
	mod    []*num.Modulus

	// leastFac is the smallest prime factor of the cyclotomic polynomial.
	leastFac int
	// isTrivial is true if the reduction is trivial.
	isTrivial bool

	// redDeg is the degree of the intermediate reducing polynomial Q_sp.
	redDeg int
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

	// ambModLen is the optimal length of ambient modulus for each non-NTT friendly modulus.
	ambModLen []int
	// ambMod is the ambient modulus.
	ambMod []*num.Modulus
	// embedder is the embedder from ambient modulus.
	embedder []*Embedder

	// diffDegNextAmbNTT is diffDegNextNTT for ambient modulus.
	diffDegNextAmbNTT []dft.Transformer
	// degNextAmbNTT is the degNextAmbNTT for ambient modulus.
	degNextAmbNTT []dft.Transformer

	// cycloPoly is the cyclotomic polynomial modulo the modulus.
	cycloPoly [][][]uint64
	// divPoly is rounding of a monomial over the cyclotomic polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/\Phi_m(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	divPoly [][][]uint64

	buf reducerBuffer
}

// NewCyclotomicReducer creates a new [CyclotomicReducer].
func NewCyclotomicReducer(params dft.RingParameters, mod []*num.Modulus) *CyclotomicReducer {
	cycloOrd, rank := params.CycloOrder(), params.Rank()
	primes, _ := num.Factor(cycloOrd)

	var redDeg, leastFactor int
	if cycloOrd%2 == 1 {
		leastFactor = primes[0]
		redDeg = cycloOrd - cycloOrd/leastFactor
	} else {
		leastFactor = primes[1]
		redDeg = cycloOrd/2 - (cycloOrd/2)/leastFactor
	}

	isTrivial := redDeg == rank

	var ambModLen []int
	var ambMod []*num.Modulus
	var embedder []*Embedder
	var diffDeg, diffDegNext, degNext int
	var diffDegNextNTT, degNextNTT []dft.Transformer
	var diffDegNextAmbNTT, degNextAmbNTT []dft.Transformer
	var cycloPoly, divPoly [][][]uint64

	if !isTrivial {
		degNext = num.NextProdPower(rank, []int{2})
		diffDeg = redDeg - rank
		diffDegNext = num.NextProdPower(2*diffDeg+1, []int{2})

		diffDegNextParams := dft.NewCyclicParameters(diffDegNext)
		degNextParams := dft.NewCyclicParameters(degNext)

		ambModLen = make([]int, len(mod))
		for i := range mod {
			if dft.IsNTTFriendly(diffDegNextParams, mod[i]) && dft.IsNTTFriendly(degNextParams, mod[i]) {
				continue
			}
			maxBits := num.Log2(max(diffDegNext, degNext)) + 2*num.Log2(mod[i].Value())
			ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
		}

		ambMod = dft.MustFindPrevNTTPrimes(params, num.MaxModulusBits, vec.Max(ambModLen))
		diffDegNextAmbNTT = make([]dft.Transformer, len(ambMod))
		degNextAmbNTT = make([]dft.Transformer, len(ambMod))
		for i := range ambMod {
			diffDegNextAmbNTT[i] = dft.NewTransformer(diffDegNextParams, ambMod[i])
			degNextAmbNTT[i] = dft.NewTransformer(degNextParams, ambMod[i])
		}

		embedder = make([]*Embedder, len(mod))
		diffDegNextNTT = make([]dft.Transformer, len(mod))
		degNextNTT = make([]dft.Transformer, len(mod))
		for i := range mod {
			if ambModLen[i] == 0 {
				diffDegNextNTT[i] = dft.NewTransformer(diffDegNextParams, mod[i])
				degNextNTT[i] = dft.NewTransformer(degNextParams, mod[i])
			} else {
				embedder[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
			}
		}

		cycloPolySigned := dft.CyclotomicPolynomial(cycloOrd)
		dividend := NewPoly(redDeg+1, 1)
		dividend.Coeffs[0][redDeg] = 1

		cycloPoly = make([][][]uint64, len(mod))
		divPoly = make([][][]uint64, len(mod))
		for i := range mod {
			divReducer := NewLongDivReducer(redDeg+1, []*num.Modulus{mod[i]}, cycloPolySigned)
			divPolyRef := divReducer.Quotient(dividend).Coeffs[0]
			if ambModLen[i] == 0 {
				cycloPoly[i] = [][]uint64{make([]uint64, degNext)}
				vec.ReduceTo(cycloPoly[i][0][:rank+1], cycloPolySigned, mod[i])
				divPoly[i] = [][]uint64{append(divPolyRef, make([]uint64, diffDegNext-len(divPolyRef))...)}
			} else {
				cycloPoly[i] = make([][]uint64, ambModLen[i])
				cycloPoly[i][0] = make([]uint64, degNext)
				vec.ReduceTo(cycloPoly[i][0][:rank+1], cycloPolySigned, mod[i])
				for j := 1; j < ambModLen[i]; j++ {
					cycloPoly[i][j] = make([]uint64, degNext)
					copy(cycloPoly[i][j], cycloPoly[i][0])
				}

				divPoly[i] = make([][]uint64, ambModLen[i])
				divPoly[i][0] = append(divPolyRef, make([]uint64, diffDegNext-len(divPolyRef))...)
				for j := 1; j < ambModLen[i]; j++ {
					divPoly[i][j] = make([]uint64, diffDegNext)
					copy(divPoly[i][j], divPoly[i][0])
				}
			}
		}

		for i := range mod {
			if ambModLen[i] == 0 {
				diffDegNextNTT[i].ForwardTo(divPoly[i][0], divPoly[i][0])
				degNextNTT[i].ForwardTo(cycloPoly[i][0], cycloPoly[i][0])
			} else {
				for j := 0; j < ambModLen[i]; j++ {
					diffDegNextAmbNTT[j].ForwardTo(divPoly[i][j], divPoly[i][j])
					degNextAmbNTT[j].ForwardTo(cycloPoly[i][j], cycloPoly[i][j])
				}
			}
		}
	}

	return &CyclotomicReducer{
		params: params,
		mod:    mod,

		leastFac:  leastFactor,
		isTrivial: isTrivial,

		redDeg:      redDeg,
		diffDeg:     diffDeg,
		diffDegNext: diffDegNext,
		degNext:     degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		embedder:  embedder,

		diffDegNextAmbNTT: diffDegNextAmbNTT,
		degNextAmbNTT:     degNextAmbNTT,

		cycloPoly: cycloPoly,
		divPoly:   divPoly,

		buf: newReducerBuffer(max(1, vec.Max(ambModLen)), cycloOrd, diffDegNext, degNext),
	}
}

// reduceTo reduces p to pOut with the idx-th modulus.
func (r *CyclotomicReducer) reduceTo(pOut, p []uint64, idx int) {
	cycloOrd, rank := r.params.CycloOrder(), r.params.Rank()

	copy(r.buf.pIn[0], p)
	if cycloOrd%2 == 1 {
		skip := cycloOrd / r.leastFac

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				r.buf.pIn[0][i*skip+j] = num.Sub(r.buf.pIn[0][i*skip+j], r.buf.pIn[0][cycloOrd-skip+j], r.mod[idx])
			}
			r.buf.pIn[0][cycloOrd-skip+j] = 0
		}
	} else {
		skip := (cycloOrd / 2) / r.leastFac

		for i := 0; i < cycloOrd/2; i++ {
			r.buf.pIn[0][i] = num.Sub(r.buf.pIn[0][i], r.buf.pIn[0][cycloOrd/2+i], r.mod[idx])
			r.buf.pIn[0][cycloOrd/2+i] = 0
		}

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				if i%2 == 0 {
					r.buf.pIn[0][i*skip+j] = num.Sub(r.buf.pIn[0][i*skip+j], r.buf.pIn[0][cycloOrd/2-skip+j], r.mod[idx])
				} else {
					r.buf.pIn[0][i*skip+j] = num.Add(r.buf.pIn[0][i*skip+j], r.buf.pIn[0][cycloOrd/2-skip+j], r.mod[idx])
				}
			}
			r.buf.pIn[0][cycloOrd/2-skip+j] = 0
		}
	}

	if !r.isTrivial {
		// pQuo = floor(pIn / X^deg)
		for i := 0; i < max(1, r.ambModLen[idx]); i++ {
			clear(r.buf.pQuo[i])
			for j := 0; j < r.diffDeg; j++ {
				r.buf.pQuo[i][j] = r.buf.pIn[0][rank+j]
			}
		}

		// pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
		if r.ambModLen[idx] == 0 {
			r.diffDegNextNTT[idx].ForwardTo(r.buf.pQuo[0], r.buf.pQuo[0])
			vec.MMulLazyTo(r.buf.pQuo[0], r.buf.pQuo[0], r.divPoly[idx][0], r.mod[idx])
			r.diffDegNextNTT[idx].InverseTo(r.buf.pQuo[0], r.buf.pQuo[0])
		} else {
			for i := 0; i < r.ambModLen[idx]; i++ {
				r.diffDegNextAmbNTT[i].ForwardTo(r.buf.pQuo[i], r.buf.pQuo[i])
				vec.MMulLazyTo(r.buf.pQuo[i], r.buf.pQuo[i], r.divPoly[idx][i], r.ambMod[i])
				r.diffDegNextAmbNTT[i].InverseTo(r.buf.pQuo[i], r.buf.pQuo[i])
			}
			r.embedder[idx].EmbedVecTo(r.buf.pQuo[:1], r.buf.pQuo[:r.ambModLen[idx]])
		}

		// pRem = floor(pQuo / X^diffDeg) % (X^degNext - 1)
		for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.diffDeg {
					break
				}
				r.buf.pQuo[0][r.diffDeg+j] = num.Add(r.buf.pQuo[0][r.diffDeg+j], r.buf.pQuo[0][r.diffDeg+i*r.degNext+j], r.mod[idx])
				r.buf.pQuo[0][r.diffDeg+i*r.degNext+j] = 0
			}
		}

		for i := 0; i < max(1, r.ambModLen[idx]); i++ {
			clear(r.buf.pRem[i])
			for j := 0; j < min(r.diffDeg, r.degNext); j++ {
				r.buf.pRem[i][j] = r.buf.pQuo[0][r.diffDeg+j]
			}
		}

		// pRem = pRem * cycloPoly % (X^degNext - 1)
		if r.ambModLen[idx] == 0 {
			r.degNextNTT[idx].ForwardTo(r.buf.pRem[0], r.buf.pRem[0])
			vec.MMulLazyTo(r.buf.pRem[0], r.buf.pRem[0], r.cycloPoly[idx][0], r.mod[idx])
			r.degNextNTT[idx].InverseTo(r.buf.pRem[0], r.buf.pRem[0])
		} else {
			for i := 0; i < r.ambModLen[idx]; i++ {
				r.degNextAmbNTT[i].ForwardTo(r.buf.pRem[i], r.buf.pRem[i])
				vec.MMulLazyTo(r.buf.pRem[i], r.buf.pRem[i], r.cycloPoly[idx][i], r.ambMod[i])
				r.degNextAmbNTT[i].InverseTo(r.buf.pRem[i], r.buf.pRem[i])
			}
			r.embedder[idx].EmbedVecTo(r.buf.pRem[:1], r.buf.pRem[:r.ambModLen[idx]])
		}

		// pIn = pIn % X^degNext - 1
		for i := 1; i <= num.DivCeil(r.redDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				r.buf.pIn[0][j] = num.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.mod[idx])
				r.buf.pIn[0][i*r.degNext+j] = 0
			}
		}

		// pOut = pIn - pRem
		for i := 0; i < rank; i++ {
			pOut[i] = num.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.mod[idx])
		}
	} else {
		copy(pOut, r.buf.pIn[0][:rank])
	}
}

// Reduce reduces p.
//
// Panics when p is in NTT form, or the rank of p is larger than CycloOrd.
func (r *CyclotomicReducer) Reduce(p *Poly) *Poly {
	pOut := NewPoly(r.params.Rank(), p.ModLen())
	r.ReduceTo(pOut, p)
	return pOut
}

// ReduceTo reduces p to pOut.
//
// Panics when p is in NTT form, or the rank of p is larger than CycloOrd.
func (r *CyclotomicReducer) ReduceTo(pOut, p *Poly) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.params.CycloOrder() {
		panic("rank must be less than or equal to cycloOrd")
	} else if pOut.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	for i := range r.mod {
		r.reduceTo(pOut.Coeffs[i], p.Coeffs[i], i)
	}
}

// Modulus returns the modulus.
func (r *CyclotomicReducer) Modulus() []*num.Modulus {
	return r.mod
}

// SubReducer returns a reducer for modulus of given indices.
func (r *CyclotomicReducer) SubReducer(idx ...int) *CyclotomicReducer {
	modCopy := make([]*num.Modulus, len(idx))

	for i := range idx {
		modCopy[i] = r.mod[idx[i]]
	}

	var ambModLenCopy []int
	var ambModCopy []*num.Modulus
	var embedderCopy []*Embedder
	var diffDegNextNTTCopy, degNextNTTCopy []dft.Transformer
	var diffDegNextAmbNTTCopy, degNextAmbNTTCopy []dft.Transformer
	var cycloPolyCopy, divPolyCopy [][][]uint64

	if !r.isTrivial {
		ambModLenCopy = make([]int, len(idx))
		cycloPolyCopy = make([][][]uint64, len(idx))
		divPolyCopy = make([][][]uint64, len(idx))
		for i := range idx {
			ambModLenCopy[i] = r.ambModLen[idx[i]]
			cycloPolyCopy[i] = r.cycloPoly[idx[i]]
			divPolyCopy[i] = r.divPoly[idx[i]]
		}

		embedderCopy = make([]*Embedder, len(idx))
		diffDegNextNTTCopy = make([]dft.Transformer, len(idx))
		degNextNTTCopy = make([]dft.Transformer, len(idx))
		for i := range idx {
			if r.embedder[idx[i]] != nil {
				embedderCopy[i] = r.embedder[idx[i]].SafeCopy()
			}
			if r.diffDegNextNTT[idx[i]] != nil {
				diffDegNextNTTCopy[i] = r.diffDegNextNTT[idx[i]].SafeCopy()
			}
			if r.degNextNTT[idx[i]] != nil {
				degNextNTTCopy[i] = r.degNextNTT[idx[i]].SafeCopy()
			}
		}

		maxAmbModLen := vec.Max(ambModLenCopy)
		diffDegNextAmbNTTCopy = make([]dft.Transformer, maxAmbModLen)
		degNextAmbNTTCopy = make([]dft.Transformer, maxAmbModLen)
		for i := 0; i < maxAmbModLen; i++ {
			diffDegNextAmbNTTCopy[i] = r.diffDegNextAmbNTT[i].SafeCopy()
			degNextAmbNTTCopy[i] = r.degNextAmbNTT[i].SafeCopy()
		}

		ambModCopy = r.ambMod[:maxAmbModLen]
	}

	return &CyclotomicReducer{
		params: r.params,
		mod:    modCopy,

		leastFac:  r.leastFac,
		isTrivial: r.isTrivial,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTTCopy,
		degNextNTT:     degNextNTTCopy,

		ambModLen: ambModLenCopy,
		ambMod:    ambModCopy,
		embedder:  embedderCopy,

		diffDegNextAmbNTT: diffDegNextAmbNTTCopy,
		degNextAmbNTT:     degNextAmbNTTCopy,

		cycloPoly: cycloPolyCopy,
		divPoly:   divPolyCopy,

		buf: newReducerBuffer(max(1, vec.Max(ambModLenCopy)), r.params.CycloOrder(), r.diffDegNext, r.degNext),
	}
}

// SafeCopy returns a thread-safe copy.
func (r *CyclotomicReducer) SafeCopy() *CyclotomicReducer {
	var embedderCopy []*Embedder
	var diffDegNextNTTCopy, degNextNTTCopy []dft.Transformer

	if !r.isTrivial {
		embedderCopy = make([]*Embedder, len(r.mod))
		diffDegNextNTTCopy = make([]dft.Transformer, len(r.mod))
		degNextNTTCopy = make([]dft.Transformer, len(r.mod))
		for i := range r.mod {
			if r.embedder[i] != nil {
				embedderCopy[i] = r.embedder[i].SafeCopy()
			}
			if r.diffDegNextNTT[i] != nil {
				diffDegNextNTTCopy[i] = r.diffDegNextNTT[i].SafeCopy()
			}
			if r.degNextNTT[i] != nil {
				degNextNTTCopy[i] = r.degNextNTT[i].SafeCopy()
			}
		}
	}

	diffDegNextAmbNTTCopy := make([]dft.Transformer, len(r.diffDegNextAmbNTT))
	for i := range r.diffDegNextAmbNTT {
		diffDegNextAmbNTTCopy[i] = r.diffDegNextAmbNTT[i].SafeCopy()
	}
	degNextAmbNTTCopy := make([]dft.Transformer, len(r.degNextAmbNTT))
	for i := range r.degNextAmbNTT {
		degNextAmbNTTCopy[i] = r.degNextAmbNTT[i].SafeCopy()
	}

	return &CyclotomicReducer{
		params: r.params,
		mod:    r.mod,

		leastFac:  r.leastFac,
		isTrivial: r.isTrivial,

		redDeg:      r.redDeg,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTTCopy,
		degNextNTT:     degNextNTTCopy,

		ambModLen: r.ambModLen,
		ambMod:    r.ambMod,
		embedder:  embedderCopy,

		diffDegNextAmbNTT: diffDegNextAmbNTTCopy,
		degNextAmbNTT:     degNextAmbNTTCopy,

		cycloPoly: r.cycloPoly,
		divPoly:   r.divPoly,

		buf: newReducerBuffer(max(1, vec.Max(r.ambModLen)), r.params.CycloOrder(), r.diffDegNext, r.degNext),
	}
}

// Reducer reduces a polynomial modulo modulus polynomial.
// In other words, it computes p(X) mod f(X), where the rank of p is at most maxRank.
// It uses the Barrett reduction for polynomial, from https://eprint.iacr.org/2015/818.
type Reducer struct {
	params dft.RingParameters
	mod    []*num.Modulus

	// maxRank is the maximum rank of the input polynomial.
	maxRank int
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

	// ambModLen is the optimal length of ambient modulus for each non-NTT friendly modulus.
	ambModLen []int
	// ambMod is the ambient modulus.
	ambMod []*num.Modulus
	// embedder is the embedder from ambient modulus.
	embedder []*Embedder

	// diffDegNextAmbNTT is diffDegNextNTT for ambient modulus.
	diffDegNextAmbNTT []dft.Transformer
	// degNextAmbNTT is the degNextAmbNTT for ambient modulus.
	degNextAmbNTT []dft.Transformer

	// modPoly is the polynomial we target to reduce to.
	modPoly [][][]uint64
	// divPoly is rounding of a monomial over the mod polynomial modulo the modulus.
	// Precisely, it is floor(X^d_qs/modPoly(X)) modulo the modulus, where d_qs is the degree of the quotient polynomial Q_sp.
	divPoly [][][]uint64

	buf reducerBuffer
}

// NewReducer creates a new [Reducer].
func NewReducer(maxRank int, mod []*num.Modulus, modPoly []int64) *Reducer {
	if maxRank < len(modPoly)-1 {
		panic("maxRank must be greater than or equal to modPoly degree")
	} else if modPoly[len(modPoly)-1] != 1 {
		panic("modPoly must be monic")
	}

	rank := len(modPoly) - 1

	degNext := num.NextProdPower(rank, []int{2})
	diffDeg := maxRank - rank - 1
	diffDegNext := num.NextProdPower(2*diffDeg+1, []int{2})

	diffDegNextParams := dft.NewCyclicParameters(diffDegNext)
	degNextParams := dft.NewCyclicParameters(degNext)

	ambParams := dft.NewCyclicParameters(2 * max(diffDegNext, degNext))
	ambModLen := make([]int, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(diffDegNextParams, mod[i]) && dft.IsNTTFriendly(degNextParams, mod[i]) {
			continue
		}
		maxBits := num.Log2(max(diffDegNext, degNext)) + 2*num.Log2(mod[i].Value())
		ambModLen[i] = int(math.Ceil(maxBits / num.MaxModulusBits))
	}

	ambMod := dft.MustFindPrevNTTPrimes(ambParams, num.MaxModulusBits, vec.Max(ambModLen))
	diffDegNextAmbNTT := make([]dft.Transformer, len(ambMod))
	degNextAmbNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		diffDegNextAmbNTT[i] = dft.NewTransformer(diffDegNextParams, ambMod[i])
		degNextAmbNTT[i] = dft.NewTransformer(degNextParams, ambMod[i])
	}

	embedder := make([]*Embedder, len(mod))
	diffDegNextNTT := make([]dft.Transformer, len(mod))
	degNextNTT := make([]dft.Transformer, len(mod))
	for i := range mod {
		if ambModLen[i] == 0 {
			diffDegNextNTT[i] = dft.NewTransformer(diffDegNextParams, mod[i])
			degNextNTT[i] = dft.NewTransformer(degNextParams, mod[i])
		} else {
			embedder[i] = NewEmbedder([]*num.Modulus{mod[i]}, ambMod[:ambModLen[i]])
		}
	}

	dividend := NewPoly(maxRank, 1)
	dividend.Coeffs[0][maxRank-1] = 1

	modPolyRed := make([][][]uint64, len(mod))
	divPoly := make([][][]uint64, len(mod))
	for i := range mod {
		divReducer := NewLongDivReducer(maxRank, []*num.Modulus{mod[i]}, modPoly)
		divPolyRef := divReducer.Quotient(dividend).Coeffs[0]

		if ambModLen[i] == 0 {
			modPolyRed[i] = [][]uint64{make([]uint64, degNext)}
			vec.ReduceTo(modPolyRed[i][0][:rank+1], modPoly, mod[i])

			divPoly[i] = [][]uint64{append(divPolyRef, make([]uint64, diffDegNext-len(divPolyRef))...)}
		} else {
			modPolyRed[i] = make([][]uint64, ambModLen[i])
			modPolyRed[i][0] = make([]uint64, degNext)
			vec.ReduceTo(modPolyRed[i][0][:rank+1], modPoly, mod[i])
			for j := 1; j < ambModLen[i]; j++ {
				modPolyRed[i][j] = make([]uint64, degNext)
				copy(modPolyRed[i][j], modPolyRed[i][0])
			}

			divPoly[i] = make([][]uint64, ambModLen[i])
			divPoly[i][0] = append(divPolyRef, make([]uint64, diffDegNext-len(divPolyRef))...)
			for j := 1; j < ambModLen[i]; j++ {
				divPoly[i][j] = make([]uint64, diffDegNext)
				copy(divPoly[i][j], divPoly[i][0])
			}
		}
	}

	for i := range mod {
		if ambModLen[i] == 0 {
			diffDegNextNTT[i].ForwardTo(divPoly[i][0], divPoly[i][0])
			degNextNTT[i].ForwardTo(modPolyRed[i][0], modPolyRed[i][0])
		} else {
			for j := 0; j < ambModLen[i]; j++ {
				diffDegNextAmbNTT[j].ForwardTo(divPoly[i][j], divPoly[i][j])
				degNextAmbNTT[j].ForwardTo(modPolyRed[i][j], modPolyRed[i][j])
			}
		}
	}

	return &Reducer{
		params: dft.NewOtherParameters(modPoly),
		mod:    mod,

		maxRank:     maxRank,
		diffDeg:     diffDeg,
		diffDegNext: diffDegNext,
		degNext:     degNext,

		diffDegNextNTT: diffDegNextNTT,
		degNextNTT:     degNextNTT,

		ambModLen: ambModLen,
		ambMod:    ambMod,
		embedder:  embedder,

		diffDegNextAmbNTT: diffDegNextAmbNTT,
		degNextAmbNTT:     degNextAmbNTT,

		modPoly: modPolyRed,
		divPoly: divPoly,

		buf: newReducerBuffer(max(1, vec.Max(ambModLen)), maxRank, diffDegNext, degNext),
	}
}

// reduceTo reduces p to pOut with the idx-th modulus.
func (r *Reducer) reduceTo(pOut, p []uint64, idx int) {
	rank := r.params.Rank()

	copy(r.buf.pIn[0], p)

	// pQuo = floor(pIn / X^deg)
	for i := 0; i < max(1, r.ambModLen[idx]); i++ {
		clear(r.buf.pQuo[i])
		for j := 0; j <= r.diffDeg; j++ {
			r.buf.pQuo[i][j] = r.buf.pIn[0][rank+j]
		}
	}

	// pQuo = pQuo * floor(X^(deg+diffDeg)/modPoly(X))
	if r.ambModLen[idx] == 0 {
		r.diffDegNextNTT[idx].ForwardTo(r.buf.pQuo[0], r.buf.pQuo[0])
		vec.MMulLazyTo(r.buf.pQuo[0], r.buf.pQuo[0], r.divPoly[idx][0], r.mod[idx])
		r.diffDegNextNTT[idx].InverseTo(r.buf.pQuo[0], r.buf.pQuo[0])
	} else {
		for i := 0; i < r.ambModLen[idx]; i++ {
			r.diffDegNextAmbNTT[i].ForwardTo(r.buf.pQuo[i], r.buf.pQuo[i])
			vec.MMulLazyTo(r.buf.pQuo[i], r.buf.pQuo[i], r.divPoly[idx][i], r.ambMod[i])
			r.diffDegNextAmbNTT[i].InverseTo(r.buf.pQuo[i], r.buf.pQuo[i])
		}
		r.embedder[idx].EmbedVecTo(r.buf.pQuo[:1], r.buf.pQuo[:r.ambModLen[idx]])
	}

	// pRem = floor(pQuo / X^diffDeg) % (X^degNext - 1)
	for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			r.buf.pQuo[0][r.diffDeg+j] = num.Add(r.buf.pQuo[0][r.diffDeg+j], r.buf.pQuo[0][r.diffDeg+i*r.degNext+j], r.mod[idx])
			r.buf.pQuo[0][r.diffDeg+i*r.degNext+j] = 0
		}
	}

	for i := 0; i < max(1, r.ambModLen[idx]); i++ {
		clear(r.buf.pRem[i])
		for j := 0; j < min(r.diffDeg+1, r.degNext); j++ {
			r.buf.pRem[i][j] = r.buf.pQuo[0][r.diffDeg+j]
		}
	}

	// pRem = pRem * modPoly % (X^degNext - 1)
	if r.ambModLen[idx] == 0 {
		r.degNextNTT[idx].ForwardTo(r.buf.pRem[0], r.buf.pRem[0])
		vec.MMulLazyTo(r.buf.pRem[0], r.buf.pRem[0], r.modPoly[idx][0], r.mod[idx])
		r.degNextNTT[idx].InverseTo(r.buf.pRem[0], r.buf.pRem[0])
	} else {
		for i := 0; i < r.ambModLen[idx]; i++ {
			r.degNextAmbNTT[i].ForwardTo(r.buf.pRem[i], r.buf.pRem[i])
			vec.MMulLazyTo(r.buf.pRem[i], r.buf.pRem[i], r.modPoly[idx][i], r.ambMod[i])
			r.degNextAmbNTT[i].InverseTo(r.buf.pRem[i], r.buf.pRem[i])
		}
		r.embedder[idx].EmbedVecTo(r.buf.pRem[:1], r.buf.pRem[:r.ambModLen[idx]])
	}

	// pIn = pIn % (X^degNext - 1)
	for i := 1; i <= num.DivCeil(r.maxRank-1, r.degNext); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j >= r.maxRank {
				break
			}
			r.buf.pIn[0][j] = num.Add(r.buf.pIn[0][j], r.buf.pIn[0][i*r.degNext+j], r.mod[idx])
			r.buf.pIn[0][i*r.degNext+j] = 0
		}
	}

	// pOut = pIn - pRem
	for i := 0; i < rank; i++ {
		pOut[i] = num.Sub(r.buf.pIn[0][i], r.buf.pRem[0][i], r.mod[idx])
	}
}

// Reduce reduces p.
//
// Panics when p is in NTT form, or the rank of p is larger than MaxRank.
func (r *Reducer) Reduce(p *Poly) *Poly {
	pOut := NewPoly(r.params.Rank(), p.ModLen())
	r.ReduceTo(pOut, p)
	return pOut
}

// ReduceTo reduces p to pOut.
//
// Panics when p is in NTT form, or the rank of p is larger than MaxRank.
func (r *Reducer) ReduceTo(pOut, p *Poly) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.maxRank {
		panic("rank must be less than or equal to maxRank")
	} else if pOut.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	for i := range r.mod {
		r.reduceTo(pOut.Coeffs[i], p.Coeffs[i], i)
	}
}

// Params returns the ring parameters.
func (r *Reducer) Params() dft.RingParameters {
	return r.params
}

// MaxRank returns the maximum possible rank of the input polynomial.
func (r *Reducer) MaxRank() int {
	return r.maxRank
}

// Modulus returns the modulus.
func (r *Reducer) Modulus() []*num.Modulus {
	return r.mod
}

// SubReducer returns a reducer for modulus of given indices.
func (r *Reducer) SubReducer(idx ...int) *Reducer {
	modCopy := make([]*num.Modulus, len(idx))
	ambModLenCopy := make([]int, len(idx))
	modPolyCopy := make([][][]uint64, len(idx))
	divPolyCopy := make([][][]uint64, len(idx))
	for i := range idx {
		modCopy[i] = r.mod[idx[i]]
		ambModLenCopy[i] = r.ambModLen[idx[i]]
		modPolyCopy[i] = r.modPoly[idx[i]]
		divPolyCopy[i] = r.divPoly[idx[i]]
	}

	embedderCopy := make([]*Embedder, len(idx))
	diffDegNextNTTCopy := make([]dft.Transformer, len(idx))
	degNextNTTCopy := make([]dft.Transformer, len(idx))
	for i := range idx {
		if r.embedder[idx[i]] != nil {
			embedderCopy[i] = r.embedder[idx[i]].SafeCopy()
		}
		if r.diffDegNextNTT[idx[i]] != nil {
			diffDegNextNTTCopy[i] = r.diffDegNextNTT[idx[i]].SafeCopy()
		}
		if r.degNextNTT[idx[i]] != nil {
			degNextNTTCopy[i] = r.degNextNTT[idx[i]].SafeCopy()
		}
	}

	maxAmbModLen := vec.Max(ambModLenCopy)
	diffDegNextAmbNTTCopy := make([]dft.Transformer, maxAmbModLen)
	degNextAmbNTTCopy := make([]dft.Transformer, maxAmbModLen)
	for i := 0; i < maxAmbModLen; i++ {
		diffDegNextAmbNTTCopy[i] = r.diffDegNextAmbNTT[i].SafeCopy()
		degNextAmbNTTCopy[i] = r.degNextAmbNTT[i].SafeCopy()
	}

	return &Reducer{
		params: r.params,
		mod:    modCopy,

		maxRank:     r.maxRank,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTTCopy,
		degNextNTT:     degNextNTTCopy,

		ambModLen: ambModLenCopy,
		ambMod:    r.ambMod[:maxAmbModLen],
		embedder:  embedderCopy,

		diffDegNextAmbNTT: diffDegNextAmbNTTCopy,
		degNextAmbNTT:     degNextAmbNTTCopy,

		modPoly: modPolyCopy,
		divPoly: divPolyCopy,

		buf: newReducerBuffer(max(1, maxAmbModLen), r.maxRank, r.diffDegNext, r.degNext),
	}

}

// SafeCopy returns a thread-safe copy.
func (r *Reducer) SafeCopy() *Reducer {
	embedderCopy := make([]*Embedder, len(r.mod))
	diffDegNextNTTCopy := make([]dft.Transformer, len(r.mod))
	degNextNTTCopy := make([]dft.Transformer, len(r.mod))
	for i := range r.mod {
		if r.embedder[i] != nil {
			embedderCopy[i] = r.embedder[i].SafeCopy()
		}
		if r.diffDegNextNTT[i] != nil {
			diffDegNextNTTCopy[i] = r.diffDegNextNTT[i].SafeCopy()
		}
		if r.degNextNTT[i] != nil {
			degNextNTTCopy[i] = r.degNextNTT[i].SafeCopy()
		}
	}

	diffDegNextAmbNTTCopy := make([]dft.Transformer, len(r.diffDegNextAmbNTT))
	for i := range r.diffDegNextAmbNTT {
		diffDegNextAmbNTTCopy[i] = r.diffDegNextAmbNTT[i].SafeCopy()
	}
	degNextAmbNTTCopy := make([]dft.Transformer, len(r.degNextAmbNTT))
	for i := range r.degNextAmbNTT {
		degNextAmbNTTCopy[i] = r.degNextAmbNTT[i].SafeCopy()
	}

	return &Reducer{
		params: r.params,
		mod:    r.mod,

		maxRank:     r.maxRank,
		diffDeg:     r.diffDeg,
		diffDegNext: r.diffDegNext,
		degNext:     r.degNext,

		diffDegNextNTT: diffDegNextNTTCopy,
		degNextNTT:     degNextNTTCopy,

		ambModLen: r.ambModLen,
		ambMod:    r.ambMod,
		embedder:  embedderCopy,

		diffDegNextAmbNTT: diffDegNextAmbNTTCopy,
		degNextAmbNTT:     degNextAmbNTTCopy,

		modPoly: r.modPoly,
		divPoly: r.divPoly,

		buf: newReducerBuffer(max(1, vec.Max(r.ambModLen)), r.maxRank, r.diffDegNext, r.degNext),
	}
}
