package crt

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// LongDivReducer reduces a polynomial modulo modulus polynomial using long division.
// In cases where the modulus polynomial is small or sparse, this might be more efficient than [Reducer].
type LongDivReducer struct {
	params dft.RingParameters
	mod    []*num.Modulus

	// maxRank is the maximum rank of the input polynomial.
	maxRank int

	// modPoly is the polynomial we target to reduce to.
	modPoly [][]uint64

	pool *sync.Pool
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

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, max(maxRank, maxRank-len(modPoly)+1))
				return &v
			},
		},
	}
}

// quoRemTo computes quotient and remainder of p to pQuo, pRem with the idx-th modulus.
func (r *LongDivReducer) quoRemTo(pQuo, pRem, p []uint64, idx int) {
	pInPtr := r.pool.Get().(*[]uint64)
	pIn := (*pInPtr)[:r.maxRank]
	defer r.pool.Put(pInPtr)

	clear(pQuo)
	clear(pIn)
	copy(pIn, p)

	for i := 0; i <= len(p)-len(r.modPoly[idx]); i++ {
		if pIn[len(pIn)-i-1] != 0 {
			pQuo[len(pQuo)-i-1] = pIn[len(pIn)-i-1]
			vec.MulSubScalarTo(
				pIn[len(pIn)-i-len(r.modPoly[idx]):len(pIn)-i],
				r.modPoly[idx],
				pQuo[len(pQuo)-i-1],
				r.mod[idx],
			)
		}
	}

	copy(pRem, pIn[:r.params.Rank()])
}

// Reduce reduces p.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) Reduce(p *Element) *Element {
	pOut := NewPoly(r.params.Rank(), p.ModLen())
	r.ReduceTo(pOut, p)
	return pOut
}

// ReduceTo reduces p to pOut.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) ReduceTo(pOut, p *Element) {
	if p.Type() != TypePoly {
		panic("input(s) must be polynomial")
	} else if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.maxRank {
		panic("rank must be less than or equal to maxRank")
	} else if pOut.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	pQuoPtr := r.pool.Get().(*[]uint64)
	pQuo := (*pQuoPtr)[:r.maxRank-len(r.modPoly[0])+1]
	defer r.pool.Put(pQuoPtr)

	for i := range r.mod {
		r.quoRemTo(pQuo, pOut.Coeffs[i], p.Coeffs[i], i)
	}
}

// Quotient returns p / modPoly.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) Quotient(p *Element) *Element {
	pOut := NewPoly(p.Rank()-r.params.Rank(), p.ModLen())
	r.QuotientTo(pOut, p)
	return pOut
}

// QuotientTo computes pOut = p / modPoly.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) QuotientTo(pOut, p *Element) {
	if p.IsNTT {
		panic("input(s) must be in standard form")
	} else if p.Rank() > r.maxRank {
		panic("rank must be less than or equal to maxRank")
	} else if pOut.ModLen() != len(r.mod) || p.ModLen() != len(r.mod) {
		panic("input(s) not consistent")
	}

	pRemPtr := r.pool.Get().(*[]uint64)
	pRem := (*pRemPtr)[:r.params.Rank()]
	defer r.pool.Put(pRemPtr)

	for i := range r.mod {
		r.quoRemTo(pOut.Coeffs[i], pRem, p.Coeffs[i], i)
	}
}

// QuoRem returns p / modPoly and p % modPoly.
//
// Panics when p is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) QuoRem(p *Element) (pQuo, pRem *Element) {
	pQuo = NewPoly(p.Rank()-r.params.Rank(), p.ModLen())
	pRem = NewPoly(r.params.Rank(), p.ModLen())
	r.QuoRemTo(pQuo, pRem, p)
	return
}

// QuoRemTo computes pQuo = p / modPoly and pRem = p % modPoly.
//
// Panics when p, pQuo or pRem is in NTT form, or the rank of p is larger than maxRank.
func (r *LongDivReducer) QuoRemTo(pQuo, pRem, p *Element) {
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

		pool: r.pool,
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

	pool *sync.Pool
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

		maxBits := make([]float64, len(mod))
		for i := range mod {
			if dft.IsNTTFriendly(diffDegNextParams, mod[i]) && dft.IsNTTFriendly(degNextParams, mod[i]) {
				continue
			}
			maxBits[i] = float64(num.Log2(max(diffDegNext, degNext)) + 2*num.Log2(mod[i].Value()))
		}

		ambMod = dft.MustFindAmbientPrimes(params, vec.Max(maxBits))
		diffDegNextAmbNTT = make([]dft.Transformer, len(ambMod))
		degNextAmbNTT = make([]dft.Transformer, len(ambMod))
		for i := range ambMod {
			diffDegNextAmbNTT[i] = dft.NewTransformer(diffDegNextParams, ambMod[i])
			degNextAmbNTT[i] = dft.NewTransformer(degNextParams, ambMod[i])
		}

		ambModLen = make([]int, len(mod))
		for i := range mod {
			if maxBits[i] == 0 {
				continue
			}
			currBits := 0.0
			for j := range ambMod {
				currBits += num.Log2(ambMod[j].Value())
				if currBits >= maxBits[i] {
					ambModLen[i] = j + 1
					break
				}
			}
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

		cycloPolySigned := params.ModulusPoly()
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

		pool: &sync.Pool{
			New: func() any {
				p := make([][]uint64, max(1, vec.Max(ambModLen)))
				for i := range p {
					p[i] = make([]uint64, max(cycloOrd, diffDegNext, degNext))
				}
				return &p
			},
		},
	}
}

// reduceTo reduces p to pOut with the idx-th modulus.
func (r *CyclotomicReducer) reduceTo(pOut, p []uint64, idx int) {
	cycloOrd, rank := r.params.CycloOrder(), r.params.Rank()

	pInPtr := r.pool.Get().(*[][]uint64)
	pIn := (*pInPtr)[0][:cycloOrd]
	defer r.pool.Put(pInPtr)

	copy(pIn, p)
	if cycloOrd%2 == 1 {
		skip := cycloOrd / r.leastFac

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				pIn[i*skip+j] = num.Sub(pIn[i*skip+j], pIn[cycloOrd-skip+j], r.mod[idx])
			}
			pIn[cycloOrd-skip+j] = 0
		}
	} else {
		skip := (cycloOrd / 2) / r.leastFac

		for i := 0; i < cycloOrd/2; i++ {
			pIn[i] = num.Sub(pIn[i], pIn[cycloOrd/2+i], r.mod[idx])
			pIn[cycloOrd/2+i] = 0
		}

		for j := 0; j < skip; j++ {
			for i := 0; i < r.leastFac-1; i++ {
				if i%2 == 0 {
					pIn[i*skip+j] = num.Sub(pIn[i*skip+j], pIn[cycloOrd/2-skip+j], r.mod[idx])
				} else {
					pIn[i*skip+j] = num.Add(pIn[i*skip+j], pIn[cycloOrd/2-skip+j], r.mod[idx])
				}
			}
			pIn[cycloOrd/2-skip+j] = 0
		}
	}

	if !r.isTrivial {
		pQuoPtr := r.pool.Get().(*[][]uint64)
		pQuo := *pQuoPtr
		for i := range pQuo {
			pQuo[i] = pQuo[i][:r.diffDegNext]
		}
		defer r.pool.Put(pQuoPtr)

		// pQuo = floor(pIn / X^deg)
		for i := 0; i < max(1, r.ambModLen[idx]); i++ {
			clear(pQuo[i])
			for j := 0; j < r.diffDeg; j++ {
				pQuo[i][j] = pIn[rank+j]
			}
		}

		// pQuo = pQuo * floor(X^(deg+diffDeg)/\Phi_m(X))
		if r.ambModLen[idx] == 0 {
			r.diffDegNextNTT[idx].ForwardTo(pQuo[0], pQuo[0])
			vec.MMulLazyTo(pQuo[0], pQuo[0], r.divPoly[idx][0], r.mod[idx])
			r.diffDegNextNTT[idx].InverseTo(pQuo[0], pQuo[0])
		} else {
			for i := 0; i < r.ambModLen[idx]; i++ {
				r.diffDegNextAmbNTT[i].ForwardTo(pQuo[i], pQuo[i])
				vec.MMulLazyTo(pQuo[i], pQuo[i], r.divPoly[idx][i], r.ambMod[i])
				r.diffDegNextAmbNTT[i].InverseTo(pQuo[i], pQuo[i])
			}
			r.embedder[idx].EmbedVecTo(pQuo[:1], pQuo[:r.ambModLen[idx]])
		}

		pRemPtr := r.pool.Get().(*[][]uint64)
		pRem := *pRemPtr
		for i := range pRem {
			pRem[i] = pRem[i][:r.degNext]
		}
		defer r.pool.Put(pRemPtr)

		// pRem = floor(pQuo / X^diffDeg) % (X^degNext - 1)
		for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.diffDeg {
					break
				}
				pQuo[0][r.diffDeg+j] = num.Add(pQuo[0][r.diffDeg+j], pQuo[0][r.diffDeg+i*r.degNext+j], r.mod[idx])
				pQuo[0][r.diffDeg+i*r.degNext+j] = 0
			}
		}

		for i := 0; i < max(1, r.ambModLen[idx]); i++ {
			clear(pRem[i])
			for j := 0; j < min(r.diffDeg, r.degNext); j++ {
				pRem[i][j] = pQuo[0][r.diffDeg+j]
			}
		}

		// pRem = pRem * cycloPoly % (X^degNext - 1)
		if r.ambModLen[idx] == 0 {
			r.degNextNTT[idx].ForwardTo(pRem[0], pRem[0])
			vec.MMulLazyTo(pRem[0], pRem[0], r.cycloPoly[idx][0], r.mod[idx])
			r.degNextNTT[idx].InverseTo(pRem[0], pRem[0])
		} else {
			for i := 0; i < r.ambModLen[idx]; i++ {
				r.degNextAmbNTT[i].ForwardTo(pRem[i], pRem[i])
				vec.MMulLazyTo(pRem[i], pRem[i], r.cycloPoly[idx][i], r.ambMod[i])
				r.degNextAmbNTT[i].InverseTo(pRem[i], pRem[i])
			}
			r.embedder[idx].EmbedVecTo(pRem[:1], pRem[:r.ambModLen[idx]])
		}

		// pIn = pIn % X^degNext - 1
		for i := 1; i <= num.DivCeil(r.redDeg, r.degNext); i++ {
			for j := 0; j < r.degNext; j++ {
				if i*r.degNext+j > r.redDeg {
					break
				}
				pIn[j] = num.Add(pIn[j], pIn[i*r.degNext+j], r.mod[idx])
				pIn[i*r.degNext+j] = 0
			}
		}

		// pOut = pIn - pRem
		for i := 0; i < rank; i++ {
			pOut[i] = num.Sub(pIn[i], pRem[0][i], r.mod[idx])
		}
	} else {
		copy(pOut, pIn[:rank])
	}
}

// Reduce reduces p.
//
// Panics when p is in NTT form, or the rank of p is larger than CycloOrd.
func (r *CyclotomicReducer) Reduce(p *Element) *Element {
	pOut := NewPoly(r.params.Rank(), p.ModLen())
	r.ReduceTo(pOut, p)
	return pOut
}

// ReduceTo reduces p to pOut.
//
// Panics when p is in NTT form, or the rank of p is larger than CycloOrd.
func (r *CyclotomicReducer) ReduceTo(pOut, p *Element) {
	if p.Type() != TypePoly {
		panic("input(s) must be polynomial")
	} else if p.IsNTT {
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
			embedderCopy[i] = r.embedder[idx[i]]
			diffDegNextNTTCopy[i] = r.diffDegNextNTT[idx[i]]
			degNextNTTCopy[i] = r.degNextNTT[idx[i]]
		}

		maxAmbModLen := vec.Max(ambModLenCopy)
		ambModCopy = r.ambMod[:maxAmbModLen]
		diffDegNextAmbNTTCopy = r.diffDegNextAmbNTT[:maxAmbModLen]
		degNextAmbNTTCopy = r.degNextAmbNTT[:maxAmbModLen]
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

		pool: r.pool,
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

	pool *sync.Pool
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

	maxBits := make([]float64, len(mod))
	for i := range mod {
		if dft.IsNTTFriendly(diffDegNextParams, mod[i]) && dft.IsNTTFriendly(degNextParams, mod[i]) {
			continue
		}
		maxBits[i] = float64(num.Log2(max(diffDegNext, degNext)) + 2*num.Log2(mod[i].Value()))
	}

	ambMod := dft.MustFindAmbientPrimes(ambParams, vec.Max(maxBits))
	diffDegNextAmbNTT := make([]dft.Transformer, len(ambMod))
	degNextAmbNTT := make([]dft.Transformer, len(ambMod))
	for i := range ambMod {
		diffDegNextAmbNTT[i] = dft.NewTransformer(diffDegNextParams, ambMod[i])
		degNextAmbNTT[i] = dft.NewTransformer(degNextParams, ambMod[i])
	}

	ambModLen := make([]int, len(mod))
	for i := range mod {
		if maxBits[i] == 0 {
			continue
		}
		currBits := 0.0
		for j := range ambMod {
			currBits += num.Log2(ambMod[j].Value())
			if currBits >= maxBits[i] {
				ambModLen[i] = j + 1
				break
			}
		}
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

		pool: &sync.Pool{
			New: func() any {
				p := make([][]uint64, max(1, vec.Max(ambModLen)))
				for i := range p {
					p[i] = make([]uint64, max(maxRank, diffDegNext, degNext))
				}
				return &p
			},
		},
	}
}

// reduceTo reduces p to pOut with the idx-th modulus.
func (r *Reducer) reduceTo(pOut, p []uint64, idx int) {
	rank := r.params.Rank()

	pInPtr := r.pool.Get().(*[][]uint64)
	pIn := (*pInPtr)[0][:r.maxRank]
	defer r.pool.Put(pInPtr)

	copy(pIn, p)

	pQuoPtr := r.pool.Get().(*[][]uint64)
	pQuo := *pQuoPtr
	for i := range pQuo {
		pQuo[i] = pQuo[i][:r.diffDegNext]
	}
	defer r.pool.Put(pQuoPtr)

	// pQuo = floor(pIn / X^deg)
	for i := 0; i < max(1, r.ambModLen[idx]); i++ {
		clear(pQuo[i])
		for j := 0; j <= r.diffDeg; j++ {
			pQuo[i][j] = pIn[rank+j]
		}
	}

	// pQuo = pQuo * floor(X^(deg+diffDeg)/modPoly(X))
	if r.ambModLen[idx] == 0 {
		r.diffDegNextNTT[idx].ForwardTo(pQuo[0], pQuo[0])
		vec.MMulLazyTo(pQuo[0], pQuo[0], r.divPoly[idx][0], r.mod[idx])
		r.diffDegNextNTT[idx].InverseTo(pQuo[0], pQuo[0])
	} else {
		for i := 0; i < r.ambModLen[idx]; i++ {
			r.diffDegNextAmbNTT[i].ForwardTo(pQuo[i], pQuo[i])
			vec.MMulLazyTo(pQuo[i], pQuo[i], r.divPoly[idx][i], r.ambMod[i])
			r.diffDegNextAmbNTT[i].InverseTo(pQuo[i], pQuo[i])
		}
		r.embedder[idx].EmbedVecTo(pQuo[:1], pQuo[:r.ambModLen[idx]])
	}

	pRemPtr := r.pool.Get().(*[][]uint64)
	pRem := *pRemPtr
	for i := range pRem {
		pRem[i] = pRem[i][:r.degNext]
	}
	defer r.pool.Put(pRemPtr)

	// pRem = floor(pQuo / X^diffDeg) % (X^degNext - 1)
	for i := 1; i <= num.DivCeil(r.diffDeg, r.degNext); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j > r.diffDeg {
				break
			}
			pQuo[0][r.diffDeg+j] = num.Add(pQuo[0][r.diffDeg+j], pQuo[0][r.diffDeg+i*r.degNext+j], r.mod[idx])
			pQuo[0][r.diffDeg+i*r.degNext+j] = 0
		}
	}

	for i := 0; i < max(1, r.ambModLen[idx]); i++ {
		clear(pRem[i])
		for j := 0; j < min(r.diffDeg+1, r.degNext); j++ {
			pRem[i][j] = pQuo[0][r.diffDeg+j]
		}
	}

	// pRem = pRem * modPoly % (X^degNext - 1)
	if r.ambModLen[idx] == 0 {
		r.degNextNTT[idx].ForwardTo(pRem[0], pRem[0])
		vec.MMulLazyTo(pRem[0], pRem[0], r.modPoly[idx][0], r.mod[idx])
		r.degNextNTT[idx].InverseTo(pRem[0], pRem[0])
	} else {
		for i := 0; i < r.ambModLen[idx]; i++ {
			r.degNextAmbNTT[i].ForwardTo(pRem[i], pRem[i])
			vec.MMulLazyTo(pRem[i], pRem[i], r.modPoly[idx][i], r.ambMod[i])
			r.degNextAmbNTT[i].InverseTo(pRem[i], pRem[i])
		}
		r.embedder[idx].EmbedVecTo(pRem[:1], pRem[:r.ambModLen[idx]])
	}

	// pIn = pIn % (X^degNext - 1)
	for i := 1; i <= num.DivCeil(r.maxRank-1, r.degNext); i++ {
		for j := 0; j < r.degNext; j++ {
			if i*r.degNext+j >= r.maxRank {
				break
			}
			pIn[j] = num.Add(pIn[j], pIn[i*r.degNext+j], r.mod[idx])
			pIn[i*r.degNext+j] = 0
		}
	}

	// pOut = pIn - pRem
	for i := 0; i < rank; i++ {
		pOut[i] = num.Sub(pIn[i], pRem[0][i], r.mod[idx])
	}
}

// Reduce reduces p.
//
// Panics when p is in NTT form, or the rank of p is larger than MaxRank.
func (r *Reducer) Reduce(p *Element) *Element {
	pOut := NewPoly(r.params.Rank(), p.ModLen())
	r.ReduceTo(pOut, p)
	return pOut
}

// ReduceTo reduces p to pOut.
//
// Panics when p is in NTT form, or the rank of p is larger than MaxRank.
func (r *Reducer) ReduceTo(pOut, p *Element) {
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
		embedderCopy[i] = r.embedder[idx[i]]
		diffDegNextNTTCopy[i] = r.diffDegNextNTT[idx[i]]
		degNextNTTCopy[i] = r.degNextNTT[idx[i]]
	}

	maxAmbModLen := vec.Max(ambModLenCopy)
	diffDegNextAmbNTTCopy := make([]dft.Transformer, maxAmbModLen)
	degNextAmbNTTCopy := make([]dft.Transformer, maxAmbModLen)
	for i := 0; i < maxAmbModLen; i++ {
		diffDegNextAmbNTTCopy[i] = r.diffDegNextAmbNTT[i]
		degNextAmbNTTCopy[i] = r.degNextAmbNTT[i]
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

		pool: r.pool,
	}
}
