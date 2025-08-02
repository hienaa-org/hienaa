package crt

import (
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Reducer reduces polynomials.
type Reducer struct {
	reducer reducer
}

// NewReducer creates a new [Reducer] for arbitrary polynomials.
func NewReducer(maxDeg int, mod *num.Modulus, modPoly []uint64) *Reducer {
	return &Reducer{
		reducer: newReducer(maxDeg, modPoly, mod),
	}
}

// Reduce reduces p and returns the result.
func (r *Reducer) Reduce(p *Poly) *Poly {
	if p.isNTT {
		panic("Reduce: cannot reduce NTT polynomials")
	}

	pOut := NewPoly(p.Rank(), p.ModLen())
	for i := range p.Coeffs {
		r.reducer.reduceTo(pOut.Coeffs[i], p.Coeffs[i])
	}
	return pOut
}

// ReduceTo reduces p to pOut.
// p and pOut must be Standard form.
func (r *Reducer) ReduceTo(pOut, p *Poly) {
	if p.isNTT || pOut.isNTT {
		panic("ReduceTo: cannot reduce NTT polynomials")
	}

	for i := range pOut.Coeffs {
		r.reducer.reduceTo(pOut.Coeffs[i], p.Coeffs[i])
	}
}

// SafeCopy returns a thread-safe copy.
func (r *Reducer) SafeCopy() *Reducer {
	return &Reducer{
		reducer: r.reducer.safeCopy(),
	}
}

// reducer is an interface for reducing single polynomials.
type reducer interface {
	// reduceTo reduces the polynomial p to pOut.
	reduceTo(pOut, p []uint64)
	// safeCopy returns a thread-safe copy.
	safeCopy() reducer
}

// newReducer creates a new [reducer].
func newReducer(maxDeg int, modPoly []uint64, mod *num.Modulus) reducer {
	deg := len(modPoly) - 1
	degNext := int(num.NextProdPower(uint64(deg), []uint64{2}))
	diffDegNext := int(num.NextProdPower(2*uint64(maxDeg-deg)+1, []uint64{2}))

	if dft.IsNTTFriendly(dft.NewCyclicParameters(max(degNext, diffDegNext)), mod) {
		return newReducerNTTModulus(maxDeg, mod, modPoly)
	}
	return newReducerAnyModulus(maxDeg, mod, modPoly)
}

// reducerBuffer is a buffer for [cyclotomicReducerNTTModulus].
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
