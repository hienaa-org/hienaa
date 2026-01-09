package crt

import (
	"math/big"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// polyBaseEvaluator is the base evaluator for all rings.
type polyBaseEvaluator struct {
	params  dft.RingParameters
	mod     []*num.Modulus
	modPoly []int64

	isNTTFriendly []bool
	ntt           []dft.Transformer
}

// newPolyBaseEvaluator creates a new [polyBaseEvaluator].
func newPolyBaseEvaluator(params dft.RingParameters, mod []*num.Modulus, modPoly []int64) polyBaseEvaluator {
	isNTTFriendly := make([]bool, len(mod))
	ntt := make([]dft.Transformer, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
		if isNTTFriendly[i] {
			ntt[i] = dft.NewTransformer(params, mod[i])
		}
	}

	return polyBaseEvaluator{
		params:  params,
		mod:     mod,
		modPoly: modPoly,

		isNTTFriendly: isNTTFriendly,
		ntt:           ntt,
	}
}

// NewPoly returns a new [Poly].
func (e *polyBaseEvaluator) NewPoly() *Poly {
	return NewPolyCustom(e.params.Rank(), len(e.mod), false)
}

// NewNTTPoly returns a new [Poly] in NTT form.
func (e *polyBaseEvaluator) NewNTTPoly() *Poly {
	return NewPolyCustom(e.params.Rank(), len(e.mod), true)
}

// NewPolyCustom creates a new [Poly] with the given parameters.
func (e *polyBaseEvaluator) NewPolyCustom(isNTT bool) *Poly {
	return NewPolyCustom(e.params.Rank(), len(e.mod), isNTT)
}

// Params returns the ring parameters.
func (e *polyBaseEvaluator) Params() dft.RingParameters {
	return e.params
}

// Modulus returns the modulus.
func (e *polyBaseEvaluator) Modulus() []*num.Modulus {
	return e.mod
}

// ModulusPoly returns the quotient polynomial of the ring.
func (e *polyBaseEvaluator) ModulusPoly() []int64 {
	return e.modPoly
}

// FwdNTT returns FwdNTT(p).
func (e *polyBaseEvaluator) FwdNTT(p *Poly) *Poly {
	pOut := e.NewPoly()
	e.FwdNTTTo(pOut, p)
	return pOut
}

// FwdNTTTo computes pOut = NTT(p).
func (e *polyBaseEvaluator) FwdNTTTo(pOut, p *Poly) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("NTTTo: inputs not consistent")
	case p.IsNTT:
		panic("NTTTo: already in NTT form")
	}

	for i := range e.ntt {
		if e.ntt[i] != nil {
			e.ntt[i].ForwardTo(pOut.Coeffs[i], p.Coeffs[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
		}
	}

	pOut.IsNTT = true
}

// InvNTT returns InvNTT(p).
func (e *polyBaseEvaluator) InvNTT(p *Poly) *Poly {
	pOut := e.NewPoly()
	e.InvNTTTo(pOut, p)
	return pOut
}

// InvNTTTo computes pOut = InvNTT(p).
func (e *polyBaseEvaluator) InvNTTTo(pOut, p *Poly) {
	switch {
	case !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p):
		panic("InvNTTTo: inputs not consistent")
	case !p.IsNTT:
		panic("InvNTTTo: already in Standard form")
	}

	for i := range e.ntt {
		if e.ntt[i] != nil {
			e.ntt[i].InverseTo(pOut.Coeffs[i], p.Coeffs[i])
		} else {
			copy(pOut.Coeffs[i], p.Coeffs[i])
		}
	}

	pOut.IsNTT = false
}

// Add returns p0 + p1.
func (e *polyBaseEvaluator) Add(p0, p1 *Poly) *Poly {
	pOut := e.NewPoly()
	e.AddTo(pOut, p0, p1)
	return pOut
}

// AddTo computes pOut = p0 + p1.
func (e *polyBaseEvaluator) AddTo(pOut, p0, p1 *Poly) {
	if !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1) {
		panic("AddTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.AddTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
	}

	pOut.IsNTT = p0.IsNTT
}

// Sub returns p0 - p1.
func (e *polyBaseEvaluator) Sub(p0, p1 *Poly) *Poly {
	pOut := e.NewPoly()
	e.SubTo(pOut, p0, p1)
	return pOut
}

// SubTo computes pOut = p0 - p1.
func (e *polyBaseEvaluator) SubTo(pOut, p0, p1 *Poly) {
	if !isTernaryToOperable(e.params.Rank(), len(e.mod), pOut, p0, p1) {
		panic("SubTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.SubTo(pOut.Coeffs[i], p0.Coeffs[i], p1.Coeffs[i], e.mod[i])
	}

	pOut.IsNTT = p0.IsNTT
}

// Neg returns -p.
func (e *polyBaseEvaluator) Neg(p *Poly) *Poly {
	pOut := NewPolyCustom(e.params.Rank(), len(e.mod), p.IsNTT)
	e.NegTo(pOut, p)
	return pOut
}

// NegTo computes pOut = -p.
func (e *polyBaseEvaluator) NegTo(pOut, p *Poly) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) {
		panic("NegTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.NegTo(pOut.Coeffs[i], p.Coeffs[i], e.mod[i])
	}

	pOut.IsNTT = p.IsNTT
}

// ScalarMul returns p * c.
func (e *polyBaseEvaluator) ScalarMul(p *Poly, c Scalar) *Poly {
	pOut := e.NewPoly()
	e.ScalarMulTo(pOut, p, c)
	return pOut
}

// ScalarMulTo computes pOut = p * c.
func (e *polyBaseEvaluator) ScalarMulTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarMulTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.ScalarMulTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
	}

	pOut.IsNTT = p.IsNTT
}

// ScalarMulAddTo computes pOut += p * c.
func (e *polyBaseEvaluator) ScalarMulAddTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarMulAddTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.ScalarMulAddTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
	}

	pOut.IsNTT = p.IsNTT
}

// ScalarMulSubTo computes pOut -= p * c.
func (e *polyBaseEvaluator) ScalarMulSubTo(pOut, p *Poly, c Scalar) {
	if !isBinaryToOperable(e.params.Rank(), len(e.mod), pOut, p) || !isScalarToOperable(len(e.mod), c) {
		panic("ScalarMulAddTo: inputs not consistent")
	}

	for i := range e.mod {
		vec.ScalarMulSubTo(pOut.Coeffs[i], p.Coeffs[i], c[i], e.mod[i])
	}

	pOut.IsNTT = p.IsNTT
}

// AsBig returns p as *[big.Int] vector.
func (e *polyBaseEvaluator) AsBig(p *Poly) []*big.Int {
	if !isConsistent(e.params.Rank(), len(e.mod), p) {
		panic("AsBig: input not consistent")
	} else if p.IsNTT {
		panic("input is in NTT form")
	}

	modBig := make([]*big.Int, len(e.mod))
	modProd := big.NewInt(1)
	for i := range e.mod {
		modBig[i] = new(big.Int).SetUint64(e.mod[i].Value())
		modProd.Mul(modProd, modBig[i])
	}
	modProdHalf := new(big.Int).Rsh(modProd, 1)

	gadget := make([]*big.Int, len(e.mod))
	for i := range modBig {
		qStar := new(big.Int).Div(modProd, modBig[i])
		qStarInv := new(big.Int).ModInverse(qStar, modBig[i])
		gadget[i] = new(big.Int).Mul(qStar, qStarInv)
		gadget[i].Mod(gadget[i], modProd)
	}

	pBig := make([]*big.Int, e.params.Rank())
	for j := 0; j < e.params.Rank(); j++ {
		pBig[j] = big.NewInt(0)
		for i := range e.mod {
			c := new(big.Int).SetUint64(p.Coeffs[i][j])
			c.Mul(c, gadget[i])
			pBig[j].Add(pBig[j], c)
		}
		pBig[j].Mod(pBig[j], modProd)
		if pBig[j].Cmp(modProdHalf) > 0 {
			pBig[j].Sub(pBig[j], modProd)
		}
	}

	return pBig
}

func (e *polyBaseEvaluator) subEvaluator(idx ...int) polyBaseEvaluator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	nttCopy := make([]dft.Transformer, len(idx))
	for i := range idx {
		modCopy[i] = e.mod[idx[i]]
		isNTTFriendlyCopy[i] = e.isNTTFriendly[idx[i]]
		if e.ntt[idx[i]] != nil {
			nttCopy[i] = e.ntt[idx[i]].SafeCopy()
		}
	}

	return polyBaseEvaluator{
		params:  e.params,
		mod:     modCopy,
		modPoly: e.modPoly,

		isNTTFriendly: isNTTFriendlyCopy,
		ntt:           nttCopy,
	}
}

func (e *polyBaseEvaluator) safeCopy() polyBaseEvaluator {
	nttCopy := make([]dft.Transformer, len(e.ntt))
	for i := range e.ntt {
		if e.ntt[i] != nil {
			nttCopy[i] = e.ntt[i].SafeCopy()
		}
	}

	return polyBaseEvaluator{
		params:  e.params,
		mod:     e.mod,
		modPoly: e.modPoly,

		isNTTFriendly: e.isNTTFriendly,
		ntt:           nttCopy,
	}
}
