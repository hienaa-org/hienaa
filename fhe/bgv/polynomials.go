package bgv

import (
	"math"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	"github.com/hienaa-org/hienaa/fhe/polyutils"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
)

type Polynomial struct {
	coeffs   map[int]Plaintext
	polyType polyutils.PolynomialType
}

// NewPolynomial creates a new polynomial from the given coefficients and type.
func NewPolynomial(c []uint64, pType polyutils.PolynomialType, msgMod *num.Modulus) *Polynomial {
	if len(c) == 0 {
		panic("polynomial must have at least one coefficient")
	}

	switch pType {
	case polyutils.Monomial:
		coeffs := make(map[int]Plaintext)
		for i, c := range c {
			if num.Reduce(c, msgMod) != 0 {
				coeffs[i] = NewScalarFrom(c, msgMod)
			}
		}
		return &Polynomial{
			coeffs:   coeffs,
			polyType: pType,
		}
	default:
		panic("unsupported polynomial type")
	}
}

// Degree returns the degree of the polynomial.
func (p *Polynomial) Degree() int {
	max := 0
	for k := range p.coeffs {
		if k > max {
			max = k
		}
	}
	return max
}

// Type returns the type of the polynomial.
func (p *Polynomial) Type() polyutils.PolynomialType {
	return p.polyType
}

// Coeffs returns the coefficients of the polynomial.
func (p *Polynomial) Coeffs() map[int]Plaintext {
	return p.coeffs
}

// TODO: Move to the right place.
// TODO: Change EvaluatePoly and evalRecurse to a in-place algorithm.
// EvaluatePoly evaluates the polynomial at the given ciphertext.
func (op *Operator) EvaluatePoly(p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) *Ciphertext {
	// Handle the edge case.
	if p.Degree() == 0 {
		res := NewCiphertextCustom(op.params.Rank(), ct.ModLen(), true)
		op.AddPlainTo(res, res, p.coeffs[0], true)
		return res
	}

	// First compute the basis.
	basis := make(map[int]*Ciphertext)
	basisLazy := make(map[int]*rlwe.Vector)

	// Hyperparameters
	ctLen := ct.ModLen()
	deg := p.Degree()
	maxLevel := int(math.Ceil(num.Log2(deg + 1)))
	maxDeg := (1 << maxLevel) - 1
	babyLevel := (maxLevel + 1) / 2
	babyDeg := 1 << babyLevel

	// BABYSTEP COMPUTATION
	bs := make(map[int]*Ciphertext)
	bsLazy := make(map[int]*rlwe.Vector)

	// Compute the power-of-two monomials.
	basis[1] = op.ctPool.Get()
	defer op.ctPool.Put(basis[1])
	basis[1] = basis[1].WithModLen(ctLen)
	if ct.IsNTT() {
		basis[1].CopyFrom(ct)
	} else {
		op.FwdNTTTo(basis[1], ct)
	}
	bs[1] = basis[1]

	for i := 1; i < maxLevel; i++ {
		basis[1<<i] = op.ctPool.Get()
		defer op.ctPool.Put(basis[1<<i])
		basis[1<<i] = basis[1<<i].WithModLen(ctLen)
		bs[1<<i] = basis[1<<i].WithModLen(ctLen)

		op.MulTo(bs[1<<i], bs[1<<(i-1)], bs[1<<(i-1)], rlk, true)
	}

	// Compute other monomials.
	for i := 1; i < babyDeg; i++ {
		if num.IsPowerOfTwo(i) {
			continue
		}

		// Check if i-th basis is needed in the babystep generation.
		check := false
		for j := int(math.Ceil(num.Log2(i))); j < babyLevel; j++ {
			for k := 0; k < deg; k += babyDeg {
				_, checkIdx := p.coeffs[k+1<<j+i]
				check = check || checkIdx
			}
		}

		if check {
			basis[i] = op.ctPool.Get()
			defer op.ctPool.Put(basis[i])
			basis[i] = basis[i].WithModLen(ctLen)
			bs[i] = basis[i].WithModLen(ctLen)

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			op.MulTo(bs[i], bs[idx1], bs[idx2], rlk, true)

			continue
		}

		// Check if i-th basis is directly needed in the evaluation.
		check = false
		for j := 0; j < deg; j += babyDeg {
			_, checkIdx := p.coeffs[j+i]
			check = check || checkIdx
		}

		if check {
			basisLazy[i] = op.vPool.Get()
			defer op.vPool.Put(basisLazy[i])
			basisLazy[i] = basisLazy[i].WithModLen(ctLen, 0)

			basis[i] = op.ctPool.Get()
			defer op.ctPool.Put(basis[i])
			basis[i] = basis[i].WithModLen(ctLen)
			bs[i] = basis[i].WithModLen(ctLen)

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			tarLen, auxIdx, auxMod := op.noise.getAuxMod(bs[idx1], bs[idx2])
			bsLazy[i] = basisLazy[i].WithModLen(tarLen, 0)

			op.tensorTo(bsLazy[i], bs[idx1], bs[idx2], tarLen, auxIdx, auxMod, true)

			bs[i] = bs[i].WithModLen(tarLen)
			bs[i].noise = op.noise.tensorTo(bs[idx1], bs[idx2], tarLen, auxIdx, auxMod)
		}
	}

	// Scale the basis to allow the lazy relin BSGS algorithm.
	for i := range basis {
		if bsLazy[i] != nil {
			for j := 0; j < 3; j++ {
				op.rlweOp.PlainOperator().ScaleTo(basisLazy[i].Value[j], bsLazy[i].Value[j], ctLen, true)
			}
			op.noise.ModSwitchTo(basis[i], bs[i], ctLen)
		} else {
			op.ModSwitchTo(basis[i], bs[i], ctLen, true)
		}
	}

	// Evaluate the polynomial.
	ctOut, vOut := op.evalRecurse(0, maxDeg, p, basis, basisLazy, rlk)

	if vOut != nil {
		ctOut = ctOut.WithModLen(vOut.BaseModLen())
		op.rlweOp.RelinTo(ctOut.Value, vOut, rlk, true)
		ctOut.noise += op.noise.noise.GadgetProd(ctLen)
	}

	op.RescaleTo(ctOut, ctOut, isNTT)

	return ctOut
}

// TODO: Change to an in-place algorithm.
// TODO: Optimise later using lazy relin BSGS algorithm.
func (op *Operator) evalRecurse(lo, hi int, p *Polynomial, basis map[int]*Ciphertext, basisLazy map[int]*rlwe.Vector, rlk *rlwe.RelinKey) (*Ciphertext, *rlwe.Vector) {
	// Hyperparameters.
	deg := p.Degree()
	maxLevel := int(math.Ceil(num.Log2(deg + 1)))
	maxDeg := (1 << maxLevel) - 1
	babyLevel := (maxLevel + 1) / 2
	babyDeg := 1 << babyLevel
	ctLen := basis[1].ModLen()

	// Current degree.
	curDeg := hi - lo + 1
	if !num.IsPowerOfTwo(curDeg) {
		panic("Current degree is not a power of two")
	}

	pOp := op.rlweOp.PlainOperator()
	if curDeg <= babyDeg && (hi != maxDeg || curDeg == 2) {
		// Babystep computation.
		if lo > deg {
			return nil, nil
		}

		var ctOut *Ciphertext
		var vOut *rlwe.Vector

		sc := op.ePool.Get(crt.TypeScalar)
		defer op.ePool.Put(sc)
		sc = sc.WithModLen(ctLen, 0)

		if val, check := p.coeffs[lo]; check {
			ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
			op.AddPlainTo(ctOut, ctOut, val, true)
		}

		for i := 1; i <= min(hi, deg)-lo; i++ {
			if val, check := p.coeffs[lo+i]; check {
				if ctOut == nil {
					ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
				}

				op.intOp.Encoder().EncodeTo(sc, val, false)

				if basisLazy[i] != nil {
					if vOut == nil {
						vOut = rlwe.NewVectorCustom(op.params.Rank(), ctLen, 0, 3, true)
						vOut.Value[0].CopyFrom(ctOut.Value.Body)
						vOut.Value[1].CopyFrom(ctOut.Value.Mask)
					}

					pOp.MulAddTo(vOut.Value[0], basisLazy[i].Value[0], sc)
					pOp.MulAddTo(vOut.Value[1], basisLazy[i].Value[1], sc)
					pOp.MulAddTo(vOut.Value[2], basisLazy[i].Value[2], sc)
				} else if vOut == nil {
					pOp.MulAddTo(ctOut.Value.Body, basis[i].Value.Body, sc)
					pOp.MulAddTo(ctOut.Value.Mask, basis[i].Value.Mask, sc)
				} else {
					pOp.MulAddTo(vOut.Value[0], basis[i].Value.Body, sc)
					pOp.MulAddTo(vOut.Value[1], basis[i].Value.Mask, sc)
				}

				switch op.noise.estimType {
				case heint.VarianceType:
					ctOut.noise += basis[i].noise * float64(val[0]) * float64(val[0])
				case heint.WorstCaseType:
					ctOut.noise += basis[i].noise * float64(val[0])
				}
			}
		}

		return ctOut, vOut
	} else {
		// Giantstep computation.
		halfDeg := curDeg >> 1

		ctLo, vLo := op.evalRecurse(lo, lo+halfDeg-1, p, basis, basisLazy, rlk)
		ctHi, vHi := op.evalRecurse(lo+halfDeg, hi, p, basis, basisLazy, rlk)

		if ctHi == nil {
			return ctLo, vLo
		}

		// Multiply basis[halfDeg] to ctHi.
		if vHi == nil {
			tarLen, auxIdx, auxMod := op.noise.getAuxMod(ctHi, basis[halfDeg])
			vHi = rlwe.NewVectorCustom(op.params.Rank(), tarLen, 0, 3, true)
			op.tensorTo(vHi, ctHi, basis[halfDeg], tarLen, auxIdx, auxMod, true)
			ctHi.noise = op.noise.tensorTo(ctHi, basis[halfDeg], tarLen, auxIdx, auxMod)
		} else {
			ctHi = ctHi.WithModLen(vHi.BaseModLen())
			op.rlweOp.RelinTo(ctHi.Value, vHi, rlk, true)
			ctHi.noise += op.noise.noise.GadgetProd(ctLen)

			tarLen, auxIdx, auxMod := op.noise.getAuxMod(ctHi, basis[halfDeg])
			vHi = vHi.WithModLen(tarLen, 0)
			op.tensorTo(vHi, ctHi, basis[halfDeg], tarLen, auxIdx, auxMod, true)
			ctHi.noise = op.noise.tensorTo(ctHi, basis[halfDeg], tarLen, auxIdx, auxMod)
		}

		if ctLo == nil {
			return ctHi, vHi
		}

		if vLo == nil {
			tarLen := vHi.BaseModLen()
			ctLoTar := ctLo.WithModLen(tarLen)
			op.ModSwitchTo(ctLoTar, ctLo, tarLen, true)

			pOp.AddTo(vHi.Value[0], vHi.Value[0], ctLoTar.Value.Body)
			pOp.AddTo(vHi.Value[1], vHi.Value[1], ctLoTar.Value.Mask)

			ctHi.noise += ctLoTar.noise
		} else {
			tarLen := vHi.BaseModLen()
			vLoTar := vLo.WithModLen(tarLen, 0)

			pOp.ScaleTo(vLoTar.Value[0], vLo.Value[0], tarLen, true)
			pOp.ScaleTo(vLoTar.Value[1], vLo.Value[1], tarLen, true)
			pOp.ScaleTo(vLoTar.Value[2], vLo.Value[2], tarLen, true)

			for i := tarLen; i < ctLen; i++ {
				switch op.noise.estimType {
				case heint.VarianceType:
					ctLo.noise /= float64(op.params.BaseModulus()[i].Value()) * float64(op.params.BaseModulus()[i].Value())
				case heint.WorstCaseType:
					ctLo.noise /= float64(op.params.BaseModulus()[i].Value())
				}
			}
			if tarLen != ctLen {
				ctLo.noise += op.noise.noise.RoundNoise()
			}

			pOp.AddTo(vHi.Value[0], vHi.Value[0], vLoTar.Value[0])
			pOp.AddTo(vHi.Value[1], vHi.Value[1], vLoTar.Value[1])
			pOp.AddTo(vHi.Value[2], vHi.Value[2], vLoTar.Value[2])

			ctHi.noise += ctLo.noise
		}

		return ctHi, vHi
	}
}
