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
	switch pType {
	case polyutils.Monomial:
		coeffs := make(map[int]Plaintext)
		for i, c := range c {
			if num.Reduce(c, msgMod) != 0 {
				coeffs[i] = NewScalarFrom(c, msgMod)
			}
		}

		if len(coeffs) == 0 {
			coeffs[0] = NewScalarFrom(0, msgMod)
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

// EvaluatePoly evaluates the polynomial at the given ciphertext.
func (op *Operator) EvaluatePoly(p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) *Ciphertext {
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.ModLen(), true)
	op.EvaluatePolyTo(ctOut, p, ct, rlk, isNTT)
	return ctOut
}

// EvaluatePoly evaluates the polynomial at the given ciphertext.
func (op *Operator) EvaluatePolyTo(ctOut *Ciphertext, p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) {
	// First compute the basis.
	basis := make(map[int]*Ciphertext)
	basisLazy := make(map[int]*rlwe.Vector)

	// Hyperparameters
	ctLen := ct.ModLen()
	deg := p.Degree()
	maxLevel := polyutils.MaxLevel(deg)
	maxDeg := polyutils.MaxDeg(deg)
	babyLevel := polyutils.BabyLevel(deg)
	babyDeg := polyutils.BabyDeg(deg)

	ctOut.Resize(ctLen)
	// Handle the edge case.
	if p.Degree() == 0 {
		ctOut.Clear()
		op.AddPlainTo(ctOut, ctOut, p.coeffs[0], isNTT)
		return
	}

	// Compute the power-of-two monomials.
	basis[1] = op.ctPool.Get()
	defer op.ctPool.Put(basis[1])
	basis[1] = basis[1].WithModLen(ctLen)

	if ct.IsNTT() {
		basis[1].CopyFrom(ct)
	} else {
		op.FwdNTTTo(basis[1], ct)
	}

	for i := 1; i < maxLevel; i++ {
		basis[1<<i] = op.ctPool.Get()
		defer op.ctPool.Put(basis[1<<i])
		basis[1<<i] = basis[1<<i].WithModLen(ctLen)

		if i < babyLevel {
			bsTmp := basis[1<<i].WithModLen(ctLen)
			op.MulTo(bsTmp, basis[1<<(i-1)], basis[1<<(i-1)], rlk, true)
			op.ModSwitchTo(basis[1<<i], bsTmp, ctLen, true)
		} else {
			op.MulTo(basis[1<<i], basis[1<<(i-1)], basis[1<<(i-1)], rlk, true)
		}
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
			bsTmp := basis[i].WithModLen(ctLen)

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			op.MulTo(bsTmp, basis[idx1], basis[idx2], rlk, true)
			op.ModSwitchTo(basis[i], bsTmp, ctLen, true)

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

			basis[i] = &Ciphertext{
				Value: &rlwe.Ciphertext{
					Body: basisLazy[i].Value[0],
					Mask: basisLazy[i].Value[1],
				},
			}

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			tarLen, auxIdx, auxMod := op.noise.getAuxMod(basis[idx1], basis[idx2])

			lazyTmp := basisLazy[i].WithModLen(tarLen, 0)
			op.tensorTo(lazyTmp, basis[idx1], basis[idx2], tarLen, auxIdx, auxMod, true)

			pOp := op.rlweOp.PlainOperator()
			pOp.ScaleTo(basisLazy[i].Value[0], lazyTmp.Value[0], ctLen, true)
			pOp.ScaleTo(basisLazy[i].Value[1], lazyTmp.Value[1], ctLen, true)
			pOp.ScaleTo(basisLazy[i].Value[2], lazyTmp.Value[2], ctLen, true)

			basis[i].noise = op.noise.noise.ModSwitch(op.noise.tensor(basis[idx1], basis[idx2], tarLen, auxIdx, auxMod), tarLen, ctLen)
		}
	}

	// Evaluate the polynomial.
	ctOut.Resize(ctLen)
	vOut := op.vPool.Get()
	defer op.vPool.Put(vOut)
	vOut = vOut.WithModLen(ctLen, 0)

	isTensor := op.evalRecurseTo(ctOut, vOut, 0, maxDeg, p, basis, basisLazy, rlk)

	if isTensor {
		op.rlweOp.RelinTo(ctOut.Value, vOut, rlk, true)
		ctOut.noise += op.noise.noise.GadgetProd(ctLen)
	}

	op.RescaleTo(ctOut, ctOut, isNTT)
}

// evalRecurseTo evaluates the polynomial at the given range using the Paterson-Stockmeyer algorithm.
func (op *Operator) evalRecurseTo(ctOut *Ciphertext, vOut *rlwe.Vector, lo, hi int, p *Polynomial, basis map[int]*Ciphertext, basisLazy map[int]*rlwe.Vector, rlk *rlwe.RelinKey) bool {
	// Hyperparameters.
	deg := p.Degree()
	maxDeg := polyutils.MaxDeg(deg)
	babyDeg := polyutils.BabyDeg(deg)
	ctLen := basis[1].ModLen()

	// Current degree.
	curDeg := hi - lo + 1
	if !num.IsPowerOfTwo(curDeg) {
		panic("Current degree is not a power of two")
	}

	// Initialize the output.
	ctOut.Clear()
	vOut.Clear()
	isTensor := false

	// Paterson-Stockmeyer algorithm.
	pOp := op.rlweOp.PlainOperator()
	if curDeg <= babyDeg && (hi != maxDeg || curDeg == 2) {
		// Babystep computation.
		if lo > deg {
			return false
		}

		scBuf := op.ePool.Get(crt.TypeScalar)
		defer op.ePool.Put(scBuf)
		scBuf = scBuf.WithModLen(ctLen, 0)

		if val, check := p.coeffs[lo]; check {
			op.AddPlainTo(ctOut, ctOut, val, true)
		}

		for i := 1; i <= min(hi, deg)-lo; i++ {
			if val, check := p.coeffs[lo+i]; check {
				op.intOp.Encoder().EncodeTo(scBuf, val, true)

				if basisLazy[i] != nil {
					if !isTensor {
						vOut.Value[0].CopyFrom(ctOut.Value.Body)
						vOut.Value[1].CopyFrom(ctOut.Value.Mask)
						vOut.Value[2].Clear()
						isTensor = true
					}

					pOp.MulAddTo(vOut.Value[0], basisLazy[i].Value[0], scBuf)
					pOp.MulAddTo(vOut.Value[1], basisLazy[i].Value[1], scBuf)
					pOp.MulAddTo(vOut.Value[2], basisLazy[i].Value[2], scBuf)
				} else if !isTensor {
					pOp.MulAddTo(ctOut.Value.Body, basis[i].Value.Body, scBuf)
					pOp.MulAddTo(ctOut.Value.Mask, basis[i].Value.Mask, scBuf)
				} else {
					pOp.MulAddTo(vOut.Value[0], basis[i].Value.Body, scBuf)
					pOp.MulAddTo(vOut.Value[1], basis[i].Value.Mask, scBuf)
				}

				var valNorm float64
				if val[0] < op.msgMod.Value()>>1 {
					valNorm = float64(val[0])
				} else {
					valNorm = float64(op.msgMod.Value() - val[0])
				}
				switch op.noise.estimType {
				case heint.VarianceType:
					ctOut.noise += basis[i].noise * valNorm * valNorm
				case heint.WorstCaseType:
					ctOut.noise += basis[i].noise * valNorm
				}
			}
		}

		return isTensor
	} else {
		// Giantstep computation.
		halfDeg := curDeg >> 1

		isTensorLo := op.evalRecurseTo(ctOut, vOut, lo, lo+halfDeg-1, p, basis, basisLazy, rlk)
		if lo+halfDeg > deg {
			return isTensorLo
		} else {
			ctHi := op.ctPool.Get()
			defer op.ctPool.Put(ctHi)
			ctHi = ctHi.WithModLen(ctLen)

			vHi := op.vPool.Get()
			defer op.vPool.Put(vHi)
			vHi = vHi.WithModLen(ctLen, 0)

			isTensorHi := op.evalRecurseTo(ctHi, vHi, lo+halfDeg, hi, p, basis, basisLazy, rlk)
			if isTensorHi {
				// Relin the output and tensor.
				ctHi = ctHi.WithModLen(vHi.BaseModLen())
				op.rlweOp.RelinTo(ctHi.Value, vHi, rlk, true)
				ctHi.noise += op.noise.noise.GadgetProd(ctLen)
			}
			tarLen, auxIdx, auxMod := op.noise.getAuxMod(ctHi, basis[halfDeg])
			vHi = vHi.WithModLen(tarLen, 0)

			op.tensorTo(vHi, ctHi, basis[halfDeg], tarLen, auxIdx, auxMod, true)
			ctHi.noise = op.noise.tensor(ctHi, basis[halfDeg], tarLen, auxIdx, auxMod)

			if !isTensorLo {
				// Add the output ciphertext to the result.
				op.ModSwitchTo(ctOut, ctOut, tarLen, true)
				vOut.Resize(tarLen, 0)

				pOp.AddTo(vOut.Value[0], vHi.Value[0], ctOut.Value.Body)
				pOp.AddTo(vOut.Value[1], vHi.Value[1], ctOut.Value.Mask)
				vOut.Value[2].CopyFrom(vHi.Value[2])

			} else {
				// Add the output vector to the result.
				vOutTar := vOut.WithModLen(tarLen, 0)
				pOp.ScaleTo(vOutTar.Value[0], vOut.Value[0], tarLen, true)
				pOp.ScaleTo(vOutTar.Value[1], vOut.Value[1], tarLen, true)
				pOp.ScaleTo(vOutTar.Value[2], vOut.Value[2], tarLen, true)
				vOut.Resize(tarLen, 0)

				ctOut.Resize(tarLen)
				ctOut.noise = op.noise.noise.ModSwitch(ctHi.noise, ctLen, tarLen)

				pOp.AddTo(vOut.Value[0], vHi.Value[0], vOut.Value[0])
				pOp.AddTo(vOut.Value[1], vHi.Value[1], vOut.Value[1])
				pOp.AddTo(vOut.Value[2], vHi.Value[2], vOut.Value[2])
			}

			ctOut.noise += ctHi.noise

			return true
		}
	}
}
