package ckks

import (
	"math"

	"github.com/hienaa-org/hienaa/fhe/polyutils.go"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

// Polynomial is a struct for polynomial in CKKS scheme.
type Polynomial struct {
	coeffs map[int]*Plaintext
	// bsgsCoeffs is the coefficients for the Paterson-Stockmeyer algorithm for Chebyshev basis.
	bsgsCoeffs map[int]*Plaintext
	polyType   polyutils.PolynomialType
}

// NewPolynomial creates a new polynomial from the given coefficients and type.
func NewPolynomial(c []float64, pType polyutils.PolynomialType, vType ValueType) *Polynomial {
	// TODO: Allow complex coefficient polynomials.
	switch pType {
	case polyutils.Monomial:
		coeffs := make(map[int]*Plaintext)
		for i, c := range c {
			coeffs[i] = NewScalarFrom(c, vType)
		}

		return &Polynomial{
			coeffs:     coeffs,
			bsgsCoeffs: coeffs,
			polyType:   pType,
		}
	case polyutils.Chebyshev:
		deg := len(c) - 1
		maxLevel := int(math.Ceil(num.Log2(deg + 1)))
		babyLevel := (maxLevel + 1) / 2
		cBSGS := make([]float64, 1<<maxLevel)
		copy(cBSGS[:deg+1], c)

		for i := maxLevel - 1; i > 0; i-- {
			curDeg := 1 << i

			if i >= babyLevel {
				for j := 0; j < 1<<(maxLevel-i-1); j++ {
					hi := 2 * curDeg * (j + 1)
					for k := 0; k < curDeg-1; k++ {
						cBSGS[hi-curDeg-k-1] -= cBSGS[hi-curDeg+k+1]
						cBSGS[hi-curDeg+k+1] *= 2
					}
				}
			} else {
				for j := 0; j < curDeg-1; j++ {
					cBSGS[1<<maxLevel-curDeg-j-1] -= cBSGS[1<<maxLevel-curDeg+j+1]
					cBSGS[1<<maxLevel-curDeg+j+1] *= 2
				}
			}
		}
		cBSGS = cBSGS[:deg+1]

		coeffs := make(map[int]*Plaintext)
		bsgsCoeffs := make(map[int]*Plaintext)
		for i := 0; i <= deg; i++ {
			coeffs[i] = NewScalarFrom(c[i], vType)
			bsgsCoeffs[i] = NewScalarFrom(cBSGS[i], vType)
		}

		return &Polynomial{
			coeffs:     coeffs,
			bsgsCoeffs: bsgsCoeffs,
			polyType:   pType,
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
func (p *Polynomial) Coeffs() map[int]*Plaintext {
	return p.coeffs
}

// TODO: Move to the right place.
// TODO: Change Evaluate and evalRecurse to a in-place algorithm.
// Evaluate evaluates the polynomial at the given ciphertext.
func (op *Operator) Evaluate(p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) *Ciphertext {
	// Compute the basis for babystep.
	basis := op.computeBasis(p, ct, rlk)

	// Evaluate the polynomial.
	deg := p.Degree()
	maxLevel := int(math.Ceil(num.Log2(deg + 1)))
	maxDeg := (1 << maxLevel) - 1
	ctOut := op.evalRecurse(0, maxDeg, p, basis, rlk)

	// Inverse NTT if needed.
	if !isNTT {
		op.InvNTTTo(ctOut, ctOut)
	}

	// clean up the basis.
	for _, ct := range basis {
		op.ctPool.Put(ct)
	}

	return ctOut
}

// computeBasis computes the basis for the polynomial evaluation.
func (op *Operator) computeBasis(p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey) map[int]*Ciphertext {
	basis := make(map[int]*Ciphertext)
	deg := p.Degree()
	maxLevel := int(math.Ceil(num.Log2(deg + 1)))
	maxDeg := (1 << maxLevel) - 1
	babyLevel := (maxLevel + 1) / 2
	babyDeg := 1 << babyLevel

	switch p.polyType {
	case polyutils.Monomial:
		// Compute the power-of-two monomials.
		basis[1] = op.ctPool.Get()
		basis[1] = basis[1].WithModLen(ct.ModLen())
		if ct.IsNTT() {
			basis[1].CopyFrom(ct)
		} else {
			op.FwdNTTTo(basis[1], ct)
		}

		for i := 1; i < maxLevel; i++ {
			basis[1<<i] = op.ctPool.Get()
			op.MulTo(basis[1<<i], basis[1<<(i-1)], basis[1<<(i-1)], rlk, true)
		}

		// Compute other monomials.
		for i := 1; i < babyLevel; i++ {
			for j := 1; j < 1<<i; j++ {
				// Check if the polynomial is sparse.
				val, check := p.coeffs[1<<i+j]
				check = check && math.Abs(val.Value[0]) >= 1/op.scFac

				// Check if (1<<i+j)-th basis is needed in the babystep generation.
				for idx := i + 1; idx < babyLevel; idx++ {
					val, checkIdx := p.coeffs[1<<idx+1<<i+j]
					check = check || (checkIdx && math.Abs(val.Value[0]) >= 1/op.scFac)
				}

				// Check if (1<<i+j)-th basis is needed in the final evaluation.
				for idx := 0; idx < maxDeg/babyDeg-1; idx++ {
					val, checkIdx := p.coeffs[idx*babyDeg+j]
					check = check || (checkIdx && math.Abs(val.Value[0]) >= 1/op.scFac)
				}

				// If the monomial is needed, compute the monomial.
				if check {
					basis[1<<i+j] = op.ctPool.Get()
					op.MulTo(basis[1<<i+j], basis[1<<i], basis[j], rlk, true)
				}
			}
		}
	case polyutils.Chebyshev:
		// Precompute constants.
		one := NewScalarFrom(1, TypeInt)
		two := NewScalarFrom(2, TypeInt)

		// Compute the power-of-two bases.
		// T_2^i(x) = 2 * T_2^{i-1}(x) * T_2^{i-1}(x) - 1
		basis[1] = op.ctPool.Get()
		basis[1] = basis[1].WithModLen(ct.ModLen())
		if ct.IsNTT() {
			basis[1].CopyFrom(ct)
		} else {
			op.FwdNTTTo(basis[1], ct)
		}

		for i := 1; i < maxLevel; i++ {
			basis[1<<i] = op.ctPool.Get()
			op.MulTo(basis[1<<i], basis[1<<(i-1)], basis[1<<(i-1)], rlk, true)
			op.MulPlainTo(basis[1<<i], basis[1<<i], two, true)
			op.SubPlainTo(basis[1<<i], basis[1<<i], one, true)
		}

		// Compute other bases.
		for i := 1; i < babyLevel; i++ {
			for j := 1; j < 1<<i; j++ {
				// Checks if the polynomial is sparse.
				val, check := p.bsgsCoeffs[1<<i+j]
				check = check || math.Abs(val.Value[0]) >= 1/op.scFac

				// Check if (1<<i+j)-th basis is needed in the babystep generation.
				for idx := i + 1; idx < babyLevel; idx++ {
					val, checkIdx := p.bsgsCoeffs[1<<idx+1<<i+j]
					check = check || (checkIdx && math.Abs(val.Value[0]) >= 1/op.scFac)
				}

				// Check if (1<<i+j)-th basis is needed in the final evaluation.
				for idx := 0; idx < maxDeg/babyDeg-1; idx++ {
					val, checkIdx := p.bsgsCoeffs[idx*babyDeg+j]
					check = check || (checkIdx && math.Abs(val.Value[0]) >= 1/op.scFac)
				}

				// If the basis is needed, compute the basis.
				if check {
					basis[1<<i+j] = op.ctPool.Get()
					op.MulTo(basis[1<<i+j], basis[1<<i], basis[j], rlk, true)
					op.MulPlainTo(basis[1<<i+j], basis[1<<i+j], two, true)
					op.SubTo(basis[1<<i+j], basis[1<<i+j], basis[1<<i-j], true)
				}
			}
		}
	default:
		panic("unsupported polynomial type")
	}

	return basis
}

// TODO: Change to an in-place algorithm.
// TODO: Optimise later using lazy relin BSGS algorithm.
func (op *Operator) evalRecurse(lo, hi int, p *Polynomial, basis map[int]*Ciphertext, rlk *rlwe.RelinKey) *Ciphertext {
	deg := p.Degree()
	maxLevel := int(math.Ceil(num.Log2(deg + 1)))
	maxDeg := (1 << maxLevel) - 1
	babyLevel := (maxLevel + 1) / 2
	babyDeg := 1 << babyLevel
	ctLen := basis[1].ModLen()
	scFac := basis[1].scFac

	// Current degree.
	curDeg := hi - lo + 1
	if !num.IsPowerOfTwo(curDeg) {
		panic("Current degree is not a power of two")
	}

	// Corner cases.
	if curDeg == 1 {
		ctOut := NewCiphertextCustom(op.params.Rank(), ctLen, true)
		ctOut.SetScalingFactor(scFac)
		if val, check := p.coeffs[lo]; check && math.Abs(val.Value[0]) >= 1/op.scFac {
			op.AddPlainTo(ctOut, ctOut, val, false)
		}
		return ctOut
	}

	// Temporary variables.
	tmpCt := op.ctPool.Get()
	defer op.ctPool.Put(tmpCt)
	tmpCt = tmpCt.WithModLen(ctLen)

	if curDeg == babyDeg {
		if lo > deg {
			return nil
		}

		var ctOut *Ciphertext
		if hi == maxDeg {
			// Babystep computation, at the highest degree.
			// We need to perform BSGS to the smallest degree to minimise the level consumption.
			for i := 0; i <= babyLevel; i++ {
				babyDeg := 1 << i
				halfBabyDeg := babyDeg >> 1
				if maxDeg-babyDeg >= deg {
					continue
				}

				if val, check := p.bsgsCoeffs[maxDeg-babyDeg+1]; check && math.Abs(val.Value[0]) >= 1/op.scFac {
					if ctOut == nil {
						ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
						ctOut.SetScalingFactor(scFac)
					}

					op.AddPlainTo(ctOut, ctOut, val, false)
				}

				for j := 1; j < min(halfBabyDeg, deg+babyDeg-maxDeg); j++ {
					if val, check := p.bsgsCoeffs[maxDeg-babyDeg+j+1]; check && math.Abs(val.Value[0]) >= 1/op.scFac {
						if ctOut == nil {
							ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
							ctOut.SetScalingFactor(scFac)
						}

						op.MulPlainTo(tmpCt, basis[j], val, false)
						op.AddTo(ctOut, ctOut, tmpCt, false)
					}
				}

				if i < babyLevel {
					if ctOut != nil {
						op.MulTo(ctOut, ctOut, basis[1<<i], rlk, true)
					}
				}
			}
		} else {
			// Babystep computation.
			if val, check := p.bsgsCoeffs[lo]; check && math.Abs(val.Value[0]) >= 1/op.scFac {
				if ctOut == nil {
					ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
					ctOut.SetScalingFactor(scFac)
				}

				op.AddPlainTo(ctOut, ctOut, val, false)
			}
			for i := 1; i <= min(curDeg-1, deg-lo); i++ {
				if val, check := p.bsgsCoeffs[lo+i]; check && math.Abs(val.Value[0]) >= 1/op.scFac {
					if ctOut == nil {
						ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
						ctOut.SetScalingFactor(scFac)
					}

					op.MulPlainTo(tmpCt, basis[i], val, false)
					op.AddTo(ctOut, ctOut, tmpCt, false)
				}
			}
		}

		return ctOut
	} else {
		// Giantstep computation.
		halfDeg := curDeg >> 1

		lo_res := op.evalRecurse(lo, lo+halfDeg-1, p, basis, rlk)
		hi_res := op.evalRecurse(lo+halfDeg, hi, p, basis, rlk)

		if hi_res == nil {
			return lo_res
		} else {
			if lo_res == nil {
				op.MulTo(hi_res, hi_res, basis[halfDeg], rlk, true)
				return hi_res
			} else {
				op.MulTo(tmpCt, hi_res, basis[halfDeg], rlk, true)
				op.AddTo(lo_res, lo_res, tmpCt, true)
				return lo_res
			}
		}
	}
}
