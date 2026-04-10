package bgv

import (
	"math"

	"github.com/hienaa-org/hienaa/fhe/polyutils.go"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
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

func (op *Operator) computeBasis(p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey) map[int]*Ciphertext {
	basis := make(map[int]*Ciphertext)
	deg := p.Degree()
	maxLevel := int(math.Ceil(num.Log2(deg + 1)))
	maxDeg := (1 << maxLevel) - 1
	babyLevel := (maxLevel + 1) / 2
	babyDeg := 1 << babyLevel

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
		for j := 1; j < (1 << i); j++ {
			// Check if the polynomial is sparse.
			_, check := p.coeffs[1<<i+j]

			// Check if (1<<i+j)-th basis is needed in the babystep generation.
			for idx := i + 1; idx < babyLevel; idx++ {
				_, checkIdx := p.coeffs[1<<idx+1<<i+j]
				check = check || checkIdx
			}

			// Check if (1<<i+j)-th basis is needed in the final evaluation.
			for idx := 0; idx < maxDeg/babyDeg-1; idx++ {
				_, checkIdx := p.coeffs[idx*babyDeg+j]
				check = check || checkIdx
			}

			// If the monomial is needed, compute the monomial.
			if check {
				basis[1<<i+j] = op.ctPool.Get()
				op.MulTo(basis[1<<i+j], basis[1<<i], basis[j], rlk, true)
			}
		}
	}

	return basis
}

// TODO: Change to an in-place algorithm.
// TODO: Optimise later using lazy relin BSGS algorithm.
func (op *Operator) evalRecurse(lo, hi int, p *Polynomial, basis map[int]*Ciphertext, rlk *rlwe.RelinKey) *Ciphertext {
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

	// Corner case (Does anybody really want to evaluate a constant polynomial?).
	if curDeg == 1 {
		ctOut := NewCiphertextCustom(op.params.Rank(), ctLen, true)
		if val, check := p.coeffs[lo]; check {
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

				if val, check := p.coeffs[maxDeg-babyDeg+1]; check {
					if ctOut == nil {
						ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
					}

					op.AddPlainTo(ctOut, ctOut, val, true)
				}

				for j := 1; j < min(halfBabyDeg, deg+babyDeg-maxDeg); j++ {
					if val, check := p.coeffs[maxDeg-babyDeg+j+1]; check {
						if ctOut == nil {
							ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
						}

						op.MulPlainTo(tmpCt, basis[j], val, true)
						op.AddTo(ctOut, ctOut, tmpCt, true)
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
			if val, check := p.coeffs[lo]; check {
				if ctOut == nil {
					ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
				}

				op.AddPlainTo(ctOut, ctOut, val, true)
			}
			for i := 1; i <= min(babyDeg, deg-lo+1); i++ {
				if val, check := p.coeffs[lo+i]; check {
					if ctOut == nil {
						ctOut = NewCiphertextCustom(op.params.Rank(), ctLen, true)
					}

					op.MulPlainTo(tmpCt, basis[i], val, true)
					op.AddTo(ctOut, ctOut, tmpCt, true)
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
