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

		op.MulTo(basis[1<<i], basis[1<<(i-1)], basis[1<<(i-1)], rlk, true)
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

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			op.MulTo(basis[i], basis[idx1], basis[idx2], rlk, true)

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

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			tarLen, auxIdx, auxMod := op.noise.getAuxMod(basis[idx1], basis[idx2])
			basisLazy[i] = basisLazy[i].WithModLen(tarLen, 0)
			basis[i] = &Ciphertext{
				Value: &rlwe.Ciphertext{
					Body: basisLazy[i].Value[0],
					Mask: basisLazy[i].Value[1],
				},
			}

			op.tensorTo(basisLazy[i], basis[idx1], basis[idx2], tarLen, auxIdx, auxMod, true)

			basis[i].noise = op.noise.tensorTo(basis[idx1], basis[idx2], tarLen, auxIdx, auxMod)
		}
	}

	// Evaluate the polynomial.
	ctOut.Resize(ctLen)
	vOut := op.vPool.Get()
	defer op.vPool.Put(vOut)
	vOut = vOut.WithModLen(ctLen, 0)

	isTensor := op.evalRecurseTo(ctOut, vOut, 0, maxDeg, p, basis, basisLazy, rlk)

	if isTensor {
		ctOut.Resize(vOut.BaseModLen())
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

	// Define required instances.
	ctBuf := op.ctPool.Get()
	defer op.ctPool.Put(ctBuf)
	ctBuf = ctBuf.WithModLen(ctLen)

	pBuf := op.ePool.Get(crt.TypePoly)
	defer op.ePool.Put(pBuf)
	pBuf = pBuf.WithModLen(ctLen, 0)

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

					pOp.ScaleTo(pBuf, basisLazy[i].Value[0], ctLen, true)
					pOp.MulAddTo(vOut.Value[0], pBuf, scBuf)
					pOp.ScaleTo(pBuf, basisLazy[i].Value[1], ctLen, true)
					pOp.MulAddTo(vOut.Value[1], pBuf, scBuf)
					pOp.ScaleTo(pBuf, basisLazy[i].Value[2], ctLen, true)
					pOp.MulAddTo(vOut.Value[2], pBuf, scBuf)

					op.noise.ModSwitchTo(ctBuf, basis[i], ctLen)
				} else if !isTensor {
					op.ModSwitchTo(ctBuf, basis[i], ctLen, true)
					pOp.MulAddTo(ctOut.Value.Body, ctBuf.Value.Body, scBuf)
					pOp.MulAddTo(ctOut.Value.Mask, ctBuf.Value.Mask, scBuf)
				} else {
					op.ModSwitchTo(ctBuf, basis[i], ctLen, true)
					pOp.MulAddTo(vOut.Value[0], ctBuf.Value.Body, scBuf)
					pOp.MulAddTo(vOut.Value[1], ctBuf.Value.Mask, scBuf)
				}

				switch op.noise.estimType {
				case heint.VarianceType:
					ctOut.noise += ctBuf.noise * float64(val[0]) * float64(val[0])
				case heint.WorstCaseType:
					ctOut.noise += ctBuf.noise * float64(val[0])
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
			vHiTar := vHi.WithModLen(tarLen, 0)
			ctHiTar := ctHi.WithModLen(tarLen)

			op.tensorTo(vHiTar, ctHi, basis[halfDeg], tarLen, auxIdx, auxMod, true)
			pOp.ScaleTo(vHi.Value[0], vHiTar.Value[0], ctLen, true)
			pOp.ScaleTo(vHi.Value[1], vHiTar.Value[1], ctLen, true)
			pOp.ScaleTo(vHi.Value[2], vHiTar.Value[2], ctLen, true)

			ctHiTar.noise = op.noise.tensorTo(ctHi, basis[halfDeg], tarLen, auxIdx, auxMod)
			op.noise.ModSwitchTo(ctHi, ctHiTar, ctLen)

			if !isTensorLo {
				// Add the output ciphertext to the result.
				pOp.AddTo(vOut.Value[0], vHi.Value[0], ctOut.Value.Body)
				pOp.AddTo(vOut.Value[1], vHi.Value[1], ctOut.Value.Mask)
				vOut.Value[2].CopyFrom(vHi.Value[2])

				ctOut.noise = ctOut.noise + ctHi.noise
			} else {
				pOp.AddTo(vOut.Value[0], vHi.Value[0], vOut.Value[0])
				pOp.AddTo(vOut.Value[1], vHi.Value[1], vOut.Value[1])
				pOp.AddTo(vOut.Value[2], vHi.Value[2], vOut.Value[2])

				ctOut.noise += ctHi.noise
			}

			return true
		}
	}
}
