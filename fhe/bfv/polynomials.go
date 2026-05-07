package bfv

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
	ctOut := NewCiphertextCustom(op.params.Rank(), ct.ModLen(), isNTT)
	op.EvaluatePolyTo(ctOut, p, ct, rlk, isNTT)
	return ctOut
}

// Evaluate evaluates the polynomial at the given ciphertext.
func (op *Operator) EvaluatePolyTo(ctOut *Ciphertext, p *Polynomial, ct *Ciphertext, rlk *rlwe.RelinKey, isNTT bool) {
	// Handle the edge case.
	if p.Degree() == 0 {
		ctOut.Clear()
		op.AddPlainTo(ctOut, ctOut, p.coeffs[0], isNTT)
	}

	// First compute the basis.
	basis := make(map[int]*Ciphertext)
	basisLazy := make(map[int]*rlwe.Vector)

	// Hyperparameters
	ctLen := ct.ModLen()
	ambLen := op.noise.getAmbLen(ct, ct)
	deg := p.Degree()
	maxLevel := polyutils.MaxLevel(deg)
	maxDeg := polyutils.MaxDeg(deg)
	babyLevel := polyutils.BabyLevel(deg)
	babyDeg := polyutils.BabyDeg(deg)

	// BABYSTEP COMPUTATION
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
			basisLazy[i] = basisLazy[i].WithModLen(ctLen+ambLen, 0)

			basis[i] = &Ciphertext{
				Value: &rlwe.Ciphertext{
					Body: basisLazy[i].Value[0].WithModLen(ctLen, 0),
					Mask: basisLazy[i].Value[1].WithModLen(ctLen, 0),
				},
			}

			idx1 := 1 << int(math.Floor(num.Log2(i)))
			idx2 := i - idx1

			op.liftAndTensorTo(basisLazy[i], basis[idx1], basis[idx2], ambLen, true)
			basis[i].noise = op.noise.liftAndTensor(basis[idx1], basis[idx2], ambLen)
		}
	}

	// Evaluate the polynomial.
	ctOut.Resize(ctLen)
	vOut := op.vPool.Get()
	defer op.vPool.Put(vOut)
	vOut = vOut.WithModLen(ctLen+ambLen, 0)

	isTensor := op.evalRecurseTo(ctOut, vOut, 0, maxDeg, p, basis, basisLazy, rlk)

	if isTensor {
		vOutBase := vOut.WithModLen(ctLen, 0)
		op.divRoundTo(vOutBase, vOut, true)
		op.rlweOp.RelinTo(ctOut.Value, vOutBase, rlk, true)
		ctOut.noise += op.noise.noise.GadgetProd(ctLen)
	}
}

// evalRecurseTo evaluates the polynomial at the given range using the Paterson-Stockmeyer algorithm.
func (op *Operator) evalRecurseTo(ctOut *Ciphertext, vOut *rlwe.Vector, lo, hi int, p *Polynomial, basis map[int]*Ciphertext, basisLazy map[int]*rlwe.Vector, rlk *rlwe.RelinKey) bool {
	// Hyperparameters.
	deg := p.Degree()
	maxDeg := polyutils.MaxDeg(deg)
	babyDeg := polyutils.BabyDeg(deg)
	ctLen := basis[1].ModLen()
	ambLen := op.noise.getAmbLen(basis[1], basis[1])

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
	ctBuf = ctBuf.WithModLen(ctLen + ambLen)

	// Paterson-Stockmeyer algorithm.
	pOp := op.ambOp.PlainOperator()
	if curDeg <= babyDeg && (hi != maxDeg || curDeg == 2) {
		// Babystep computation.
		if lo > deg {
			return false
		}

		outMod := op.ambOp.Params.FullModulus()[:ctLen+ambLen]
		emb := crt.NewVecEmbedder(outMod, []*num.Modulus{op.msgMod})

		scBuf := op.ePool.Get(crt.TypeScalar)
		defer op.ePool.Put(scBuf)
		scBufAmb := scBuf.WithModLen(ctLen+ambLen, 0)
		scBuf = scBuf.WithModLen(ctLen, 0)

		if val, check := p.coeffs[lo]; check {
			op.AddPlainTo(ctOut, ctOut, val, true)
		}

		for i := 1; i <= min(hi, deg)-lo; i++ {
			if val, check := p.coeffs[lo+i]; check {
				scBufAmb.Value.Coeffs[0][0] = val[0]
				emb.EmbedTo(scBufAmb.Value.Coeffs, scBufAmb.Value.Coeffs[:1])

				if basisLazy[i] != nil {
					if !isTensor {
						pOp.ScaleTo(vOut.Value[0], ctOut.Value.Body, ctLen+ambLen, true)
						pOp.ScaleTo(vOut.Value[1], ctOut.Value.Mask, ctLen+ambLen, true)
						vOut.Value[2].Clear()
						isTensor = true
					}

					pOp.MulAddTo(vOut.Value[0], basisLazy[i].Value[0], scBufAmb)
					pOp.MulAddTo(vOut.Value[1], basisLazy[i].Value[1], scBufAmb)
					pOp.MulAddTo(vOut.Value[2], basisLazy[i].Value[2], scBufAmb)
				} else if !isTensor {
					pOp.MulAddTo(ctOut.Value.Body, basis[i].Value.Body, scBuf)
					pOp.MulAddTo(ctOut.Value.Mask, basis[i].Value.Mask, scBuf)
				} else {
					op.ambOp.ScaleTo(ctBuf.Value, basis[i].Value, ctLen+ambLen, true)
					pOp.MulAddTo(vOut.Value[0], ctBuf.Value.Body, scBufAmb)
					pOp.MulAddTo(vOut.Value[1], ctBuf.Value.Mask, scBufAmb)
				}

				switch op.noise.estimType {
				case heint.VarianceType:
					ctOut.noise += basis[i].noise * float64(val[0]) * float64(val[0])
				case heint.WorstCaseType:
					ctOut.noise += basis[i].noise * float64(val[0])
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
			vHi = vHi.WithModLen(ctLen+ambLen, 0)

			isTensorHi := op.evalRecurseTo(ctHi, vHi, lo+halfDeg, hi, p, basis, basisLazy, rlk)
			if isTensorHi {
				// Relin the output and tensor.
				vHiDiv := vHi.WithModLen(ctLen, 0)
				op.divRoundTo(vHiDiv, vHi, true)
				op.rlweOp.RelinTo(ctHi.Value, vHiDiv, rlk, true)

				ctHi.noise += op.noise.noise.GadgetProd(ctLen) + op.noise.noise.RoundNoise()
			}
			op.liftAndTensorTo(vHi, ctHi, basis[halfDeg], ambLen, true)
			ctHi.noise = op.noise.liftAndTensor(ctHi, basis[halfDeg], ambLen)

			if !isTensorLo {
				// Add the output ciphertext to the result.
				op.ambOp.ScaleTo(ctBuf.Value, ctOut.Value, ctLen+ambLen, true)
				pOp.AddTo(vOut.Value[0], vHi.Value[0], ctBuf.Value.Body)
				pOp.AddTo(vOut.Value[1], vHi.Value[1], ctBuf.Value.Mask)
				vOut.Value[2].CopyFrom(vHi.Value[2])

				ctOut.noise += ctHi.noise + op.noise.noise.RoundNoise()
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
