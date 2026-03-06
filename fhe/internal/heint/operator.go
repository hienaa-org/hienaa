package heint

import (
	"cmp"
	"math"
	"slices"
	"sync"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// Operator evaluates operations over [*Ciphertext] and [*Plaintext].
type Operator struct {
	params rlwe.Parameters
	msgMod *num.Modulus

	rlweOp  *rlwe.Operator
	encoder *Encoder
	scFacs  []*rlwe.Element

	sPool  *sync.Pool
	ptPool *sync.Pool
	ctPool *sync.Pool
}

// NewOperator creates a new [Operator].
func NewOperator(params rlwe.Parameters, msgMod *num.Modulus) *Operator {
	rlweOp := rlwe.NewOperator(params)

	ringParams := params.RingParams()
	extraBits := float64(0)
	for _, mod := range params.BaseModulus() {
		extraBits += num.Log2(mod.Value())
	}
	extraBits = extraBits + num.Log2(ringParams.ExpFactor())
	extraLen := int(math.Ceil(extraBits / num.MaxModulusBits))
	baseMod := params.BaseModulus()

	var extraMod []*num.Modulus
	var bitlen float64
	for {
		extraMod = make([]*num.Modulus, extraLen)

		gap, err := dft.NTTPrimeGap(ringParams)
		if err != nil {
			panic(err)
		}

		bitlen = extraBits
		start := (uint64(math.Floor(num.MaxModulus/float64(gap))))*gap + 1
		if start > num.MaxModulus {
			for start > num.MaxModulus {
				start -= gap
			}
		}
		prime := num.MustPrevPrime(start, gap)

		cnt := 0
		for cnt < extraLen {
			primemod := num.NewModulus(prime)
			for {
				_, ok := slices.BinarySearchFunc(baseMod, primemod, func(a, b *num.Modulus) int {
					return cmp.Compare(a.Value(), b.Value())
				})

				if !ok {
					break
				} else {
					prime = num.MustPrevPrime(prime, gap)
					primemod = num.NewModulus(prime)
				}
			}
			extraMod[cnt] = primemod
			bitlen -= num.Log2(primemod.Value())
			prime = num.MustPrevPrime(prime, gap)
			cnt++
		}

		if bitlen > 0 {
			extraLen++
		} else {
			break
		}
	}

	ambBaseMod := append(baseMod, extraMod...)
	encoder := NewEncoder(params, msgMod)

	return &Operator{
		params: params,
		msgMod: msgMod,

		rlweOp:  rlweOp,
		encoder: encoder,
		scFacs:  computeScalingFactor(params.BaseModulus(), msgMod),

		sPool: &sync.Pool{
			New: func() any {
				return rlwe.NewScalar(params, true)
			},
		},
		ptPool: &sync.Pool{
			New: func() any {
				return rlwe.NewElement(params.Rank(), len(ambBaseMod), len(params.AuxModulus()), true)
			},
		},
		ctPool: &sync.Pool{
			New: func() any {
				return rlwe.NewCiphertextCustom(params.Rank(), len(ambBaseMod), 0, true)
			},
		},
	}
}

// Parameters returns the parameters.
func (op *Operator) Parameters() rlwe.Parameters {
	return op.params
}

// Encoder returns the encoder.
func (op *Operator) Encoder() *Encoder {
	return op.encoder
}

// NegTo computes ctOut = -ct.
func (op *Operator) NegTo(ctOut, ct *rlwe.Ciphertext, isNTT bool) {
	op.rlweOp.NegTo(ctOut, ct)
	if isNTT && !ctOut.IsNTT() {
		op.rlweOp.FwdNTTTo(ctOut, ctOut)
	} else if !isNTT && ctOut.IsNTT() {
		op.rlweOp.InvNTTTo(ctOut, ctOut)
	}
}

// AddTo computes ctOut = ct0 + ct1.
func (op *Operator) AddTo(ctOut, ct0, ct1 *rlwe.Ciphertext, isNTT bool) {
	tarLen := min(ct0.BaseModLen(), ct1.BaseModLen())

	c0 := op.ctPool.Get().(*rlwe.Ciphertext)
	c1 := op.ctPool.Get().(*rlwe.Ciphertext)
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)

	c0 = c0.WithModLen(tarLen, 0)
	c1 = c1.WithModLen(tarLen, 0)

	op.rlweOp.ScaleTo(c0, ct0, tarLen, isNTT)
	op.rlweOp.ScaleTo(c1, ct1, tarLen, isNTT)

	op.rlweOp.AddTo(ctOut, c0, c1)
}

// AddPlainTo computes ctOut = ct + pt.
func (op *Operator) AddPlainTo(ctOut, ct *rlwe.Ciphertext, pt *crt.Element, isNTT bool) {
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()

	if pt.Rank() == 1 { // Input is a scalar
		e := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(e)
		e = e.WithModLen(ct.BaseModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.BaseModLen()-1], e)
		rlweOp.AddElementTo(ctOut, ct, e)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut, ctOut)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut, ctOut)
		}
	} else { // Input is a polynomial
		e := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(e)
		e = e.WithModLen(ct.BaseModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.BaseModLen()-1], e)
		if isNTT && !ct.IsNTT() {
			rlweOp.AddElementTo(ctOut, ct, e)
			rlweOp.FwdNTTTo(ctOut, ctOut)
		} else if !isNTT && ct.IsNTT() {
			rlweOp.InvNTTTo(ctOut, ct)
			rlweOp.AddElementTo(ctOut, ctOut, e)
		} else {
			if isNTT {
				pOp.FwdNTTTo(e, e)
			}
			rlweOp.AddElementTo(ctOut, ct, e)
		}
	}
}

// AddElementTo computes ctOut = ct + e.
func (op *Operator) AddElementTo(ctOut, ct *rlwe.Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.BaseModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if e.Rank() == 1 { // Input is a scalar
		eEcd := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.BaseModLen()-1], e)
		rlweOp.AddElementTo(ctOut, ct, eEcd)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut, ctOut)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut, ctOut)
		}
	} else { // Input is a polynomial
		eEcd := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.BaseModLen()-1], e)

		if ct.IsNTT() == e.IsNTT() { // First add, then perform NTT/iNTT needed.
			rlweOp.AddElementTo(ctOut, ct, eEcd)
			if ct.IsNTT() && !isNTT {
				rlweOp.InvNTTTo(ctOut, ctOut)
			} else if !ct.IsNTT() && isNTT {
				rlweOp.FwdNTTTo(ctOut, ctOut)
			}
		} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then add.
			if e.IsNTT() {
				pOp.InvNTTTo(eEcd, eEcd)
			} else {
				pOp.FwdNTTTo(eEcd, eEcd)
			}
			rlweOp.AddElementTo(ctOut, ct, eEcd)
		} else { // First perform NTT/iNTT on ct, then add.
			if ct.IsNTT() {
				rlweOp.InvNTTTo(ctOut, ct)
			} else {
				rlweOp.FwdNTTTo(ctOut, ct)
			}
			rlweOp.AddElementTo(ctOut, ctOut, eEcd)
		}
	}
}

// SubTo computes ctOut = ct0 - ct1.
func (op *Operator) SubTo(ctOut, ct0, ct1 *rlwe.Ciphertext, isNTT bool) {
	tarLen := min(ct0.BaseModLen(), ct1.BaseModLen())

	c0 := op.ctPool.Get().(*rlwe.Ciphertext)
	c1 := op.ctPool.Get().(*rlwe.Ciphertext)
	defer op.ctPool.Put(c0)
	defer op.ctPool.Put(c1)

	c0 = c0.WithModLen(tarLen, 0)
	c1 = c1.WithModLen(tarLen, 0)

	op.rlweOp.ScaleTo(c0, ct0, tarLen, isNTT)
	op.rlweOp.ScaleTo(c1, ct1, tarLen, isNTT)

	op.rlweOp.SubTo(ctOut, c0, c1)
}

// SubPlainTo computes ctOut = ct - pt.
func (op *Operator) SubPlainTo(ctOut, ct *rlwe.Ciphertext, pt *crt.Element, isNTT bool) {
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if pt.Rank() == 1 { // Input is a scalar
		e := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(e)
		e = e.WithModLen(ct.BaseModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.BaseModLen()-1], e)
		rlweOp.SubElementTo(ctOut, ct, e)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut, ctOut)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut, ctOut)
		}
	} else { // Input is a polynomial
		e := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(e)
		e = e.WithModLen(ct.BaseModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		pOp.MulTo(e, op.scFacs[ct.BaseModLen()-1], e)
		if isNTT && !ct.IsNTT() {
			rlweOp.SubElementTo(ctOut, ct, e)
			rlweOp.FwdNTTTo(ctOut, ctOut)
		} else if !isNTT && ct.IsNTT() {
			rlweOp.InvNTTTo(ctOut, ct)
			rlweOp.SubElementTo(ctOut, ctOut, e)
		} else {
			if isNTT {
				pOp.FwdNTTTo(e, e)
			}
			rlweOp.SubElementTo(ctOut, ct, e)
		}
	}
}

// SubElementTo computes ctOut = ct - e.
func (op *Operator) SubElementTo(ctOut, ct *rlwe.Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.BaseModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if e.Rank() == 1 { // Input is a scalar
		eEcd := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.BaseModLen()-1], e)
		rlweOp.SubElementTo(ctOut, ct, eEcd)

		if isNTT && !ctOut.IsNTT() {
			rlweOp.FwdNTTTo(ctOut, ctOut)
		} else if !isNTT && ctOut.IsNTT() {
			rlweOp.InvNTTTo(ctOut, ctOut)
		}
	} else { // Input is a polynomial
		eEcd := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(eEcd)
		eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

		pOp.MulTo(eEcd, op.scFacs[ct.BaseModLen()-1], e)

		if ct.IsNTT() == e.IsNTT() { // First sub, then perform NTT/iNTT needed.
			rlweOp.SubElementTo(ctOut, ct, eEcd)
			if ct.IsNTT() && !isNTT {
				rlweOp.InvNTTTo(ctOut, ctOut)
			} else if !ct.IsNTT() && isNTT {
				rlweOp.FwdNTTTo(ctOut, ctOut)
			}
		} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then sub.
			if e.IsNTT() {
				pOp.InvNTTTo(eEcd, eEcd)
			} else {
				pOp.FwdNTTTo(eEcd, eEcd)
			}
			rlweOp.SubElementTo(ctOut, ct, eEcd)
		} else { // First perform NTT/iNTT on ct, then sub.
			if ct.IsNTT() {
				rlweOp.InvNTTTo(ctOut, ct)
			} else {
				rlweOp.FwdNTTTo(ctOut, ct)
			}
			rlweOp.SubElementTo(ctOut, ctOut, eEcd)
		}
	}
}

// MulPlainTo computes ctOut = ct * pt.
func (op *Operator) MulPlainTo(ctOut, ct *rlwe.Ciphertext, pt *crt.Element, isNTT bool) {
	rlweOp := op.rlweOp
	if pt.Rank() == 1 { // Input is a scalar
		e := op.sPool.Get().(*rlwe.Element)
		defer op.sPool.Put(e)
		e = e.WithModLen(ct.BaseModLen(), 0)

		op.encoder.EncodeTo(e, pt, false)
		rlweOp.MulElementTo(ctOut, ct, e)
		if ct.IsNTT() && !isNTT {
			rlweOp.InvNTTTo(ctOut, ctOut)
		} else if !ct.IsNTT() && isNTT {
			rlweOp.FwdNTTTo(ctOut, ctOut)
		}
	} else { // Input is a polynomial
		e := op.ptPool.Get().(*rlwe.Element)
		defer op.ptPool.Put(e)
		e = e.WithModLen(ct.BaseModLen(), 0)

		op.encoder.EncodeTo(e, pt, true)
		if ct.IsNTT() {
			op.rlweOp.MulElementTo(ctOut, ct, e)
		} else {
			op.rlweOp.FwdNTTTo(ctOut, ct)
			op.rlweOp.MulElementTo(ctOut, ctOut, e)
		}

		if !isNTT {
			op.rlweOp.InvNTTTo(ctOut, ctOut)
		}
	}
}

// MulElementTo computes ctOut = ct * e.
func (op *Operator) MulElementTo(ctOut, ct *rlwe.Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.BaseModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()
	if e.Rank() == 1 || (ct.IsNTT() && e.IsNTT()) {
		rlweOp.MulElementTo(ctOut, ct, e)
	} else {
		if ct.IsNTT() {
			eNTT := op.ptPool.Get().(*rlwe.Element)
			defer op.ptPool.Put(eNTT)
			eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

			pOp.FwdNTTTo(eNTT, e)
			rlweOp.MulElementTo(ctOut, ct, eNTT)
		} else {
			rlweOp.FwdNTTTo(ctOut, ct)

			if e.IsNTT() {
				rlweOp.MulElementTo(ctOut, ctOut, e)
			} else {
				eNTT := op.ptPool.Get().(*rlwe.Element)
				defer op.ptPool.Put(eNTT)
				eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

				pOp.FwdNTTTo(eNTT, e)
				rlweOp.MulElementTo(ctOut, ctOut, eNTT)
			}
		}
	}

	if ctOut.IsNTT() && !isNTT {
		rlweOp.InvNTTTo(ctOut, ctOut)
	} else if !ctOut.IsNTT() && isNTT {
		rlweOp.FwdNTTTo(ctOut, ctOut)
	}
}
