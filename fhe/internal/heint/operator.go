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

	ePool  *rlwe.ElementPool
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

		ePool: rlwe.NewElementPool(params, true, true),
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
func (op *Operator) AddPlainTo(ctOut, ct *rlwe.Ciphertext, pt []uint64, isNTT bool) {
	rlweOp := op.rlweOp
	pOp := rlweOp.PlainOperator()

	var e *rlwe.Element
	switch len(pt) {
	case 1:
		e = op.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		e = op.ePool.Get(crt.TypePoly)
	}
	defer op.ePool.Put(e)
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

// AddElementTo computes ctOut = ct + e.
func (op *Operator) AddElementTo(ctOut, ct *rlwe.Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.BaseModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.rlweOp.PlainOperator()

	eEcd := op.ePool.Get(e.Type())
	defer op.ePool.Put(eEcd)
	eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

	pOp.MulTo(eEcd, op.scFacs[ct.BaseModLen()-1], e)

	if ct.IsNTT() == e.IsNTT() { // First add, then perform NTT/iNTT needed.
		op.rlweOp.AddElementTo(ctOut, ct, eEcd)
		if ct.IsNTT() && !isNTT {
			op.rlweOp.InvNTTTo(ctOut, ctOut)
		} else if !ct.IsNTT() && isNTT {
			op.rlweOp.FwdNTTTo(ctOut, ctOut)
		}
	} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then add.
		if e.IsNTT() {
			pOp.InvNTTTo(eEcd, eEcd)
		} else {
			pOp.FwdNTTTo(eEcd, eEcd)
		}
		op.rlweOp.AddElementTo(ctOut, ct, eEcd)
	} else { // First perform NTT/iNTT on ct, then add.
		if ct.IsNTT() {
			op.rlweOp.InvNTTTo(ctOut, ct)
		} else {
			op.rlweOp.FwdNTTTo(ctOut, ct)
		}
		op.rlweOp.AddElementTo(ctOut, ctOut, eEcd)
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
func (op *Operator) SubPlainTo(ctOut, ct *rlwe.Ciphertext, pt []uint64, isNTT bool) {
	pOp := op.rlweOp.PlainOperator()

	var e *rlwe.Element
	switch len(pt) {
	case 1:
		e = op.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		e = op.ePool.Get(crt.TypePoly)
	}
	defer op.ePool.Put(e)
	e = e.WithModLen(ct.BaseModLen(), 0)

	op.encoder.EncodeTo(e, pt, false)
	pOp.MulTo(e, op.scFacs[ct.BaseModLen()-1], e)
	if isNTT && !ct.IsNTT() {
		op.rlweOp.SubElementTo(ctOut, ct, e)
		op.rlweOp.FwdNTTTo(ctOut, ctOut)
	} else if !isNTT && ct.IsNTT() {
		op.rlweOp.InvNTTTo(ctOut, ct)
		op.rlweOp.SubElementTo(ctOut, ctOut, e)
	} else {
		if isNTT {
			pOp.FwdNTTTo(e, e)
		}
		op.rlweOp.SubElementTo(ctOut, ct, e)
	}
}

// SubElementTo computes ctOut = ct - e.
func (op *Operator) SubElementTo(ctOut, ct *rlwe.Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.BaseModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.rlweOp.PlainOperator()

	eEcd := op.ePool.Get(e.Type())
	defer op.ePool.Put(eEcd)
	eEcd = eEcd.WithModLen(e.BaseModLen(), 0)

	pOp.MulTo(eEcd, op.scFacs[ct.BaseModLen()-1], e)

	if ct.IsNTT() == e.IsNTT() { // First sub, then perform NTT/iNTT needed.
		op.rlweOp.SubElementTo(ctOut, ct, eEcd)
		if ct.IsNTT() && !isNTT {
			op.rlweOp.InvNTTTo(ctOut, ctOut)
		} else if !ct.IsNTT() && isNTT {
			op.rlweOp.FwdNTTTo(ctOut, ctOut)
		}
	} else if e.IsNTT() != isNTT { // First perform NTT/iNTT on e, then sub.
		if e.IsNTT() {
			pOp.InvNTTTo(eEcd, eEcd)
		} else {
			pOp.FwdNTTTo(eEcd, eEcd)
		}
		op.rlweOp.SubElementTo(ctOut, ct, eEcd)
	} else { // First perform NTT/iNTT on ct, then sub.
		if ct.IsNTT() {
			op.rlweOp.InvNTTTo(ctOut, ct)
		} else {
			op.rlweOp.FwdNTTTo(ctOut, ct)
		}
		op.rlweOp.SubElementTo(ctOut, ctOut, eEcd)
	}
}

// MulPlainTo computes ctOut = ct * pt.
func (op *Operator) MulPlainTo(ctOut, ct *rlwe.Ciphertext, pt []uint64, isNTT bool) {
	var e *rlwe.Element
	switch len(pt) {
	case 1:
		e = op.ePool.Get(crt.TypeScalar)
	case ct.Rank():
		e = op.ePool.Get(crt.TypePoly)
	}
	defer op.ePool.Put(e)
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

// MulElementTo computes ctOut = ct * e.
func (op *Operator) MulElementTo(ctOut, ct *rlwe.Ciphertext, e *rlwe.Element, isNTT bool) {
	if ct.BaseModLen() != e.BaseModLen() || e.AuxModLen() > 0 {
		panic("inconsistent input(s)")
	}

	pOp := op.rlweOp.PlainOperator()
	if e.Type() == crt.TypeScalar || (ct.IsNTT() && e.IsNTT()) {
		op.rlweOp.MulElementTo(ctOut, ct, e)
	} else {
		if ct.IsNTT() {
			eNTT := op.ePool.Get(e.Type())
			defer op.ePool.Put(eNTT)
			eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

			pOp.FwdNTTTo(eNTT, e)
			op.rlweOp.MulElementTo(ctOut, ct, eNTT)
		} else {
			op.rlweOp.FwdNTTTo(ctOut, ct)

			if e.IsNTT() {
				op.rlweOp.MulElementTo(ctOut, ctOut, e)
			} else {
				eNTT := op.ePool.Get(e.Type())
				defer op.ePool.Put(eNTT)
				eNTT = eNTT.WithModLen(e.BaseModLen(), 0)

				pOp.FwdNTTTo(eNTT, e)
				op.rlweOp.MulElementTo(ctOut, ctOut, eNTT)
			}
		}
	}

	if ctOut.IsNTT() && !isNTT {
		op.rlweOp.InvNTTTo(ctOut, ctOut)
	} else if !ctOut.IsNTT() && isNTT {
		op.rlweOp.FwdNTTTo(ctOut, ctOut)
	}
}
