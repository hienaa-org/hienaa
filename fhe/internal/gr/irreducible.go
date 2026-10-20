package gr

import (
	"math/bits"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// irreduciblePoly returns a random irreducible polynomial of rank r modulo p.
func irreduciblePoly(p uint64, r int) []int64 {
	if r <= 0 {
		panic("r must be positive")
	} else if !num.IsPrime(p) {
		panic("p must be prime")
	}

	if r == 1 {
		return []int64{0, 1}
	}

	rSrc := csprng.NewUniformSamplerWithSeed(nil)
	q := num.NewModulus(p)
	f := make([]int64, r+1)
	f[r] = 1
	for {
		for i := 0; i < r; i++ {
			f[i] = rSrc.SampleN(-int64(p>>1), int64(p>>1))
		}
		if isIrreducible(f, q) {
			return f
		}
	}
}

// isIrreducible checks if f is irreducible modulo q.
func isIrreducible(f []int64, q *num.Modulus) bool {
	op := crt.NewOperator(dft.NewOtherParameters(f), []*num.Modulus{q})
	fRed := vec.Reduce(f, q)

	hPrev := op.NewPoly(dft.FormCoeff)
	hPrev.Coeffs[0][1] = 1
	op.FwdNTTTo(hPrev, hPrev)
	h, hInvNTT := op.NewPoly(dft.FormNTT), op.NewPoly(dft.FormCoeff)

	p0 := make([]uint64, len(f))

	qv := q.Value()
	for range op.Params().Rank() >> 1 {
		h.CopyFrom(hPrev)
		for b := bits.Len64(qv) - 2; b >= 0; b-- {
			op.MulTo(h, h, h)
			if (qv>>b)&1 == 1 {
				op.MulTo(h, h, hPrev)
			}
		}
		hPrev.CopyFrom(h)

		op.InvNTTTo(hInvNTT, h)

		copy(p0, fRed)
		hInvNTT.Coeffs[0][1] = num.Sub(hInvNTT.Coeffs[0][1], 1, q)

		if len(gcd(p0, hInvNTT.Coeffs[0], q)) != 1 {
			return false
		}
	}
	return true
}

func gcd(p0, p1 []uint64, mod *num.Modulus) []uint64 {
	p0, p1 = trim(p0), trim(p1)
	if len(p0) < len(p1) {
		p0, p1 = p1, p0
	}

	for len(p1) > 1 {
		inv := num.Inv(p1[len(p1)-1], mod)
		for i := len(p0) - len(p1); i >= 0; i-- {
			c := num.Mul(p0[i+len(p1)-1], inv, mod)
			if c != 0 {
				vec.MulSubTo(p0[i:i+len(p1)], p1, c, mod)
			}
		}
		p0, p1 = p1, trim(p0[:len(p1)-1])
	}

	if len(p1) == 1 {
		return p1
	}
	return p0
}

func trim(p []uint64) []uint64 {
	for len(p) > 0 && p[len(p)-1] == 0 {
		p = p[:len(p)-1]
	}
	return p
}
