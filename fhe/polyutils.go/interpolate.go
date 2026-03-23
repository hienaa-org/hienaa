package polyutils

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

func nu(m, p int) int {
	res := 0
	for m%p == 0 {
		m /= p
		res++
	}
	return res
}

func mu(p, e int) int {
	i := p
	cnt := 0
	for cnt < e {
		cnt += nu(i, p)
		i += p
	}
	return i - p
}

type Interpolater struct {
	p   int
	r   int
	pr  *num.Modulus
	buf *sync.Pool
}

func NewInterpolater(p, r int) *Interpolater {
	return &Interpolater{
		p:  p,
		r:  r,
		pr: num.NewModulus(num.Exp(uint64(p), uint64(r), nil)),
		buf: &sync.Pool{
			New: func() any {
				return make([]uint64, mu(p, r))
			},
		},
	}
}

// finiteDiffTo computes the finite difference of the polynomial that passes through the given points.
func (it *Interpolater) finiteDiffTo(res, x, y []uint64) {
	inLen := len(x)

	if inLen > mu(it.p, it.r) {
		panic("invalid input length")
	} else if len(res) != len(y) || len(res) != len(x) {
		panic("inconsistent input length")
	}

	tmp1 := it.buf.Get().([]uint64)
	tmp2 := it.buf.Get().([]uint64)
	defer it.buf.Put(tmp1)
	defer it.buf.Put(tmp2)
	tmp1 = tmp1[:inLen]
	tmp2 = tmp2[:inLen]

	clear(res)
	copy(tmp1, y)
	clear(tmp2)

	for i := 1; i < inLen; i++ {
		res[i-1] = tmp1[0]
		for j := 0; j < inLen-i; j++ {
			resj := num.Sub(tmp1[j+1], tmp1[j], it.pr)
			denom := num.Sub(x[j+i], x[j], it.pr)

			for denom%uint64(it.p) == 0 {
				denom /= uint64(it.p)
				if resj%uint64(it.p) != 0 {
					panic("invalid input")
				}
				resj /= uint64(it.p)
			}

			tmp2[j] = num.Mul(resj, num.Inv(denom, it.pr), it.pr)
		}
		copy(tmp1, tmp2)
	}
	res[inLen-1] = tmp1[0]
}

// Interpolate interpolates the polynomial that passes through the given points.
func (it *Interpolater) Interpolate(x, y []uint64) []uint64 {
	if len(x) != len(y) {
		panic("inconsistent input length")
	} else if len(x) > mu(it.p, it.r) {
		panic("invalid input length")
	}

	inLen := len(x)
	diff := it.buf.Get().([]uint64)
	defer it.buf.Put(diff)
	diff = diff[:inLen]

	it.finiteDiffTo(diff, x, y)

	res := make([]uint64, inLen)

	basis := it.buf.Get().([]uint64)
	defer it.buf.Put(basis)
	basis = basis[:inLen]
	clear(basis)

	basis[0] = 1
	for i := 0; i < inLen; i++ {
		vec.MulAddScalarTo(res, basis, diff[i], it.pr)

		if i < inLen-1 {
			basis[i+1] = basis[i]
			for j := i; j > 0; j-- {
				basis[j] = num.Sub(basis[j-1], num.Mul(basis[j], x[i], it.pr), it.pr)
			}
			basis[0] = num.Mul(basis[0], it.pr.Value()-x[i], it.pr)
		}
	}

	return res
}

// DigitExtractPoly returns the polynomial that extracts the i-th digit of the plaintext.
func (it *Interpolater) DigitExtractPoly(i int) []uint64 {
	switch it.p {
	case 2:
		deg := i + 2

		x := it.buf.Get().([]uint64)
		y := it.buf.Get().([]uint64)
		defer it.buf.Put(x)
		defer it.buf.Put(y)
		x = x[:deg]
		y = y[:deg]

		for j := 0; j < deg; j++ {
			x[j] = uint64(j)
			y[j] = uint64(j & 1)
		}

		res := it.Interpolate(x, y)
		for j := 1; j < deg; j += 2 {
			res[j] = 0
		}
		return res
	default:
		deg := (it.p-1)*(i-1) + 2

		x := it.buf.Get().([]uint64)
		y := it.buf.Get().([]uint64)
		defer it.buf.Put(x)
		defer it.buf.Put(y)
		x = x[:deg]
		y = y[:deg]

		for j := 0; j < deg; j++ {
			x[j] = uint64(j)
			y[j] = num.Reduce((j+it.p/2)%it.p-it.p/2, it.pr)
		}

		res := it.Interpolate(x, y)
		for j := 0; j < deg; j += 2 {
			res[j] = 0
		}
		return res
	}
}
