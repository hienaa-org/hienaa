package pack

import (
	"math"
	"sync"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type pow2CyclotomicRealPacker struct {
	params dft.RingParameters

	// packLen is the packing length.
	packLen int

	// roots are the roots of unity.
	roots []complex128
	//
	group []int

	// cube is the form of the hypercube structure.
	cube []int
	// cubeGen is the corresponding generator for the hypercube structure.
	cubeGen []uint64

	pool *sync.Pool
}

// newPow2CyclotomicRealPacker creates a new [pow2CyclotomicRealPacker].
func newPow2CyclotomicRealPacker(params dft.RingParameters) *pow2CyclotomicRealPacker {
	roots := make([]complex128, params.CycloOrder())
	for i := 0; i < params.CycloOrder(); i++ {
		cos := math.Cos(2 * math.Pi * float64(i) / float64(params.CycloOrder()))
		sin := math.Sin(2 * math.Pi * float64(i) / float64(params.CycloOrder()))
		roots[i] = complex(cos, sin)
	}

	group := make([]int, params.Rank())
	for i := 0; i < params.Rank(); i++ {
		group[i] = int(num.Exp(uint64(5), uint64(i), nil)) & (params.CycloOrder() - 1)
	}

	return &pow2CyclotomicRealPacker{
		params: params,

		packLen: params.Rank() >> 1,

		roots: roots,
		group: group,

		cube:    []int{params.CycloOrder() >> 1},
		cubeGen: []uint64{5},

		pool: &sync.Pool{
			New: func() any {
				v := make([]complex128, params.Rank()>>1)
				return &v
			},
		},
	}
}

func (p *pow2CyclotomicRealPacker) Params() dft.RingParameters {
	return p.params
}

func (p *pow2CyclotomicRealPacker) PackLen() int {
	return p.packLen
}

func (p *pow2CyclotomicRealPacker) Pack(v []complex128) []float64 {
	vPack := make([]float64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

// TODO: Optimise the algorithm. This algorithm is due to ia.cr/2018/1043.
// If we have an efficient FFT implementation, we can directly use it.
func (p *pow2CyclotomicRealPacker) PackTo(vPack []float64, v []complex128) {
	if len(vPack) != p.params.Rank() || p.packLen%len(v) != 0 {
		panic("input(s) shape not consistent")
	}

	bufPtr := p.pool.Get().(*[]complex128)
	buf := *bufPtr
	defer p.pool.Put(bufPtr)
	buf = buf[:len(v)]
	clear(buf)

	for i := 0; i < len(v); i++ {
		buf[i] = v[i] / complex(float64(len(v)), 0)
	}
	for idx := int(num.Log2(len(v))); idx >= 1; idx-- {
		for i := 0; i < len(v); i += (1 << idx) {
			lenh, lenQ := (1<<idx)>>1, (1<<idx)<<2
			gap := p.params.CycloOrder() / lenQ
			for j := 0; j < lenh; j++ {
				idx1, idx2 := i+j, i+j+lenh
				rootsi := p.roots[(lenQ-(p.group[j]&(lenQ-1)))*gap]
				buf[idx1], buf[idx2] = buf[idx1]+buf[idx2], (buf[idx1]-buf[idx2])*rootsi
			}
		}
	}
	vec.RadixReverseInPlace(buf, 2)

	clear(vPack)
	for i := 0; i < len(v); i++ {
		vPack[i*p.params.Rank()/len(v)/2] = real(buf[i])
		vPack[(i+len(v))*p.params.Rank()/len(v)/2] = imag(buf[i])
	}
}

func (p *pow2CyclotomicRealPacker) UnPack(vPack []float64) []complex128 {
	v := make([]complex128, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

func (p *pow2CyclotomicRealPacker) UnPackTo(v []complex128, vPack []float64) {
	if p.packLen%len(v) != 0 || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	bufPtr := p.pool.Get().(*[]complex128)
	buf := *bufPtr
	defer p.pool.Put(bufPtr)

	for i := 0; i < len(buf); i++ {
		buf[i] = complex(vPack[i], vPack[i+len(buf)])
	}

	vec.RadixReverseInPlace(buf, 2)
	for idx := 1; idx <= int(num.Log2(len(buf))); idx++ {
		for i := 0; i < len(buf); i += (1 << idx) {
			lenh, lenQ := (1<<idx)>>1, (1<<idx)<<2
			gap := p.params.CycloOrder() / lenQ
			for j := 0; j < lenh; j++ {
				tmp := p.roots[(p.group[j]&(lenQ-1))*gap] * buf[i+j+lenh]
				buf[i+j], buf[i+j+lenh] = buf[i+j]+tmp, buf[i+j]-tmp
			}
		}
	}

	copy(v, buf[:len(v)])
}

func (p *pow2CyclotomicRealPacker) Cube() []int {
	return p.cube
}

func (p *pow2CyclotomicRealPacker) CubeGen() []uint64 {
	return p.cubeGen
}

func (p *pow2CyclotomicRealPacker) RotIdxToAutIdx(idx []int) int {
	if len(idx) != 1 {
		panic("input(s) shape not consistent")
	}
	return int(num.Exp(5, uint64(idx[0]), nil)) & (p.params.CycloOrder() - 1)
}

type trivialRealPacker struct {
	params dft.RingParameters

	packLen int

	cube    []int
	cubeGen []uint64
}

func newTrivialRealPacker(params dft.RingParameters) *trivialRealPacker {
	return &trivialRealPacker{
		params: params,

		packLen: 1,

		cube:    []int{1},
		cubeGen: []uint64{1},
	}
}

func (p *trivialRealPacker) Params() dft.RingParameters {
	return p.params
}

func (p *trivialRealPacker) PackLen() int {
	return p.packLen
}

func (p *trivialRealPacker) Pack(v []complex128) []float64 {
	vPack := make([]float64, p.params.Rank())
	p.PackTo(vPack, v)
	return vPack
}

func (p *trivialRealPacker) PackTo(vPack []float64, v []complex128) {
	if len(vPack) != p.params.Rank() || len(v) != p.packLen {
		panic("input(s) shape not consistent")
	} else if imag(v[0]) != 0 {
		panic("input(s) shape not consistent")
	}

	clear(vPack)
	vPack[0] = real(v[0])
}

func (p *trivialRealPacker) UnPack(vPack []float64) []complex128 {
	v := make([]complex128, p.packLen)
	p.UnPackTo(v, vPack)
	return v
}

func (p *trivialRealPacker) UnPackTo(v []complex128, vPack []float64) {
	if len(v) != p.packLen || len(vPack) != p.params.Rank() {
		panic("input(s) shape not consistent")
	}

	v[0] = complex(vPack[0], 0)
}

func (p *trivialRealPacker) Cube() []int {
	return p.cube
}

func (p *trivialRealPacker) CubeGen() []uint64 {
	return p.cubeGen
}

func (p *trivialRealPacker) RotIdxToAutIdx(idx []int) int {
	return 1
}
