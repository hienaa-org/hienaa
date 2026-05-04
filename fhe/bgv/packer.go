package bgv

import (
	"math"
	"slices"

	"github.com/hienaa-org/hienaa/fhe/internal/heint"
	pack "github.com/hienaa-org/hienaa/fhe/internal/pack"
	"github.com/hienaa-org/hienaa/fhe/rlwe"
	"github.com/hienaa-org/hienaa/math/num"
)

// Packer packs a vector of uint64 into [*Plaintext].
type Packer struct {
	params rlwe.Parameters
	pack   pack.IntPacker
	ecd    *heint.Encoder
}

// NewPacker creates a new [Packer].
func NewPacker(params rlwe.Parameters, msgMod *num.Modulus) *Packer {
	return &Packer{
		params: params,
		pack:   pack.NewIntPacker(params.RingParams(), msgMod),
		ecd:    heint.NewEncoder(params, msgMod),
	}
}

// Params returns the parameters.
func (p *Packer) Params() rlwe.Parameters {
	return p.params
}

// PackLen returns the length of packable vector.
func (p *Packer) PackLen() int {
	return p.pack.PackLen()
}

// Pack packs v.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) Pack(v []uint64) Plaintext {
	ptOut := NewPoly(p.params.Rank())
	p.PackTo(ptOut, v)
	return ptOut
}

// PackTo packs v to ptOut.
// Panics when the message length does not divide the packing length,
// or when the output modulus length is larger than the number of moduli.
func (p *Packer) PackTo(ptOut Plaintext, v []uint64) {
	if p.PackLen()%len(v) != 0 {
		panic("len(v) must divide PackLen")
	} else if ptOut.Rank() != p.params.Rank() {
		panic("inconsistent input(s)")
	}

	p.pack.PackTo(ptOut, v)
}

// UnPack unpacks pt.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPack(pt Plaintext) []uint64 {
	vOut := make([]uint64, p.pack.PackLen())
	p.UnPackTo(vOut, pt)
	return vOut
}

// UnPack unpacks pt to vOut.
// Panics when the moduli length of the input plaintext is larger than the number of moduli,
// or when the message length does not divide the packing length, or when the output modulus length is larger than the number of moduli.
func (p *Packer) UnPackTo(vOut []uint64, pt Plaintext) {
	if p.PackLen()%len(vOut) != 0 {
		panic("len(vOut) must divide PackLen")
	} else if pt.Rank() != p.params.Rank() {
		panic("inconsistent input(s)")
	}

	p.pack.UnPackTo(vOut, pt)
}

// Cube returns the form of the hypercube structure.
func (p *Packer) Cube() []int {
	return p.pack.Cube()
}

// CubeGen returns the corresponding generator for the hypercube structure.
func (p *Packer) CubeGen() []uint64 {
	return p.pack.CubeGen()
}

// RotIdxToAutIdx converts a rotation index to an automorphism index.
func (p *Packer) RotIdxToAutIdx(idx []int) int {
	return p.pack.RotIdxToAutIdx(idx)
}

// Encode encodes a [Plaintext] into a [*rlwe.Element].
func (p *Packer) Encode(eIn Plaintext, hasAux, isNTT bool) *rlwe.Element {
	return p.ecd.Encode(eIn, hasAux, isNTT)
}

// EncodeCustom encodes a [Plaintext] into a [*rlwe.Element] with custom parameters.
func (p *Packer) EncodeCustom(eIn Plaintext, baseLen, auxLen int, isNTT bool) *rlwe.Element {
	return p.ecd.EncodeCustom(eIn, baseLen, auxLen, isNTT)
}

// EncodeTo encodes a [Plaintext] into a [*rlwe.Element].
func (p *Packer) EncodeTo(eOut *rlwe.Element, eIn Plaintext, isNTT bool) {
	p.ecd.EncodeTo(eOut, eIn, isNTT)
}

// Decode decodes a [*rlwe.Element] into a [Plaintext].
func (p *Packer) Decode(e *rlwe.Element) Plaintext {
	return p.ecd.Decode(e)
}

// DecodeTo decodes a [*rlwe.Element] into a [Plaintext].
func (p *Packer) DecodeTo(eOut Plaintext, e *rlwe.Element) {
	p.ecd.DecodeTo(eOut, e)
}

// GenPlainMatrix generates a plaintext matrix from a matrix.
func (p *Packer) GenPlainMatrix(mat map[[2]int]uint64, dim int, bsgsRatio float64) *rlwe.PlainMatrix {
	return p.GenPlainMatrixCustom(mat, dim, len(p.params.BaseModulus()), len(p.params.AuxModulus()), bsgsRatio)
}

// getBSGSBase returns the base for the BSGS matrix multiplication algorithm.
// If bsgsRatio = 0, we do not use the BSGS algorithm at all.
func (p *Packer) GenPlainMatrixCustom(mat map[[2]int]uint64, dim, baseLen, auxLen int, bsgsRatio float64) *rlwe.PlainMatrix {
	if bsgsRatio < 0 {
		panic("bsgsRatio must be positive")
	}

	// Compute the BSGS parameters.
	cubeLen := len(p.Cube())
	cubeRot := make([][]int, cubeLen)
	for i := range cubeRot {
		cubeRot[i] = make([]int, 0)
	}

	usedIdx := make([]int, 0)
	tmpIdx := make([]int, cubeLen)
	tmpRotIdx := make([]int, cubeLen)
	for i := 0; i < dim; i++ {
		idx := i
		for j := range tmpRotIdx {
			tmpRotIdx[j] = idx % p.Cube()[j]
			idx /= p.Cube()[j]
		}

		check := false
		for j := 0; j < dim; j++ {
			idx = j
			for k := range tmpIdx {
				tmpIdx[k] = idx % p.Cube()[k]
				idx /= p.Cube()[k]
			}

			idx = (tmpIdx[cubeLen-1] - tmpRotIdx[cubeLen-1] + p.Cube()[cubeLen-1]) % p.Cube()[cubeLen-1]
			for k := cubeLen - 2; k >= 0; k-- {
				idx *= p.Cube()[k]
				idx += (tmpIdx[k] - tmpRotIdx[k] + p.Cube()[k]) % p.Cube()[k]
			}

			_, checkIdx := mat[[2]int{j, idx}]
			if checkIdx {
				check = true
				break
			}
		}

		if check {
			for k := range tmpRotIdx {
				if !slices.Contains(cubeRot[k], tmpRotIdx[k]) {
					cubeRot[k] = append(cubeRot[k], tmpRotIdx[k])
				}
			}

			usedIdx = append(usedIdx, i)
		}
	}

	// Compute the optimal giant step base.
	gsBase := make([]int, cubeLen)
	for i := 0; i < cubeLen; i++ {
		slices.Sort(cubeRot[i])
		closestRatio := math.MaxFloat64
		for j := 2; j <= p.Cube()[i]; j++ {
			babyStep := make([]int, 0)
			giantStep := make([]int, 0)
			for k := range cubeRot[i] {
				bs, gs := cubeRot[i][k]%j, cubeRot[i][k]/j
				if !slices.Contains(babyStep, bs) {
					babyStep = append(babyStep, bs)
				}
				if !slices.Contains(giantStep, gs) {
					giantStep = append(giantStep, gs)
				}
			}
			babyRot := len(babyStep)
			giantRot := len(giantStep)

			ratio := float64(babyRot) / float64(giantRot)
			if math.Abs(ratio-bsgsRatio) < math.Abs(closestRatio-bsgsRatio) {
				closestRatio = ratio
				gsBase[i] = j
			}
		}
	}

	// Generate the diagonal and BSGS parameters.
	vTmp := make([]uint64, dim)
	packTmp := make([]uint64, p.params.Rank())
	bsTmp := make([]int, cubeLen)
	gsTmp := make([]int, cubeLen)
	bs := make([]int, 0)
	gs := make([]int, 0)
	diag := make(map[int]map[int]*rlwe.Element)
	cycloOrd := num.NewModulus(p.params.RingParams().CycloOrder())
	for _, idx := range usedIdx {
		for j := range tmpRotIdx {
			tmpRotIdx[j] = idx % p.Cube()[j]
			idx /= p.Cube()[j]
		}

		for j := 0; j < dim; j++ {
			rotIdx := j
			for k := range tmpIdx {
				tmpIdx[k] = rotIdx % p.Cube()[k]
				rotIdx /= p.Cube()[k]
			}

			rotIdx = (tmpIdx[cubeLen-1] - tmpRotIdx[cubeLen-1] + p.Cube()[cubeLen-1]) % p.Cube()[cubeLen-1]
			for k := cubeLen - 2; k >= 0; k-- {
				rotIdx *= p.Cube()[k]
				rotIdx += (tmpIdx[k] - tmpRotIdx[k] + p.Cube()[k]) % p.Cube()[k]
			}

			val, checkIdx := mat[[2]int{j, rotIdx}]
			if checkIdx {
				vTmp[j] = val
			} else {
				vTmp[j] = 0
			}
		}

		for j := range tmpRotIdx {
			bsTmp[j] = tmpRotIdx[j] % gsBase[j]
			gsTmp[j] = (tmpRotIdx[j] / gsBase[j]) * gsBase[j]
		}

		bsIdx := p.RotIdxToAutIdx(bsTmp)
		gsIdx := p.RotIdxToAutIdx(gsTmp)

		if !slices.Contains(bs, bsIdx) {
			bs = append(bs, bsIdx)
		}
		if !slices.Contains(gs, gsIdx) {
			gs = append(gs, gsIdx)
		}

		gsInv := int(num.Inv(uint64(gsIdx), cycloOrd))
		if diag[gsIdx] == nil {
			diag[gsIdx] = make(map[int]*rlwe.Element)
		}
		p.pack.PackTo(packTmp, vTmp)
		diag[gsIdx][bsIdx] = p.EncodeCustom(packTmp, baseLen, auxLen, true)
		p.ecd.PlainOperator().AutTo(diag[gsIdx][bsIdx], diag[gsIdx][bsIdx], gsInv)
	}

	return rlwe.NewPlainMatrix(diag, rlwe.BSGSParams{
		BabyStep:  bs,
		GiantStep: gs,
	})
}
