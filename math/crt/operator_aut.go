package crt

import (
	"math/bits"
	"slices"
	"sync"

	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
)

// autOperator implements automorphism operations.
type autOperator interface {
	// CanAut returns whether the automorphism can be applied.
	CanAut(idx int) bool
	// Aut returns aut(e, idx).
	// Panics when the automorphism index is invalid.
	Aut(e *Element, idx int) *Element
	// AutTo computes eOut = aut(e, idx).
	// Panics when the automorphism index is invalid.
	AutTo(eOut, e *Element, idx int)
}

// pow2CyclotomicAutOperator is a [autOperator] for power-of-two cyclotomic ring.
type pow2CyclotomicAutOperator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	pool *sync.Pool
}

// newPow2CyclotomicAutOperator creates a new [pow2CyclotomicAutOperator].
func newPow2CyclotomicAutOperator(params dft.RingParameters, mod []*num.Modulus) pow2CyclotomicAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return pow2CyclotomicAutOperator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, params.Rank())
				return &v
			},
		},
	}
}

// CanAut returns whether the given automorphism index is valid.
func (op *pow2CyclotomicAutOperator) CanAut(idx int) bool {
	cycloOrd := op.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return idx%2 == 1
}

// Aut returns aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *pow2CyclotomicAutOperator) Aut(e *Element, idx int) *Element {
	eOut := NewPoly(e.Rank(), len(op.mod))
	op.AutTo(eOut, e, idx)
	return eOut
}

// AutTo computes eOut = aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *pow2CyclotomicAutOperator) AutTo(eOut, e *Element, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	switch e.Type() {
	case TypeScalar:
		eOut.CopyFrom(e)
	case TypePoly:
		if !op.CanAut(idx) {
			panic("invalid automorphism index")
		}

		cycloOrd, rank := op.params.CycloOrder(), op.params.Rank()
		idx = (idx%cycloOrd + cycloOrd) % cycloOrd

		if idx == 1 {
			eOut.CopyFrom(e)
			return
		}

		eBufPtr := op.pool.Get().(*[]uint64)
		eBuf := *eBufPtr
		defer op.pool.Put(eBufPtr)

		for i := range op.mod {
			if e.IsNTT && op.isNTTFriendly[i] {
				copy(eBuf, e.Coeffs[i])
				revShiftBits := 64 - int(num.Log2(rank))
				for j := 0; j < rank; j++ {
					jOut := ((2*j + 1) * idx) & (cycloOrd - 1)
					idxIn := int(bits.Reverse64((uint64(jOut)-1)/2) >> revShiftBits)
					idxOut := int(bits.Reverse64(uint64(j)) >> revShiftBits)
					eOut.Coeffs[i][idxOut] = eBuf[idxIn]
				}
			} else {
				clear(eBuf)
				for j := 0; j < rank; j++ {
					idxOut := (j * idx) & (cycloOrd - 1)
					if idxOut >= rank {
						eBuf[idxOut-rank] = num.Neg(e.Coeffs[i][j], op.mod[i])
					} else {
						eBuf[idxOut] = e.Coeffs[i][j]
					}
				}
				copy(eOut.Coeffs[i], eBuf)
			}
		}

		eOut.IsNTT = e.IsNTT
	}
}

func (op *pow2CyclotomicAutOperator) subOperator(idx ...int) pow2CyclotomicAutOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
	}

	return pow2CyclotomicAutOperator{
		params:        op.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		pool: op.pool,
	}
}

// anyCyclotomicAutOperator is a [autOperator] for any cyclotomic ring.
type anyCyclotomicAutOperator struct {
	params        dft.RingParameters
	cycloOrdMod   *num.Modulus
	mod           []*num.Modulus
	isNTTFriendly []bool

	reducer *CyclotomicReducer

	// primeExpMods is the prime power factors of the cyclotomic order.
	primeExpMods []*num.Modulus
	// rootExps is the generators modulo prime power factors of the cyclotomic order.
	rootExps [][]uint64
	// dims is the dimension of the hypercube structure.
	dims []int

	pool *sync.Pool
}

// newAnyCyclotomicAutOperator creates a new [anyCyclotomicAutOperator].
func newAnyCyclotomicAutOperator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) anyCyclotomicAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	primes, exps := num.Factor(params.CycloOrder())
	pExpMods := make([]*num.Modulus, len(primes))
	rootExps := make([][]uint64, len(primes))
	dims := make([]int, len(primes))
	for i := range rootExps {
		pExp := 1
		for j := 0; j < int(exps[i]); j++ {
			pExp *= primes[i]
		}
		pExpMods[i] = num.NewModulus(pExp)
		dims[i] = pExp - pExp/primes[i]
		rootExps[i] = make([]uint64, dims[i])
		root := uint64(1)
		if pExp != 2 {
			root = num.GeneratorsWithFactors(pExpMods[i], []uint64{uint64(primes[i])}, []uint64{uint64(exps[i])})[0]
		}

		if primes[i] == 2 && exps[i] > 2 {
			root = 5
		}

		rootExps[i][0] = 1
		for j := 1; j < len(rootExps[i]); j++ {
			rootExps[i][j] = num.Mul(root, rootExps[i][j-1], pExpMods[i])
		}
	}

	if primes[0] == 2 {
		if exps[0] == 1 {
			rootExps = rootExps[1:]
			pExpMods = pExpMods[1:]
			dims = dims[1:]
		} else if exps[0] > 2 {
			rootExps = append([][]uint64{rootExps[0][:len(rootExps[0])/2], {1, pExpMods[0].Value() - 1}}, rootExps[1:]...)
			pExpMods = append([]*num.Modulus{pExpMods[0], pExpMods[0]}, pExpMods[1:]...)
			dims = append([]int{dims[0] / 2, 2}, dims[1:]...)
		}
	}

	return anyCyclotomicAutOperator{
		params:        params,
		cycloOrdMod:   num.NewModulus(params.CycloOrder()),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		reducer: reducer,

		primeExpMods: pExpMods,
		rootExps:     rootExps,
		dims:         dims,

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, params.CycloOrder())
				return &v
			},
		},
	}
}

// CanAut returns whether the given automorphism index is valid.
func (op *anyCyclotomicAutOperator) CanAut(idx int) bool {
	cycloOrd := op.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return num.GCD(idx, cycloOrd) == 1
}

// Aut returns aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *anyCyclotomicAutOperator) Aut(e *Element, idx int) *Element {
	eOut := NewPoly(e.Rank(), len(op.mod))
	op.AutTo(eOut, e, idx)
	return eOut
}

// AutTo computes eOut = aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *anyCyclotomicAutOperator) AutTo(eOut, e *Element, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	switch e.Type() {
	case TypeScalar:
		eOut.CopyFrom(e)
	case TypePoly:
		if !op.CanAut(idx) {
			panic("invalid automorphism index")
		}

		cycloOrd, rank := op.params.CycloOrder(), op.params.Rank()
		idx = (idx%cycloOrd + cycloOrd) % cycloOrd

		if idx == 1 {
			eOut.CopyFrom(e)
			return
		}

		eBufPtr := op.pool.Get().(*[]uint64)
		eBuf := *eBufPtr
		defer op.pool.Put(eBufPtr)

		for i := range op.mod {
			if e.IsNTT && op.isNTTFriendly[i] {
				idxDigits := make([]int, len(op.primeExpMods))
				for cnt := 0; cnt < len(idxDigits); {
					idxRed := num.Reduce(uint64(idx), op.primeExpMods[cnt])
					if op.primeExpMods[cnt].Value()%8 == 0 {
						if idx%4 == 1 {
							idxDigits[cnt+1] = 0
						} else {
							idxDigits[cnt+1] = 1
							idxRed = num.Neg(idxRed, op.primeExpMods[cnt])
						}
						idxDigits[cnt] = slices.Index(op.rootExps[cnt], idxRed)
						cnt += 2
					} else {
						idxDigits[cnt] = slices.Index(op.rootExps[cnt], idxRed)
						cnt += 1
					}
				}

				copy(eBuf, e.Coeffs[i])

				idxOutDigits := make([]int, len(op.primeExpMods))
				for j := 0; j < rank; j++ {
					idxIn := j
					for k := 0; k < len(idxOutDigits); k++ {
						idxOutDigits[k] = (idxIn - idxDigits[k] + op.dims[k]) % op.dims[k]
						idxIn /= op.dims[k]
					}
					idxOut := idxOutDigits[len(idxDigits)-1]
					for k := len(idxOutDigits) - 2; k >= 0; k-- {
						idxOut *= op.dims[k]
						idxOut += idxOutDigits[k]
					}
					eOut.Coeffs[i][idxOut] = eBuf[j]
				}
			} else {
				clear(eBuf)
				for j := 0; j < rank; j++ {
					idxOut := num.Mul(uint64(j), uint64(idx), op.cycloOrdMod)
					eBuf[idxOut] = e.Coeffs[i][j]
				}
				op.reducer.reduceTo(eOut.Coeffs[i], eBuf, i)
			}
		}

		eOut.IsNTT = e.IsNTT
	}
}

func (op *anyCyclotomicAutOperator) subOperator(idx ...int) anyCyclotomicAutOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
	}

	return anyCyclotomicAutOperator{
		params:        op.params,
		cycloOrdMod:   op.cycloOrdMod,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		reducer: op.reducer.SubReducer(idx...),

		primeExpMods: op.primeExpMods,
		rootExps:     op.rootExps,
		dims:         op.dims,

		pool: op.pool,
	}
}

// pow2AutFixedAutOperator is a [autOperator] for power-of-two autfixed ring.
type pow2AutFixedAutOperator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	pool *sync.Pool
}

// newPow2AutFixedAutOperator creates a new [pow2AutFixedAutOperator].
func newPow2AutFixedAutOperator(params dft.RingParameters, mod []*num.Modulus) pow2AutFixedAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return pow2AutFixedAutOperator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, params.Rank())
				return &v
			},
		},
	}
}

// CanAut returns whether the given automorphism index is valid.
func (op *pow2AutFixedAutOperator) CanAut(idx int) bool {
	cycloOrd := op.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return idx%4 == 1
}

// Aut returns aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *pow2AutFixedAutOperator) Aut(e *Element, idx int) *Element {
	eOut := NewPoly(e.Rank(), len(op.mod))
	op.AutTo(eOut, e, idx)
	return eOut
}

// AutTo computes eOut = aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *pow2AutFixedAutOperator) AutTo(eOut, e *Element, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	switch e.Type() {
	case TypeScalar:
		eOut.CopyFrom(e)
	case TypePoly:
		if !op.CanAut(idx) {
			panic("invalid automorphism index")
		}

		cycloOrd, rank := op.params.CycloOrder(), op.params.Rank()
		idx = (idx%cycloOrd + cycloOrd) % cycloOrd

		if idx == 1 {
			eOut.CopyFrom(e)
			return
		}

		eBufPtr := op.pool.Get().(*[]uint64)
		eBuf := *eBufPtr
		defer op.pool.Put(eBufPtr)

		for i := range op.mod {
			if e.IsNTT && op.isNTTFriendly[i] {
				copy(eBuf, e.Coeffs[i])
				revShiftBits := 64 - int(num.Log2(rank)+1)
				for j := 0; j < rank; j++ {
					jOut := ((2*j + 1) * idx) & (cycloOrd - 1)
					idxIn := int(bits.Reverse64((uint64(jOut)-1)/2) >> revShiftBits)
					if idxIn >= rank {
						idxIn = 2*rank - 1 - idxIn
					}
					idxOut := int(bits.Reverse64(uint64(j)) >> revShiftBits)
					if idxOut >= rank {
						idxOut = 2*rank - 1 - idxOut
					}
					eOut.Coeffs[i][idxOut] = eBuf[idxIn]
				}
			} else {
				clear(eBuf)
				for j := 0; j < rank; j++ {
					idxOut := (j * idx) & (cycloOrd - 1)
					switch {
					case idxOut >= 3*rank:
						eBuf[4*rank-idxOut] = e.Coeffs[i][j]
					case idxOut >= 2*rank:
						eBuf[idxOut-2*rank] = num.Neg(e.Coeffs[i][j], op.mod[i])
					case idxOut >= rank:
						eBuf[2*rank-idxOut] = num.Neg(e.Coeffs[i][j], op.mod[i])
					default:
						eBuf[idxOut] = e.Coeffs[i][j]
					}
				}
				copy(eOut.Coeffs[i], eBuf)
			}
		}

		eOut.IsNTT = e.IsNTT
	}
}

func (op *pow2AutFixedAutOperator) subOperator(idx ...int) pow2AutFixedAutOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
	}

	return pow2AutFixedAutOperator{
		params:        op.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		pool: op.pool,
	}
}

// primeAutFixedAutOperator is a [autOperator] for prime autfixed ring.
type primeAutFixedAutOperator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	// rootPow are the power of the generator modulo the cyclotomic order.
	rootPow []uint64
	// rootPowInv are the powers of the inverse of the generator modulo the cyclotomic order.
	rootPowInv []uint64

	pool *sync.Pool
}

// newPrimeAutFixedAutOperator creates a new [primeAutFixedAutOperator].
func newPrimeAutFixedAutOperator(params dft.RingParameters, mod []*num.Modulus) primeAutFixedAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	cycloOrdMod := num.NewModulus(params.CycloOrder())
	root := num.Generators(cycloOrdMod)[0]
	rootInv := num.Inv(root, cycloOrdMod)

	rootPow := make([]uint64, params.Rank())
	rootPowInv := make([]uint64, params.Rank())
	rootPow[0] = 1
	rootPowInv[0] = 1
	for i := 1; i < params.Rank(); i++ {
		rootPow[i] = num.Mul(rootPow[i-1], root, cycloOrdMod)
		rootPowInv[i] = num.Mul(rootPowInv[i-1], rootInv, cycloOrdMod)
	}

	return primeAutFixedAutOperator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		rootPow:    rootPow,
		rootPowInv: rootPowInv,

		pool: &sync.Pool{
			New: func() any {
				v := make([]uint64, params.Rank())
				return &v
			},
		},
	}
}

// CanAut returns whether the given automorphism index is valid.
func (op *primeAutFixedAutOperator) CanAut(idx int) bool {
	cycloOrd := op.params.CycloOrder()
	idx = (idx%cycloOrd + cycloOrd) % cycloOrd
	return slices.Contains(op.rootPow, uint64(idx)) || slices.Contains(op.rootPowInv, uint64(idx))
}

// Aut returns aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *primeAutFixedAutOperator) Aut(e *Element, idx int) *Element {
	eOut := NewPoly(e.Rank(), len(op.mod))
	op.AutTo(eOut, e, idx)
	return eOut
}

// AutTo computes eOut = aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *primeAutFixedAutOperator) AutTo(eOut, e *Element, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), eOut, e)

	switch e.Type() {
	case TypeScalar:
		eOut.CopyFrom(e)
	case TypePoly:
		if !op.CanAut(idx) {
			panic("invalid automorphism index")
		}

		cycloOrd, rank := op.params.CycloOrder(), op.params.Rank()
		idx = (idx%cycloOrd + cycloOrd) % cycloOrd

		eBufPtr := op.pool.Get().(*[]uint64)
		eBuf := *eBufPtr
		defer op.pool.Put(eBufPtr)

		rotIdx := slices.Index(op.rootPow, uint64(idx))
		if rotIdx == -1 {
			rotIdx = slices.Index(op.rootPowInv, uint64(idx))
			rotIdx = (rank - rotIdx) % rank
		}

		for i := range op.mod {
			if e.IsNTT && op.isNTTFriendly[i] {
				copy(eBuf[:rank-rotIdx], e.Coeffs[i][rotIdx:])
				copy(eBuf[rank-rotIdx:], e.Coeffs[i][:rotIdx])
				copy(eOut.Coeffs[i], eBuf)
			} else {
				copy(eBuf[:rotIdx], e.Coeffs[i][rank-rotIdx:])
				copy(eBuf[rotIdx:], e.Coeffs[i][:rank-rotIdx])
				copy(eOut.Coeffs[i], eBuf)
			}
		}

		eOut.IsNTT = e.IsNTT
	}
}

func (op *primeAutFixedAutOperator) subOperator(idx ...int) primeAutFixedAutOperator {
	modCopy := make([]*num.Modulus, len(idx))
	isNTTFriendlyCopy := make([]bool, len(idx))
	for i := range idx {
		modCopy[i] = op.mod[idx[i]]
		isNTTFriendlyCopy[i] = op.isNTTFriendly[idx[i]]
	}

	return primeAutFixedAutOperator{
		params:        op.params,
		mod:           modCopy,
		isNTTFriendly: isNTTFriendlyCopy,

		rootPow:    op.rootPow,
		rootPowInv: op.rootPowInv,

		pool: op.pool,
	}
}

// noAutOperator is a no-op [autOperator].
type noAutOperator struct{}

// CanAut returns whether the given automorphism index is valid.
func (op *noAutOperator) CanAut(idx int) bool {
	return false
}

// Aut returns aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *noAutOperator) Aut(e *Element, idx int) *Element {
	panic("automorphism not supported")
}

// AutTo computes eOut = aut(e, idx).
// Panics when the automorphism index is invalid.
func (op *noAutOperator) AutTo(eOut, e *Element, idx int) {
	panic("automorphism not supported")
}
