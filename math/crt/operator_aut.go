package crt

import (
	"math/bits"
	"slices"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

type autType byte

const (
	typePow2CyclotomicAut autType = iota
	typeAnyCyclotomicAut
	typePow2AutFixedAut
	typePrimeAutFixedAut
	typeNoAut
)

func autTypeOf(params dft.RingParameters) autType {
	switch params.RingType() {
	case dft.TypeCyclotomic:
		if num.IsPowerOfTwo(params.CycloIndex()) {
			return typePow2CyclotomicAut
		}
		return typeAnyCyclotomicAut

	case dft.TypeAutFixed:
		if num.IsPowerOfTwo(params.CycloIndex()) {
			return typePow2AutFixedAut
		}
		return typePrimeAutFixedAut
	}

	return typeNoAut
}

// pow2CyclotomicAutOperator is an AutOperator for [typePow2CyclotomicAut].
type pow2CyclotomicAutOperator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	pool *pool.Pool[*[]uint64]
}

// newPow2CyclotomicAutOperator creates a new [pow2CyclotomicAutOperator].
func newPow2CyclotomicAutOperator(params dft.RingParameters, mod []*num.Modulus) *pow2CyclotomicAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &pow2CyclotomicAutOperator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.Rank())
			return &v
		}),
	}
}

func (op *pow2CyclotomicAutOperator) canAut(idx int) bool {
	idx64 := num.Reduce(idx, uint64(op.params.CycloIndex()))
	return idx64%2 == 1
}

func (op *pow2CyclotomicAutOperator) autTo(pOut, p *Poly, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if !op.canAut(idx) {
		panic("invalid automorphism index")
	}

	cycloIdx64, rank64 := uint64(op.params.CycloIndex()), uint64(op.params.Rank())
	idx64 := num.Reduce(idx, cycloIdx64)

	if idx64 == 1 {
		pOut.CopyFrom(p)
		return
	}

	pBufPtr := op.pool.Get()
	pBuf := *pBufPtr
	defer op.pool.Put(pBufPtr)

	for i := range op.mod {
		pRow, pOutRow := p.Coeffs[i], pOut.Coeffs[i]
		if p.Form == dft.FormNTT && op.isNTTFriendly[i] {
			copy(pBuf, p.Coeffs[i])
			revShiftBits := 64 - uint64(num.Log2(rank64))
			for j := uint64(0); j < rank64; j++ {
				jOut := ((2*j + 1) * idx64) & (cycloIdx64 - 1)
				idxIn := bits.Reverse64((jOut-1)/2) >> revShiftBits
				idxOut := bits.Reverse64(j) >> revShiftBits
				pOutRow[idxOut] = pBuf[idxIn]
			}
		} else {
			clear(pBuf)
			var idxOut uint64
			for j := uint64(0); j < rank64; j++ {
				if idxOut >= rank64 {
					pBuf[idxOut-rank64] = num.Neg(pRow[j], op.mod[i])
				} else {
					pBuf[idxOut] = pRow[j]
				}
				idxOut = (idxOut + idx64) & (cycloIdx64 - 1)
			}
			copy(pOutRow, pBuf)
		}
	}

	pOut.Form = p.Form
}

func (op *pow2CyclotomicAutOperator) withModIdx(idx ...int) *pow2CyclotomicAutOperator {
	return &pow2CyclotomicAutOperator{
		params:        op.params,
		mod:           vec.Gather(op.mod, idx...),
		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),

		pool: op.pool,
	}
}

func (op *pow2CyclotomicAutOperator) slice(lo, hi int) *pow2CyclotomicAutOperator {
	return &pow2CyclotomicAutOperator{
		params:        op.params,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],

		pool: op.pool,
	}
}

func (op *pow2CyclotomicAutOperator) append(op0 *pow2CyclotomicAutOperator) *pow2CyclotomicAutOperator {
	return &pow2CyclotomicAutOperator{
		params:        op.params,
		mod:           slices.Concat(op.mod, op0.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),

		pool: op.pool,
	}
}

func (op *pow2CyclotomicAutOperator) appendTmpModulus(mod *num.Modulus) *pow2CyclotomicAutOperator {
	return &pow2CyclotomicAutOperator{
		params:        op.params,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),

		pool: op.pool,
	}
}

// anyCyclotomicAutOperator is an AutOperator for [typeAnyCyclotomicAut].
type anyCyclotomicAutOperator struct {
	params        dft.RingParameters
	cycloIdxMod   *num.Modulus
	mod           []*num.Modulus
	isNTTFriendly []bool

	reducer *CyclotomicReducer

	// primeExpMods is the prime power factors of the cyclotomic index.
	primeExpMods []*num.Modulus
	// rootExps is the generators modulo prime power factors of the cyclotomic index.
	rootExps [][]uint64
	// dims is the dimension of the hypercube structure.
	dims []int

	pool *pool.Pool[*[]uint64]
}

// newAnyCyclotomicAutOperator creates a new [anyCyclotomicAutOperator].
func newAnyCyclotomicAutOperator(params dft.RingParameters, mod []*num.Modulus, reducer *CyclotomicReducer) *anyCyclotomicAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range isNTTFriendly {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	primes, exps := num.Factor(uint64(params.CycloIndex()))
	pExpMods := make([]*num.Modulus, len(primes))
	rootExps := make([][]uint64, len(primes))
	dims := make([]int, len(primes))
	for i := range rootExps {
		pExp := uint64(1)
		for range exps[i] {
			pExp *= primes[i]
		}
		pExpMods[i] = num.NewModulus(pExp)

		dims[i] = int(pExp - pExp/primes[i])
		rootExps[i] = make([]uint64, dims[i])
		root := uint64(1)
		if pExp != 2 {
			root = num.Generators(pExp)[0]
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

	return &anyCyclotomicAutOperator{
		params:        params,
		cycloIdxMod:   num.NewModulus(params.CycloIndex()),
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		reducer: reducer,

		primeExpMods: pExpMods,
		rootExps:     rootExps,
		dims:         dims,

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.CycloIndex())
			return &v
		}),
	}
}

func (op *anyCyclotomicAutOperator) canAut(idx int) bool {
	cycloIdx64 := uint64(op.params.CycloIndex())
	idx64 := num.Reduce(idx, cycloIdx64)
	return num.GCD(idx64, cycloIdx64) == 1
}

func (op *anyCyclotomicAutOperator) autTo(pOut, p *Poly, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if !op.canAut(idx) {
		panic("invalid automorphism index")
	}

	cycloIdx64 := uint64(op.params.CycloIndex())
	idx64 := num.Reduce(idx, cycloIdx64)

	if idx64 == 1 {
		pOut.CopyFrom(p)
		return
	}

	pBufPtr := op.pool.Get()
	pBuf := *pBufPtr
	defer op.pool.Put(pBufPtr)

	for i := range op.mod {
		pRow, pOutRow := p.Coeffs[i], pOut.Coeffs[i]
		if p.Form == dft.FormNTT && op.isNTTFriendly[i] {
			idxDigits := make([]int, len(op.primeExpMods))
			for cnt := 0; cnt < len(idxDigits); {
				idxRed := num.Reduce(idx64, op.primeExpMods[cnt])
				if op.primeExpMods[cnt].Value()%8 == 0 {
					if idx64%4 == 1 {
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

			copy(pBuf, pRow)

			idxOutDigits := make([]int, len(op.primeExpMods))
			for j := 0; j < op.params.Rank(); j++ {
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
				pOutRow[idxOut] = pBuf[j]
			}
		} else {
			clear(pBuf)
			var idxOut uint64
			for j := uint64(0); j < uint64(op.params.Rank()); j++ {
				pBuf[idxOut] = pRow[j]
				idxOut = num.Add(idxOut, idx64, op.cycloIdxMod)
			}
			op.reducer.reduceTo(pOutRow, pBuf, i)
		}
	}

	pOut.Form = p.Form
}

func (op *anyCyclotomicAutOperator) withModIdx(idx ...int) *anyCyclotomicAutOperator {
	return &anyCyclotomicAutOperator{
		params:        op.params,
		cycloIdxMod:   op.cycloIdxMod,
		mod:           vec.Gather(op.mod, idx...),
		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),

		reducer: op.reducer.WithModIdx(idx...),

		primeExpMods: op.primeExpMods,
		rootExps:     op.rootExps,
		dims:         op.dims,

		pool: op.pool,
	}
}

func (op *anyCyclotomicAutOperator) slice(lo, hi int) *anyCyclotomicAutOperator {
	return &anyCyclotomicAutOperator{
		params:        op.params,
		cycloIdxMod:   op.cycloIdxMod,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],

		reducer: op.reducer.Slice(lo, hi),

		primeExpMods: op.primeExpMods,
		rootExps:     op.rootExps,
		dims:         op.dims,

		pool: op.pool,
	}
}

func (op *anyCyclotomicAutOperator) append(op0 *anyCyclotomicAutOperator) *anyCyclotomicAutOperator {
	return &anyCyclotomicAutOperator{
		params:        op.params,
		cycloIdxMod:   op.cycloIdxMod,
		mod:           slices.Concat(op.mod, op0.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),

		reducer: op.reducer.Append(op0.reducer),

		primeExpMods: op.primeExpMods,
		rootExps:     op.rootExps,
		dims:         op.dims,

		pool: op.pool,
	}
}

func (op *anyCyclotomicAutOperator) appendTmpModulus(mod *num.Modulus) *anyCyclotomicAutOperator {
	return &anyCyclotomicAutOperator{
		params:        op.params,
		cycloIdxMod:   op.cycloIdxMod,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),

		reducer: op.reducer.AppendTmpModulus(mod),

		primeExpMods: op.primeExpMods,
		rootExps:     op.rootExps,
		dims:         op.dims,

		pool: op.pool,
	}
}

// pow2AutFixedAutOperator is an AutOperator for [typePow2AutFixedAut].
type pow2AutFixedAutOperator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	pool *pool.Pool[*[]uint64]
}

// newPow2AutFixedAutOperator creates a new [pow2AutFixedAutOperator].
func newPow2AutFixedAutOperator(params dft.RingParameters, mod []*num.Modulus) *pow2AutFixedAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	return &pow2AutFixedAutOperator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.Rank())
			return &v
		}),
	}
}

func (op *pow2AutFixedAutOperator) canAut(idx int) bool {
	idx64 := num.Reduce(idx, uint64(op.params.CycloIndex()))
	return idx64%4 == 1
}

func (op *pow2AutFixedAutOperator) autTo(pOut, p *Poly, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if !op.canAut(idx) {
		panic("invalid automorphism index")
	}

	cycloIdx64, rank64 := uint64(op.params.CycloIndex()), uint64(op.params.Rank())
	idx64 := num.Reduce(idx, cycloIdx64)

	if idx64 == 1 {
		pOut.CopyFrom(p)
		return
	}

	pBufPtr := op.pool.Get()
	pBuf := *pBufPtr
	defer op.pool.Put(pBufPtr)

	for i := range op.mod {
		pRow, pOutRow := p.Coeffs[i], pOut.Coeffs[i]
		if p.Form == dft.FormNTT && op.isNTTFriendly[i] {
			copy(pBuf, pRow)
			revShiftBits := 64 - uint64(num.Log2(rank64)+1)
			for j := uint64(0); j < rank64; j++ {
				jOut := ((2*j + 1) * idx64) & (cycloIdx64 - 1)
				idxIn := bits.Reverse64((jOut-1)/2) >> revShiftBits
				if idxIn >= rank64 {
					idxIn = 2*rank64 - 1 - idxIn
				}
				idxOut := bits.Reverse64(j) >> revShiftBits
				if idxOut >= rank64 {
					idxOut = 2*rank64 - 1 - idxOut
				}
				pOutRow[idxOut] = pBuf[idxIn]
			}
		} else {
			clear(pBuf)
			var idxOut uint64
			for j := uint64(0); j < rank64; j++ {
				switch {
				case idxOut >= 3*rank64:
					pBuf[4*rank64-idxOut] = pRow[j]
				case idxOut >= 2*rank64:
					pBuf[idxOut-2*rank64] = num.Neg(pRow[j], op.mod[i])
				case idxOut >= rank64:
					pBuf[2*rank64-idxOut] = num.Neg(pRow[j], op.mod[i])
				default:
					pBuf[idxOut] = pRow[j]
				}
				idxOut = (idxOut + idx64) & (cycloIdx64 - 1)
			}
			copy(pOutRow, pBuf)
		}
	}

	pOut.Form = p.Form
}

func (op *pow2AutFixedAutOperator) withModIdx(idx ...int) *pow2AutFixedAutOperator {
	return &pow2AutFixedAutOperator{
		params:        op.params,
		mod:           vec.Gather(op.mod, idx...),
		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),

		pool: op.pool,
	}
}

func (op *pow2AutFixedAutOperator) slice(lo, hi int) *pow2AutFixedAutOperator {
	return &pow2AutFixedAutOperator{
		params:        op.params,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],

		pool: op.pool,
	}
}

func (op *pow2AutFixedAutOperator) append(op0 *pow2AutFixedAutOperator) *pow2AutFixedAutOperator {
	return &pow2AutFixedAutOperator{
		params:        op.params,
		mod:           slices.Concat(op.mod, op0.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),

		pool: op.pool,
	}
}

func (op *pow2AutFixedAutOperator) appendTmpModulus(mod *num.Modulus) *pow2AutFixedAutOperator {
	return &pow2AutFixedAutOperator{
		params:        op.params,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),

		pool: op.pool,
	}
}

// primeAutFixedAutOperator is an AutOperator for [typePrimeAutFixedAut].
type primeAutFixedAutOperator struct {
	params        dft.RingParameters
	mod           []*num.Modulus
	isNTTFriendly []bool

	// rootPow are the power of the generator modulo the cyclotomic index.
	rootPow []uint64
	// rootPowInv are the powers of the inverse of the generator modulo the cyclotomic index.
	rootPowInv []uint64

	pool *pool.Pool[*[]uint64]
}

// newPrimeAutFixedAutOperator creates a new [primeAutFixedAutOperator].
func newPrimeAutFixedAutOperator(params dft.RingParameters, mod []*num.Modulus) *primeAutFixedAutOperator {
	isNTTFriendly := make([]bool, len(mod))
	for i := range mod {
		isNTTFriendly[i] = dft.IsNTTFriendly(params, mod[i])
	}

	cycloIdx64 := uint64(params.CycloIndex())
	root := num.Generators(cycloIdx64)[0]
	rootInv := num.Inv(root, cycloIdx64)

	rootPow := make([]uint64, params.Rank())
	rootPowInv := make([]uint64, params.Rank())
	rootPow[0] = 1
	rootPowInv[0] = 1
	for i := 1; i < params.Rank(); i++ {
		rootPow[i] = num.Mul(rootPow[i-1], root, cycloIdx64)
		rootPowInv[i] = num.Mul(rootPowInv[i-1], rootInv, cycloIdx64)
	}

	return &primeAutFixedAutOperator{
		params:        params,
		mod:           mod,
		isNTTFriendly: isNTTFriendly,

		rootPow:    rootPow,
		rootPowInv: rootPowInv,

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.Rank())
			return &v
		}),
	}
}

func (op *primeAutFixedAutOperator) canAut(idx int) bool {
	idx64 := num.Reduce(idx, uint64(op.params.CycloIndex()))
	return slices.Contains(op.rootPow, idx64) || slices.Contains(op.rootPowInv, idx64)
}

func (op *primeAutFixedAutOperator) autTo(pOut, p *Poly, idx int) {
	isUnaryOperable(op.params.Rank(), len(op.mod), pOut, p)

	if !op.canAut(idx) {
		panic("invalid automorphism index")
	}

	cycloIdx64, rank := uint64(op.params.CycloIndex()), op.params.Rank()
	idx64 := num.Reduce(idx, cycloIdx64)

	pBufPtr := op.pool.Get()
	pBuf := *pBufPtr
	defer op.pool.Put(pBufPtr)

	rotIdx := slices.Index(op.rootPow, idx64)
	if rotIdx == -1 {
		rotIdx = slices.Index(op.rootPowInv, idx64)
		rotIdx = (op.params.Rank() - rotIdx) % op.params.Rank()
	}

	for i := range op.mod {
		if p.Form == dft.FormNTT && op.isNTTFriendly[i] {
			copy(pBuf[:rank-rotIdx], p.Coeffs[i][rotIdx:])
			copy(pBuf[rank-rotIdx:], p.Coeffs[i][:rotIdx])
			copy(pOut.Coeffs[i], pBuf)
		} else {
			copy(pBuf[:rotIdx], p.Coeffs[i][rank-rotIdx:])
			copy(pBuf[rotIdx:], p.Coeffs[i][:rank-rotIdx])
			copy(pOut.Coeffs[i], pBuf)
		}
	}

	pOut.Form = p.Form
}

func (op *primeAutFixedAutOperator) withModIdx(idx ...int) *primeAutFixedAutOperator {
	return &primeAutFixedAutOperator{
		params:        op.params,
		mod:           vec.Gather(op.mod, idx...),
		isNTTFriendly: vec.Gather(op.isNTTFriendly, idx...),

		rootPow:    op.rootPow,
		rootPowInv: op.rootPowInv,

		pool: op.pool,
	}
}

func (op *primeAutFixedAutOperator) slice(lo, hi int) *primeAutFixedAutOperator {
	return &primeAutFixedAutOperator{
		params:        op.params,
		mod:           op.mod[lo:hi],
		isNTTFriendly: op.isNTTFriendly[lo:hi],

		rootPow:    op.rootPow,
		rootPowInv: op.rootPowInv,

		pool: op.pool,
	}
}

func (op *primeAutFixedAutOperator) append(op0 *primeAutFixedAutOperator) *primeAutFixedAutOperator {
	return &primeAutFixedAutOperator{
		params:        op.params,
		mod:           slices.Concat(op.mod, op0.mod),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, op0.isNTTFriendly),

		rootPow:    op.rootPow,
		rootPowInv: op.rootPowInv,

		pool: op.pool,
	}
}

func (op *primeAutFixedAutOperator) appendTmpModulus(mod *num.Modulus) *primeAutFixedAutOperator {
	return &primeAutFixedAutOperator{
		params:        op.params,
		mod:           slices.Concat(op.mod, []*num.Modulus{mod}),
		isNTTFriendly: slices.Concat(op.isNTTFriendly, []bool{false}),

		rootPow:    op.rootPow,
		rootPowInv: op.rootPowInv,

		pool: op.pool,
	}
}
