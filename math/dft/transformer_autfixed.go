package dft

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// pow2AutFixedTransformer is a transformer for power-of-two conjugate invariant ring.
type pow2AutFixedTransformer struct {
	params RingParameters
	mod    *num.Modulus

	// tw is the twiddle factor for NTT.
	tw []uint64
	// twS is the Shoup form of tw.
	twS []uint64
	// twInv is the twiddle factor for InvNTT.
	twInv []uint64
	// twInvS is the Shoup form of twInv.
	twInvS []uint64

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	pool *pool.Pool[*[]uint64]
}

// newPow2AutFixedTransformer creates a new [pow2AutFixedTransformer].
func newPow2AutFixedTransformer(params RingParameters, mod *num.Modulus) *pow2AutFixedTransformer {
	root := num.Generators(mod.Value())

	twLarge := make([]uint64, 2*params.rank)
	twInvLarge := make([]uint64, 2*params.rank)
	twLarge[0], twLarge[1] = 1, num.NthRoot(uint64(params.cycloIdx), root, mod.Value())
	twInvLarge[0], twInvLarge[1] = 1, num.Inv(twLarge[1], mod)
	for i := 2; i < 2*params.rank; i++ {
		twLarge[i] = num.Mul(twLarge[i-1], twLarge[1], mod)
		twInvLarge[i] = num.Mul(twInvLarge[i-1], twInvLarge[1], mod)
	}
	vec.RadixReverseInPlace(twLarge, 2)
	vec.RadixReverseInPlace(twInvLarge, 2)

	tw := make([]uint64, params.rank)
	twInv := make([]uint64, params.rank)
	tw[0] = twLarge[1]
	twInv[0] = twInvLarge[1]
	for m := 1; m <= params.rank/2; m <<= 1 {
		copy(tw[m:2*m], twLarge[2*m:3*m])
		copy(twInv[m:2*m], twInvLarge[2*m:3*m])
	}

	return &pow2AutFixedTransformer{
		params: params,
		mod:    mod,

		tw:     tw,
		twS:    vec.SForm(tw, mod),
		twInv:  twInv,
		twInvS: vec.SForm(twInv, mod),

		rankInv: num.InvMForm(num.Inv(uint64(2*params.rank), mod), mod),

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.rank)
			return &v
		}),
	}
}

func (ntt *pow2AutFixedTransformer) forwardTo(vNTT, v []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(v))

	vBuf[0] = 0
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = wRev[7]
		wOut[1] = wRev[6]
		wOut[2] = wRev[5]
		wOut[3] = wRev[4]

		wOut[4] = wRev[3]
		wOut[5] = wRev[2]
		wOut[6] = wRev[1]
		wOut[7] = wRev[0]
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = v[ntt.params.rank-i]
	}

	tw0Neg, tw0NegS := ntt.mod.Value()-ntt.tw[0], -ntt.twS[0]-1

	copy(vNTT, v)
	vec.SMulAddLazyTo(vNTT, vBuf, tw0Neg, tw0NegS, ntt.mod)
	fwdNTTInPlacePow2(vNTT, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(vNTT, vNTT, ntt.mod)
}

func (ntt *pow2AutFixedTransformer) inverseTo(v, vNTT []uint64) {
	copy(v, vNTT)
	invNTTInPlacePow2(v, ntt.twInv, ntt.twInvS, ntt.mod.Value())

	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(v))

	vBuf[0] = 0
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = wRev[7]
		wOut[1] = wRev[6]
		wOut[2] = wRev[5]
		wOut[3] = wRev[4]

		wOut[4] = wRev[3]
		wOut[5] = wRev[2]
		wOut[6] = wRev[1]
		wOut[7] = wRev[0]
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = v[ntt.params.rank-i]
	}

	vec.SMulAddLazyTo(v, vBuf, ntt.tw[0], ntt.twS[0], ntt.mod)
	v[0] <<= 1
	vec.MulTo(v, v, ntt.rankInv, ntt.mod)
}

// primeAutFixedTransformer is a transformer for [typePrimeAutFixed].
type primeAutFixedTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT *pow235CyclicTransformer

	// isPow2 is true if the rank is a power of two.
	isPow2 bool

	// root is the sum of powers of primitive root.
	root []uint64
	// rootS is the Shoup form of root.
	rootS []uint64
	// rootInv is the sum of powers of inverse primitive root minus fold.
	rootInv []uint64
	// rootInvS is the Shoup form of rootInv.
	rootInvS []uint64

	pool *pool.Pool[*[]uint64]
}

// newPrimeAutFixedTransformer creates a new [primeAutFixedTransformer].
func newPrimeAutFixedTransformer(params RingParameters, mod *num.Modulus) *primeAutFixedTransformer {
	fold := (params.cycloIdx - 1) / params.rank

	isPow2 := num.IsPowerOfTwo(params.rank)

	cycloRoot := num.Generators(uint64(params.cycloIdx))[0]
	modRoot := num.NthRoot(uint64(params.cycloIdx), num.Generators(mod.Value()), mod.Value())

	var ambRank int
	if isPow2 {
		ambRank = params.rank
	} else {
		ambRank = num.NextProdPower(2*params.rank-1, []int{2})
	}
	ambNTT := newPow235CyclicTransformer(NewCyclicParameters(ambRank), mod)

	modRootPowSum := make([]uint64, ambRank)
	modRootPowInvSum := make([]uint64, ambRank)

	cycloRootPowRank := num.Exp(cycloRoot, uint64(params.rank), params.cycloIdxMod)
	modRootPow := modRoot
	modRootPowInv := num.Inv(modRoot, mod)

	var modRootPowNext, modRootPowInvNext uint64
	for i := 0; i < fold; i++ {
		modRootPow = num.Exp(modRootPow, cycloRootPowRank, mod)
		modRootPowInv = num.Exp(modRootPowInv, cycloRootPowRank, mod)

		modRootPowNext = modRootPow
		modRootPowInvNext = modRootPowInv

		for j := 0; j < params.rank; j++ {
			modRootPowSum[j] = num.Add(modRootPowSum[j], modRootPowNext, mod)
			modRootPowInvSum[j] = num.Add(modRootPowInvSum[j], modRootPowInvNext, mod)

			modRootPowNext = num.Exp(modRootPowNext, cycloRoot, mod)
			modRootPowInvNext = num.Exp(modRootPowInvNext, cycloRoot, mod)
		}
	}

	vec.SubTo(modRootPowInvSum[:params.rank], modRootPowInvSum[:params.rank], uint64(fold), mod)

	rootScale := num.MForm(num.Inv(uint64(ambRank), mod), mod)
	rootInvScale := num.InvMForm(num.Inv(uint64(params.cycloIdx)*uint64(ambRank), mod), mod)

	vec.MulTo(modRootPowSum, modRootPowSum, rootScale, mod)
	fwdNTTInPlacePow2(modRootPowSum, ambNTT.tw[0], ambNTT.twS[0], mod.Value())
	vec.Reduce4QTo(modRootPowSum, modRootPowSum, mod)

	vec.MulTo(modRootPowInvSum, modRootPowInvSum, rootInvScale, mod)
	fwdNTTInPlacePow2(modRootPowInvSum, ambNTT.tw[0], ambNTT.twS[0], mod.Value())
	vec.Reduce4QTo(modRootPowInvSum, modRootPowInvSum, mod)

	return &primeAutFixedTransformer{
		params: params,
		mod:    mod,

		ambNTT: ambNTT,

		isPow2: isPow2,

		root:     modRootPowSum,
		rootS:    vec.SForm(modRootPowSum, mod),
		rootInv:  modRootPowInvSum,
		rootInvS: vec.SForm(modRootPowInvSum, mod),

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, ambRank)
			return &v
		}),
	}
}

func (ntt *primeAutFixedTransformer) forwardTo(vNTT, v []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(v))

	vBuf[0] = v[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = wRev[7]
		wOut[1] = wRev[6]
		wOut[2] = wRev[5]
		wOut[3] = wRev[4]

		wOut[4] = wRev[3]
		wOut[5] = wRev[2]
		wOut[6] = wRev[1]
		wOut[7] = wRev[0]
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = v[ntt.params.rank-i]
	}

	clear(vBuf[ntt.params.rank:])

	fwdNTTInPlacePow2(vBuf, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())
	vec.SMulLazyTo(vBuf, vBuf, ntt.root, ntt.rootS, ntt.mod)
	invNTTInPlacePow2(vBuf, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())

	if ntt.isPow2 {
		vec.Reduce2QTo(vNTT, vBuf, ntt.mod)
	} else {
		vec.AddTo(vNTT, vBuf[:ntt.params.rank], vBuf[ntt.params.rank:2*ntt.params.rank], nil)
		vec.Reduce4QTo(vNTT, vNTT, ntt.mod)
	}
}

func (ntt *primeAutFixedTransformer) inverseTo(v, vNTT []uint64) {
	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(vNTT))

	vBuf[0] = vNTT[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = wRev[7]
		wOut[1] = wRev[6]
		wOut[2] = wRev[5]
		wOut[3] = wRev[4]

		wOut[4] = wRev[3]
		wOut[5] = wRev[2]
		wOut[6] = wRev[1]
		wOut[7] = wRev[0]
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = vNTT[ntt.params.rank-i]
	}

	clear(vBuf[ntt.params.rank:])

	fwdNTTInPlacePow2(vBuf, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())
	vec.SMulLazyTo(vBuf, vBuf, ntt.rootInv, ntt.rootInvS, ntt.mod)
	invNTTInPlacePow2(vBuf, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())

	if ntt.isPow2 {
		vec.Reduce2QTo(v, vBuf, ntt.mod)
	} else {
		vec.AddTo(v, vBuf[:ntt.params.rank], vBuf[ntt.params.rank:2*ntt.params.rank], nil)
		vec.Reduce4QTo(v, v, ntt.mod)
	}
}
