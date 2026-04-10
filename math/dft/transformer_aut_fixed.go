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
	// twInv is the twiddle factor for InvNTT.
	twInv []uint64
	// twS is the Shoup form of tw.
	twS []uint64
	// twInvS is the Shoup form of twInv.
	twInvS []uint64

	// rankInv is the modular inverse of the rank.
	rankInv uint64

	pool *pool.Pool[*[]uint64]
}

// newPow2AutFixedTransformer creates a new [pow2AutFixedTransformer].
func newPow2AutFixedTransformer(params RingParameters, mod *num.Modulus) *pow2AutFixedTransformer {
	root := num.Generators(mod)

	twLarge := make([]uint64, 2*params.rank)
	twInvLarge := make([]uint64, 2*params.rank)
	twLarge[0], twLarge[1] = 1, num.NthRoot(params.cycloOrd, root, mod)
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

	twS := make([]uint64, params.rank)
	twInvS := make([]uint64, params.rank)
	for i := 0; i < params.rank; i++ {
		twS[i] = num.SForm(tw[i], mod)
		twInvS[i] = num.SForm(twInv[i], mod)
	}

	return &pow2AutFixedTransformer{
		params: params,
		mod:    mod,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: num.InvMForm(num.Inv(uint64(2*params.rank), mod), mod),

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, params.rank)
			return &v
		}),
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *pow2AutFixedTransformer) ForwardTo(vNTT, v []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	tw0Neg, tw0NegS := ntt.mod.Value()-ntt.tw[0], -ntt.twS[0]-1

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(v))

	vBuf[0] = v[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = w[0] + num.SMul(wRev[7], tw0Neg, tw0NegS, ntt.mod)
		wOut[1] = w[1] + num.SMul(wRev[6], tw0Neg, tw0NegS, ntt.mod)
		wOut[2] = w[2] + num.SMul(wRev[5], tw0Neg, tw0NegS, ntt.mod)
		wOut[3] = w[3] + num.SMul(wRev[4], tw0Neg, tw0NegS, ntt.mod)

		wOut[4] = w[4] + num.SMul(wRev[3], tw0Neg, tw0NegS, ntt.mod)
		wOut[5] = w[5] + num.SMul(wRev[2], tw0Neg, tw0NegS, ntt.mod)
		wOut[6] = w[6] + num.SMul(wRev[1], tw0Neg, tw0NegS, ntt.mod)
		wOut[7] = w[7] + num.SMul(wRev[0], tw0Neg, tw0NegS, ntt.mod)
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = v[i] + num.SMul(v[ntt.params.rank-i], tw0Neg, tw0NegS, ntt.mod)
	}

	nttInPlacePow2(vBuf, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(vNTT, vBuf, ntt.mod)
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *pow2AutFixedTransformer) InverseTo(v, vNTT []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	copy(v, vNTT)
	inttInPlacePow2(v, ntt.twInv, ntt.twInvS, ntt.mod.Value())

	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(v))

	vBuf[0] = v[0] << 1
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		w := (*[8]uint64)(unsafe.Add(r, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = w[0] + num.SMul(wRev[7], ntt.tw[0], ntt.twS[0], ntt.mod)
		wOut[1] = w[1] + num.SMul(wRev[6], ntt.tw[0], ntt.twS[0], ntt.mod)
		wOut[2] = w[2] + num.SMul(wRev[5], ntt.tw[0], ntt.twS[0], ntt.mod)
		wOut[3] = w[3] + num.SMul(wRev[4], ntt.tw[0], ntt.twS[0], ntt.mod)

		wOut[4] = w[4] + num.SMul(wRev[3], ntt.tw[0], ntt.twS[0], ntt.mod)
		wOut[5] = w[5] + num.SMul(wRev[2], ntt.tw[0], ntt.twS[0], ntt.mod)
		wOut[6] = w[6] + num.SMul(wRev[1], ntt.tw[0], ntt.twS[0], ntt.mod)
		wOut[7] = w[7] + num.SMul(wRev[0], ntt.tw[0], ntt.twS[0], ntt.mod)
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = v[i] + num.SMul(v[ntt.params.rank-i], ntt.tw[0], ntt.twS[0], ntt.mod)
	}

	vec.MulScalarTo(v, vBuf, ntt.rankInv, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *pow2AutFixedTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *pow2AutFixedTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

// primeAutFixedTransformer is a transformer for prime order autfixed ring.
type primeAutFixedTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT *pow235CyclicTransformer

	// isPow2 is true of the rank is power-of-two.
	isPow2 bool

	// root is the sum of powers of primitive root.
	// Pre-transformed for a fast convolution.
	root []uint64
	// rootInv is the sum of powers of inverse primitive root.
	// Pre-transformed for a fast convolution.
	rootInv []uint64

	// fold is cyclotomic order divided by rank.
	fold uint64
	// ambRankInvM is the modular inverse of the rank of the ambient NTT in Montgomery form.
	ambRankInvM uint64
	// cycloOrdInv is the modular inverse of the cyclotomic order.
	cycloOrdInv uint64

	pool *pool.Pool[*[]uint64]
}

// newPrimeAutFixedTransformer creates a new [primeAutFixedTransformer].
func newPrimeAutFixedTransformer(params RingParameters, mod *num.Modulus) *primeAutFixedTransformer {
	fold := int((params.cycloOrd - 1) / params.rank)

	isPow2 := num.IsPowerOfTwo(params.rank)

	cycloOrdMod := num.NewModulus(params.cycloOrd)
	cycloRoot := num.Generators(cycloOrdMod)[0]
	modRoot := num.NthRoot(params.cycloOrd, num.Generators(mod), mod)

	var ambRank int
	if isPow2 {
		ambRank = params.rank
	} else {
		ambRank = num.NextProdPower(2*params.rank-1, []int{2})
	}
	ambNTT := newCyclicPow235Transformer(NewCyclicParameters(ambRank), mod)
	ambRankInv := num.MForm(num.Inv(uint64(ambRank), mod), mod)

	modRootPowSum := make([]uint64, ambRank)
	modRootPowInvSum := make([]uint64, ambRank)

	cycloRootPowRank := num.Exp(cycloRoot, uint64(params.rank), cycloOrdMod)
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

	ambNTT.ForwardTo(modRootPowSum, modRootPowSum)
	ambNTT.ForwardTo(modRootPowInvSum, modRootPowInvSum)

	return &primeAutFixedTransformer{
		params: params,
		mod:    mod,

		ambNTT: ambNTT,

		isPow2: isPow2,

		root:    modRootPowSum,
		rootInv: modRootPowInvSum,

		fold:        uint64(fold),
		ambRankInvM: ambRankInv,
		cycloOrdInv: num.Inv(uint64(params.cycloOrd), mod),

		pool: pool.NewPool(func() *[]uint64 {
			v := make([]uint64, ambRank)
			return &v
		}),
	}
}

// ForwardTo transforms the uint64 vector to NTT form.
func (ntt *primeAutFixedTransformer) ForwardTo(vNTT, v []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

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

	nttInPlacePow2(vBuf, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())
	vec.MFormTo(vBuf, vBuf, ntt.mod)

	vec.MMulLazyTo(vBuf, vBuf, ntt.root, ntt.mod)

	inttInPlacePow2(vBuf, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())
	vec.MMulScalarTo(vBuf, vBuf, ntt.ambRankInvM, ntt.mod)

	if ntt.isPow2 {
		copy(vNTT, vBuf)
	} else {
		vec.AddTo(vNTT, vBuf[:ntt.params.rank], vBuf[ntt.params.rank:2*ntt.params.rank], ntt.mod)
	}
}

// InverseTo transforms the uint64 vector to Standard form.
func (ntt *primeAutFixedTransformer) InverseTo(v, vNTT []uint64) {
	checkLength(ntt.params.rank, len(vNTT), len(v))

	vBufPtr := ntt.pool.Get()
	vBuf := *vBufPtr
	defer ntt.pool.Put(vBufPtr)

	M := ((ntt.params.rank - 1) >> 3) << 3
	L := unsafe.Sizeof(uint64(0))

	rOut := unsafe.Pointer(unsafe.SliceData(vBuf))
	r := unsafe.Pointer(unsafe.SliceData(vNTT))

	sumFold := vNTT[0]
	vBuf[0] = vNTT[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Add(rOut, uintptr(1+i)*L))
		wRev := (*[8]uint64)(unsafe.Add(r, uintptr(ntt.params.rank-(i+8))*L))

		wOut[0] = wRev[7]
		wOut[1] = wRev[6]
		wOut[2] = wRev[5]
		wOut[3] = wRev[4]

		sumFold = num.Add(sumFold, wOut[0], ntt.mod)
		sumFold = num.Add(sumFold, wOut[1], ntt.mod)
		sumFold = num.Add(sumFold, wOut[2], ntt.mod)
		sumFold = num.Add(sumFold, wOut[3], ntt.mod)

		wOut[4] = wRev[3]
		wOut[5] = wRev[2]
		wOut[6] = wRev[1]
		wOut[7] = wRev[0]

		sumFold = num.Add(sumFold, wOut[4], ntt.mod)
		sumFold = num.Add(sumFold, wOut[5], ntt.mod)
		sumFold = num.Add(sumFold, wOut[6], ntt.mod)
		sumFold = num.Add(sumFold, wOut[7], ntt.mod)
	}

	for i := M + 1; i < ntt.params.rank; i++ {
		vBuf[i] = vNTT[ntt.params.rank-i]
		sumFold = num.Add(sumFold, vBuf[i], ntt.mod)
	}

	clear(vBuf[ntt.params.rank:])

	nttInPlacePow2(vBuf, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())

	vec.MMulLazyTo(vBuf, vBuf, ntt.rootInv, ntt.mod)

	inttInPlacePow2(vBuf, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())
	vec.MulScalarLazyTo(vBuf, vBuf, ntt.ambNTT.rankInv, ntt.mod)

	if !ntt.isPow2 {
		vec.AddTo(vBuf[:ntt.params.rank-1], vBuf[:ntt.params.rank-1], vBuf[ntt.params.rank:2*ntt.params.rank-1], ntt.mod)
	}

	sumFold = num.MMul(sumFold, ntt.fold, ntt.mod)
	vec.SubScalarTo(v, vBuf[:ntt.params.rank], sumFold, ntt.mod)
	vec.MulScalarTo(v, v, ntt.cycloOrdInv, ntt.mod)
}

// Params returns the ring parameters.
func (ntt *primeAutFixedTransformer) Params() RingParameters {
	return ntt.params
}

// Modulus returns the modulus used for the transform.
func (ntt *primeAutFixedTransformer) Modulus() *num.Modulus {
	return ntt.mod
}
