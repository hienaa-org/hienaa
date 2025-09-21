package dft

import (
	"unsafe"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// autFixedPow2Transformer is a transformer for power-of-two conjugate invariant ring.
type autFixedPow2Transformer struct {
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

	buf transformerBuffer
}

// newAutFixedPow2Transformer creates a new [autFixedPow2Transformer].
func newAutFixedPow2Transformer(params RingParameters, mod *num.Modulus) *autFixedPow2Transformer {
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

	return &autFixedPow2Transformer{
		params: params,
		mod:    mod,

		tw:     tw,
		twS:    twS,
		twInv:  twInv,
		twInvS: twInvS,

		rankInv: num.InvMForm(num.Inv(uint64(2*params.rank), mod), mod),

		buf: newTransformerBuffer(params.rank),
	}
}

func (ntt *autFixedPow2Transformer) ForwardTo(vNTT, v []uint64) {
	tw0Neg, tw0NegS := ntt.mod.Value()-ntt.tw[0], -ntt.twS[0]-1

	M := ((ntt.params.rank - 1) >> 3) << 3

	ntt.buf.coeffs[0] = v[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&ntt.buf.coeffs[1+i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[1+i]))
		wRev := (*[8]uint64)(unsafe.Pointer(&v[ntt.params.rank-(i+8)]))

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
		ntt.buf.coeffs[i] = v[i] + num.SMul(v[ntt.params.rank-i], tw0Neg, tw0NegS, ntt.mod)
	}

	nttInPlacePow2(ntt.buf.coeffs, ntt.tw, ntt.twS, ntt.mod.Value())
	vec.MFormTo(vNTT, ntt.buf.coeffs, ntt.mod)
}

func (ntt *autFixedPow2Transformer) InverseTo(v, vNTT []uint64) {
	copy(v, vNTT)

	inttInPlacePow2(v, ntt.twInv, ntt.twInvS, ntt.mod.Value())

	M := ((ntt.params.rank - 1) >> 3) << 3

	ntt.buf.coeffs[0] = v[0] << 1
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&ntt.buf.coeffs[1+i]))
		w := (*[8]uint64)(unsafe.Pointer(&v[1+i]))
		wRev := (*[8]uint64)(unsafe.Pointer(&v[ntt.params.rank-(i+8)]))

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
		ntt.buf.coeffs[i] = v[i] + num.SMul(v[ntt.params.rank-i], ntt.tw[0], ntt.twS[0], ntt.mod)
	}

	vec.ScalarMulTo(v, ntt.buf.coeffs, ntt.rankInv, ntt.mod)
}

func (ntt *autFixedPow2Transformer) Params() RingParameters {
	return ntt.params
}

func (ntt *autFixedPow2Transformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *autFixedPow2Transformer) SafeCopy() Transformer {
	return &autFixedPow2Transformer{
		params: ntt.params,
		mod:    ntt.mod,

		tw:     ntt.tw,
		twS:    ntt.twS,
		twInv:  ntt.twInv,
		twInvS: ntt.twInvS,

		rankInv: ntt.rankInv,

		buf: newTransformerBuffer(ntt.params.rank),
	}
}

// autFixedPrimeTransformer is a transformer for prime order decomposition ring.
type autFixedPrimeTransformer struct {
	params RingParameters
	mod    *num.Modulus

	ambNTT *cyclicPow235Transformer

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

	buf transformerBuffer
}

// newAutFixedPrimeTransformer creates a new [autFixedPrimeTransformer].
func newAutFixedPrimeTransformer(params RingParameters, mod *num.Modulus) *autFixedPrimeTransformer {
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

	return &autFixedPrimeTransformer{
		params: params,
		mod:    mod,

		ambNTT: ambNTT,

		isPow2: isPow2,

		root:    modRootPowSum,
		rootInv: modRootPowInvSum,

		fold:        uint64(fold),
		ambRankInvM: ambRankInv,
		cycloOrdInv: num.Inv(uint64(params.cycloOrd), mod),

		buf: newTransformerBuffer(ambRank),
	}
}

func (ntt *autFixedPrimeTransformer) ForwardTo(vNTT, v []uint64) {
	M := ((ntt.params.rank - 1) >> 3) << 3

	ntt.buf.coeffs[0] = v[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&ntt.buf.coeffs[1+i]))
		wRev := (*[8]uint64)(unsafe.Pointer(&v[ntt.params.rank-(i+8)]))

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
		ntt.buf.coeffs[i] = v[ntt.params.rank-i]
	}

	clear(ntt.buf.coeffs[ntt.params.rank:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())
	vec.MFormTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.mod)

	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.root, ntt.mod)

	inttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())
	vec.ScalarMMulTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.ambRankInvM, ntt.mod)

	if ntt.isPow2 {
		copy(vNTT, ntt.buf.coeffs)
	} else {
		vec.AddTo(vNTT, ntt.buf.coeffs[:ntt.params.rank], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank], ntt.mod)
	}
}

func (ntt *autFixedPrimeTransformer) InverseTo(v, vNTT []uint64) {
	M := ((ntt.params.rank - 1) >> 3) << 3

	ntt.buf.coeffs[0] = vNTT[0]
	sumFold := vNTT[0]
	for i := 0; i < M; i += 8 {
		wOut := (*[8]uint64)(unsafe.Pointer(&ntt.buf.coeffs[1+i]))
		wRev := (*[8]uint64)(unsafe.Pointer(&vNTT[ntt.params.rank-(i+8)]))

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
		ntt.buf.coeffs[i] = vNTT[ntt.params.rank-i]
		sumFold = num.Add(sumFold, ntt.buf.coeffs[i], ntt.mod)
	}

	clear(ntt.buf.coeffs[ntt.params.rank:])

	nttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.tw[0], ntt.ambNTT.twS[0], ntt.mod.Value())

	vec.MMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.rootInv, ntt.mod)

	inttInPlacePow2(ntt.buf.coeffs, ntt.ambNTT.twInv[0], ntt.ambNTT.twInvS[0], ntt.mod.Value())
	vec.ScalarMulLazyTo(ntt.buf.coeffs, ntt.buf.coeffs, ntt.ambNTT.rankInv, ntt.mod)

	if !ntt.isPow2 {
		vec.AddTo(ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[:ntt.params.rank-1], ntt.buf.coeffs[ntt.params.rank:2*ntt.params.rank-1], ntt.mod)
	}

	sumFold = num.MMul(sumFold, ntt.fold, ntt.mod)
	vec.ScalarSubTo(v, ntt.buf.coeffs[:ntt.params.rank], sumFold, ntt.mod)
	vec.ScalarMulTo(v, v, ntt.cycloOrdInv, ntt.mod)
}

func (ntt *autFixedPrimeTransformer) Params() RingParameters {
	return ntt.params
}

func (ntt *autFixedPrimeTransformer) Modulus() *num.Modulus {
	return ntt.mod
}

func (ntt *autFixedPrimeTransformer) SafeCopy() Transformer {
	return &autFixedPrimeTransformer{
		params: ntt.params,
		mod:    ntt.mod,

		ambNTT: ntt.ambNTT,

		isPow2: ntt.isPow2,

		root:    ntt.root,
		rootInv: ntt.rootInv,

		fold:        ntt.fold,
		ambRankInvM: ntt.ambRankInvM,
		cycloOrdInv: ntt.cycloOrdInv,

		buf: newTransformerBuffer(ntt.ambNTT.params.rank),
	}
}
