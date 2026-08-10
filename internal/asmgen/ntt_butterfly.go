package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func FwdButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	vHi := ZMM()
	VPSRLQ(Imm(32), v, vHi)

	uQ := ZMM()
	VPSUBQ(twoQ, u, uQ)
	VPMINUQ(uQ, u, u)

	quo, t0, t1 := ZMM(), ZMM(), ZMM()
	Mul64HiApproxAVX512(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, t0)
	VPMULLQ(quo, q, t1)
	VPSUBQ(t1, t0, t0)

	t0Q := ZMM()
	VPSUBQ(twoQ, t0, t0Q)
	VPMINUQ(t0Q, t0, t0)

	VPSUBQ(t0, u, v)
	VPADDQ(twoQ, v, v)
	VPADDQ(t0, u, u)
}

func InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	uSubv := ZMM()
	VPSUBQ(v, u, uSubv)
	VPADDQ(v, u, u)
	VPADDQ(twoQ, uSubv, v)

	uQ := ZMM()
	VPSUBQ(twoQ, u, uQ)
	VPMINUQ(uQ, u, u)

	vHi := ZMM()
	VPSRLQ(Imm(32), v, vHi)

	quo := ZMM()
	Mul64HiApproxAVX512(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, v)
	VPMULLQ(quo, q, quo)
	VPSUBQ(quo, v, v)

	vQ := ZMM()
	VPSUBQ(twoQ, v, vQ)
	VPMINUQ(vQ, v, v)
}
