package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func FwdButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	vHi := ZMM()
	VPSRLQ(Imm(32), v, vHi)

	uSubQ := ZMM()
	VPSUBQ(twoQ, u, uSubQ)
	VPMINUQ(uSubQ, u, u)

	quo, t0, t1 := ZMM(), ZMM(), ZMM()
	Mul64HiAVX512(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, t0)
	VPMULLQ(quo, q, t1)
	VPSUBQ(t1, t0, t0)

	VPSUBQ(t0, u, v)
	VPADDQ(twoQ, v, v)
	VPADDQ(t0, u, u)
}

func InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	VPADDQ(v, u, u)
	VPADDQ(v, v, v)
	VPSUBQ(v, u, v)
	VPADDQ(twoQ, v, v)

	uSubQ := ZMM()
	VPSUBQ(twoQ, u, uSubQ)
	VPMINUQ(uSubQ, u, u)

	vHi := ZMM()
	VPSRLQ(Imm(32), v, vHi)

	quo := ZMM()
	Mul64HiAVX512(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, v)
	VPMULLQ(quo, q, quo)
	VPSUBQ(quo, v, v)
}
