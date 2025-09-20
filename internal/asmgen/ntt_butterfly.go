package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func ButterflyAVX2(u, v, w, wSwap, wS, wSHi, q, qSwap, twoQ, maskLo, maskHi, allOne reg.VecVirtual) {
	vHi := YMM()
	VPSRLQ(Imm(32), v, vHi)

	subQ := YMM()
	GreaterOrEqualThanAVX2(u, twoQ, allOne, subQ)
	VPAND(twoQ, subQ, subQ)
	VPSUBQ(subQ, u, u)

	quo, t0, t1 := YMM(), YMM(), YMM()
	Mul64HiAVX2(v, vHi, wS, wSHi, maskLo, quo)
	Mul64LoAVX2(v, w, wSwap, maskHi, t0)
	Mul64LoAVX2(quo, q, qSwap, maskHi, t1)
	VPSUBQ(t1, t0, t0)

	VPSUBQ(t0, u, v)
	VPADDQ(twoQ, v, v)
	VPADDQ(t0, u, u)
}

func ButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	vHi := ZMM()
	VPSRLQ(Imm(32), v, vHi)

	subQ, subQMask := ZMM(), K()
	VPCMPUQ(Imm(0o5), twoQ, u, subQMask)
	VMOVAPD_Z(twoQ, subQMask, subQ)
	VPSUBQ(subQ, u, u)

	quo, t0, t1 := ZMM(), ZMM(), ZMM()
	Mul64HiAVX512(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, t0)
	VPMULLQ(quo, q, t1)
	VPSUBQ(t1, t0, t0)

	VPSUBQ(t0, u, v)
	VPADDQ(twoQ, v, v)
	VPADDQ(t0, u, u)
}

func ButterflyAVX512YMM(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	vHi := YMM()
	VPSRLQ(Imm(32), v, vHi)

	subQ, subQMask := YMM(), K()
	VPCMPUQ(Imm(0o5), twoQ, u, subQMask)
	VMOVAPD_Z(twoQ, subQMask, subQ)
	VPSUBQ(subQ, u, u)

	quo, t0, t1 := YMM(), YMM(), YMM()
	Mul64HiAVX2(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, t0)
	VPMULLQ(quo, q, t1)
	VPSUBQ(t1, t0, t0)

	VPSUBQ(t0, u, v)
	VPADDQ(twoQ, v, v)
	VPADDQ(t0, u, u)
}

func ButterflyX86(u, v, w, wS, q, twoQ reg.Register) {
	subQ := GP64()
	MOVQ(u, subQ)
	SUBQ(twoQ, subQ)
	CMPQ(u, twoQ)
	CMOVQCC(subQ, u)

	quo := GP64()
	MOVQ(wS, reg.RDX)
	MULXQ(v, reg.RDX, quo)

	t := GP64()
	MOVQ(v, t)
	IMULQ(w, t)
	IMULQ(q, quo)
	SUBQ(quo, t)

	MOVQ(u, v)
	SUBQ(t, v)
	ADDQ(twoQ, v)
	ADDQ(t, u)
}
