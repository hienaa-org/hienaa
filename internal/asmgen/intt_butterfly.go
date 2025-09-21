package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func InvButterflyAVX2(u, v, w, wSwap, wS, wSHi, q, qSwap, maskLo, maskHi, maskSign, allOne reg.VecVirtual) {
	twoQ := YMM()
	VPADDQ(q, q, twoQ)

	VPADDQ(v, u, u)
	VPADDQ(v, v, v)
	VPSUBQ(v, u, v)
	VPADDQ(twoQ, v, v)

	subQ := YMM()
	GreaterOrEqualThanAVX2(u, twoQ, maskSign, allOne, subQ)
	VPAND(twoQ, subQ, subQ)
	VPSUBQ(subQ, u, u)

	vHi := YMM()
	VPSRLQ(Imm(32), v, vHi)

	quo := YMM()
	Mul64HiAVX2(v, vHi, wS, wSHi, maskLo, quo)
	Mul64LoAVX2(quo, q, qSwap, maskHi, quo)
	Mul64LoAVX2(v, w, wSwap, maskHi, v)
	VPSUBQ(quo, v, v)
}

func InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	VPADDQ(v, u, u)
	VPADDQ(v, v, v)
	VPSUBQ(v, u, v)
	VPADDQ(twoQ, v, v)

	subQ, subQMask := ZMM(), K()
	VPCMPUQ(Imm(0o5), twoQ, u, subQMask)
	VMOVAPD_Z(twoQ, subQMask, subQ)
	VPSUBQ(subQ, u, u)

	vHi := ZMM()
	VPSRLQ(Imm(32), v, vHi)

	quo := ZMM()
	Mul64HiAVX512(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, v)
	VPMULLQ(quo, q, quo)
	VPSUBQ(quo, v, v)
}

func InvButterflyAVX512YMM(u, v, w, wS, wSHi, q, twoQ, maskLo reg.VecVirtual) {
	VPADDQ(v, u, u)
	VPADDQ(v, v, v)
	VPSUBQ(v, u, v)
	VPADDQ(twoQ, v, v)

	subQ, subQMask := YMM(), K()
	VPCMPUQ(Imm(0o5), twoQ, u, subQMask)
	VMOVAPD_Z(twoQ, subQMask, subQ)
	VPSUBQ(subQ, u, u)

	vHi := YMM()
	VPSRLQ(Imm(32), v, vHi)

	quo := YMM()
	Mul64HiAVX2(v, vHi, wS, wSHi, maskLo, quo)
	VPMULLQ(v, w, v)
	VPMULLQ(quo, q, quo)
	VPSUBQ(quo, v, v)
}

func InvButterflyX86(u, v, w, wS, q, twoQ reg.Register) {
	ADDQ(v, u)
	NEGQ(v)
	ADDQ(v, v)
	ADDQ(u, v)
	ADDQ(twoQ, v)

	subQ := GP64()
	MOVQ(u, subQ)
	SUBQ(twoQ, subQ)
	CMPQ(u, twoQ)
	CMOVQCC(subQ, u)

	quo := GP64()
	MOVQ(wS, reg.RDX)
	MULXQ(v, reg.RDX, quo)

	IMULQ(w, v)
	IMULQ(q, quo)
	SUBQ(quo, v)
}
