package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func FwdButterflyAVX(avxType AVXType, u, v, w, wS, wF, negQ, twoQ, qf, qfInv, mask52, zero reg.VecVirtual) {
	switch avxType {
	case TypeAVX2, TypeAVX512:
		const CMP_LT_OQ = 0x11
		const CMP_GE_OQ = 0x1d

		vw := VMM(avxType)

		VMULPD(v, wF, vw)
		VFMSUB213PD(vw, wF, v)
		lo := v

		quo := VMM(avxType)
		VMULPD(vw, qfInv, quo)

		switch avxType {
		case TypeAVX2:
			VROUNDPD(Imm(1), quo, quo)
		case TypeAVX512:
			VRNDSCALEPD(Imm(1), quo, quo)
		}

		VFNMADD231PD(quo, qf, vw)
		VADDPD(vw, lo, vw)

		switch avxType {
		case TypeAVX2:
			vwQ := VMM(avxType)
			VADDPD(qf, vw, vwQ)
			VBLENDVPD(vw, vwQ, vw, vw)

			VSUBPD(vw, u, v)
			VADDPD(vw, u, u)

			uQ, vQ := VMM(avxType), VMM(avxType)
			VSUBPD(qf, u, uQ)
			VBLENDVPD(uQ, u, uQ, u)
			VADDPD(qf, v, vQ)
			VBLENDVPD(v, vQ, v, v)

		case TypeAVX512:
			vwQ := K()
			VCMPPD(Imm(CMP_LT_OQ), zero, vw, vwQ)
			VADDPD(qf, vw, vwQ, vw)

			VSUBPD(vw, u, v)
			VADDPD(vw, u, u)

			uQ, vQ := K(), K()
			VCMPPD(Imm(CMP_GE_OQ), qf, u, uQ)
			VSUBPD(qf, u, uQ, u)
			VCMPPD(Imm(CMP_LT_OQ), zero, v, vQ)
			VADDPD(qf, v, vQ, v)
		}
	case TypeAVX512IFMA:
		uQ := ZMM()
		VPSUBQ(twoQ, u, uQ)
		VPMINUQ(uQ, u, u)

		quo, t := ZMM(), ZMM()
		VPXORQ(quo, quo, quo)
		VPXORQ(t, t, t)

		VPMADD52HUQ(v, wS, quo)
		VPMADD52LUQ(v, w, t)
		VPMADD52LUQ(quo, negQ, t)
		VPANDQ(t, mask52, t)

		VPSUBQ(t, u, v)
		VPADDQ(twoQ, v, v)
		VPADDQ(t, u, u)
	}
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
