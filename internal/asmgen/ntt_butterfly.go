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

		t := VMM(avxType)

		VMULPD(v, wF, t)
		VFMSUB213PD(t, wF, v)
		lo := v

		quo := VMM(avxType)
		VMULPD(t, qfInv, quo)

		switch avxType {
		case TypeAVX2:
			VROUNDPD(Imm(1), quo, quo)
		case TypeAVX512:
			VRNDSCALEPD(Imm(1), quo, quo)
		}

		VFNMADD231PD(quo, qf, t)
		VADDPD(t, lo, t)

		switch avxType {
		case TypeAVX2:
			tQ := VMM(avxType)
			VADDPD(qf, t, tQ)
			VBLENDVPD(t, tQ, t, t)

			VSUBPD(t, u, v)
			VADDPD(t, u, u)

			uQ, vQ := VMM(avxType), VMM(avxType)
			VSUBPD(qf, u, uQ)
			VBLENDVPD(uQ, u, uQ, u)
			VADDPD(qf, v, vQ)
			VBLENDVPD(v, vQ, v, v)

		case TypeAVX512:
			tQ := K()
			VCMPPD(Imm(CMP_LT_OQ), zero, t, tQ)
			VADDPD(qf, t, tQ, t)

			VSUBPD(t, u, v)
			VADDPD(t, u, u)

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

func InvButterflyAVX(avxType AVXType, u, v, w, wS, wF, negQ, twoQ, qf, qfInv, mask52, zero reg.VecVirtual) {
	switch avxType {
	case TypeAVX2, TypeAVX512:
		const CMP_LT_OQ = 0x11
		const CMP_GE_OQ = 0x1d

		t := VMM(avxType)
		VSUBPD(v, u, t)
		VADDPD(v, u, u)
		v, t = t, v

		switch avxType {
		case TypeAVX2:
			uQ, vQ := VMM(avxType), VMM(avxType)
			VSUBPD(qf, u, uQ)
			VBLENDVPD(uQ, u, uQ, u)
			VADDPD(qf, v, vQ)
			VBLENDVPD(v, vQ, v, v)
		case TypeAVX512:
			uQ, vQ := K(), K()
			VCMPPD(Imm(CMP_GE_OQ), qf, u, uQ)
			VSUBPD(qf, u, uQ, u)
			VCMPPD(Imm(CMP_LT_OQ), zero, v, vQ)
			VADDPD(qf, v, vQ, v)
		}

		VMULPD(v, wF, t)
		VFMSUB213PD(t, wF, v)
		lo := v

		quo := VMM(avxType)
		VMULPD(t, qfInv, quo)

		switch avxType {
		case TypeAVX2:
			VROUNDPD(Imm(1), quo, quo)
		case TypeAVX512:
			VRNDSCALEPD(Imm(1), quo, quo)
		}

		VFNMADD231PD(quo, qf, t)
		VADDPD(t, lo, t)

		switch avxType {
		case TypeAVX2:
			tQ := VMM(avxType)
			VADDPD(qf, t, tQ)
			VBLENDVPD(t, tQ, t, t)
		case TypeAVX512:
			tQ := K()
			VCMPPD(Imm(CMP_LT_OQ), zero, t, tQ)
			VADDPD(qf, t, tQ, t)
		}

	case TypeAVX512IFMA:
		VPADDQ(v, u, u)
		VPADDQ(v, v, v)
		VPSUBQ(v, u, v)
		VPADDQ(twoQ, v, v)

		uQ := ZMM()
		VPSUBQ(twoQ, u, uQ)
		VPMINUQ(uQ, u, u)

		quo, t := ZMM(), ZMM()
		VPXORQ(quo, quo, quo)
		VPXORQ(t, t, t)

		VPMADD52HUQ(v, wS, quo)
		VPMADD52LUQ(v, w, t)
		VPMADD52LUQ(quo, negQ, t)
		VPANDQ(t, mask52, v)
	}
}
