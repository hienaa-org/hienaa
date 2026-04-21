package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func NTTConstants() {
	ConstData("MASK_LO", U64(1<<32-1))
	ConstData("MASK_52", U64(1<<52-1))

	GLOBL("PERM_00112233", RODATA|NOPTR)
	for i, idx := range []uint64{0, 0, 1, 1, 2, 2, 3, 3} {
		DATA(i*8, U64(idx))
	}

	GLOBL("PERM_02461357", RODATA|NOPTR)
	for i, idx := range []uint64{0, 2, 4, 6, 1, 3, 5, 7} {
		DATA(i*8, U64(idx))
	}

	GLOBL("PERM_04152637", RODATA|NOPTR)
	for i, idx := range []uint64{0, 4, 1, 5, 2, 6, 3, 7} {
		DATA(i*8, U64(idx))
	}
}

func FwdNTTInPlacePow2UnrollAVX512(isIFMA bool) {
	if !isIFMA {
		TEXT("fwdNTTInPlacePow2UnrollAVX512", NOSPLIT, "func(coeffs, tw, twS []uint64, q uint64)")
	} else {
		TEXT("fwdNTTInPlacePow2UnrollAVX512IFMA", NOSPLIT, "func(coeffs, tw, twS []uint64, q uint64)")
	}
	Pragma("noescape")

	var maskLo, mask52, zero reg.VecVirtual
	if !isIFMA {
		maskLo = ZMM()
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)
	} else {
		mask52 = ZMM()
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_52"), 0), mask52)
		zero = ZMM()
		VPXORQ(zero, zero, zero)
	}

	coeffs := Load(Param("coeffs").Base(), GP64())
	tw := Load(Param("tw").Base(), GP64())
	twS := Load(Param("twS").Base(), GP64())
	N := Load(Param("coeffs").Len(), GP64())

	wIdx := GP64()
	MOVQ(U64(1), wIdx)

	q, twoQ := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPADDQ(q, q, twoQ)

	w, wS := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	var wSHi reg.VecVirtual
	if !isIFMA {
		wSHi = ZMM()
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	NN := GP64()
	MOVQ(N, NN)
	SHRQ(Imm(1), NN)

	t, m := GP64(), GP64()
	MOVQ(NN, t)

	j, jt := GP64(), GP64()
	XORQ(j, j)
	MOVQ(j, jt)
	ADDQ(t, jt)

	JMP(LabelRef("first_loop_end"))
	Label("first_loop_body")

	u, v := ZMM(), ZMM()
	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	if !isIFMA {
		FwdButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		FwdButterflyAVX512IFMA(u, v, w, wS, q, twoQ, zero, mask52)
	}

	VMOVDQU64(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(8), j)
	ADDQ(Imm(8), jt)

	Label("first_loop_end")
	CMPQ(j, NN)
	JL(LabelRef("first_loop_body"))

	MOVQ(N, NN)
	SHRQ(Imm(4), NN)

	MOVQ(U64(2), m)
	JMP(LabelRef("m_loop_end"))
	Label("m_loop_body")

	SHRQ(Imm(1), t)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("i_loop_end"))
	Label("i_loop_body")

	j1, j2 := GP64(), GP64()
	MOVQ(i, j1)
	IMULQ(t, j1)
	SHLQ(Imm(1), j1)
	MOVQ(j1, j2)
	ADDQ(t, j2)

	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	MOVQ(j1, j)
	MOVQ(j2, jt)
	JMP(LabelRef("j_loop_end"))
	Label("j_loop_body")

	u, v = ZMM(), ZMM()
	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	if !isIFMA {
		FwdButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		FwdButterflyAVX512IFMA(u, v, w, wS, q, twoQ, zero, mask52)
	}

	VMOVDQU64(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(8), j)
	ADDQ(Imm(8), jt)

	Label("j_loop_end")
	CMPQ(j, j2)
	JL(LabelRef("j_loop_body"))

	ADDQ(Imm(1), i)

	Label("i_loop_end")
	CMPQ(i, m)
	JL(LabelRef("i_loop_body"))

	SHLQ(Imm(1), m)

	Label("m_loop_end")
	CMPQ(m, NN)
	JLE(LabelRef("m_loop_body"))

	XORQ(i, i)
	JMP(LabelRef("t_4_loop_end"))
	Label("t_4_loop_body")

	w0, w0S := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w0)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, w0S)
	INCQ(wIdx)

	w1, w1S := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w1)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, w1S)
	INCQ(wIdx)

	VSHUFI64X2(Imm(0b10_10_00_00), w1, w0, w)
	VSHUFI64X2(Imm(0b10_10_00_00), w1S, w0S, wS)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8}, v)

	uu, vv := ZMM(), ZMM()
	VSHUFI64X2(Imm(0b01_00_01_00), v, u, uu)
	VSHUFI64X2(Imm(0b11_10_11_10), v, u, vv)

	if !isIFMA {
		FwdButterflyAVX512(uu, vv, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		FwdButterflyAVX512IFMA(uu, vv, w, wS, q, twoQ, zero, mask52)
	}

	VSHUFI64X2(Imm(0b01_00_01_00), vv, uu, u)
	VSHUFI64X2(Imm(0b11_10_11_10), vv, uu, v)

	VMOVDQU64(u, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8})

	ADDQ(Imm(16), i)

	Label("t_4_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_4_loop_body"))

	PERM_00112233 := ZMM()
	VMOVDQU64(NewDataAddr(NewStaticSymbol("PERM_00112233"), 0), PERM_00112233)

	XORQ(i, i)
	JMP(LabelRef("t_2_loop_end"))
	Label("t_2_loop_body")

	VMOVDQU64(Mem{Base: tw, Index: wIdx, Scale: 8}, w.AsY())
	VMOVDQU64(Mem{Base: twS, Index: wIdx, Scale: 8}, wS.AsY())
	ADDQ(Imm(4), wIdx)

	VPERMQ(w, PERM_00112233, w)
	VPERMQ(wS, PERM_00112233, wS)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8}, v)

	VSHUFI64X2(Imm(0b10_00_10_00), v, u, uu)
	VSHUFI64X2(Imm(0b11_01_11_01), v, u, vv)

	if !isIFMA {
		FwdButterflyAVX512(uu, vv, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		FwdButterflyAVX512IFMA(uu, vv, w, wS, q, twoQ, zero, mask52)
	}

	VSHUFI64X2(Imm(0b01_00_01_00), vv, uu, u)
	VSHUFI64X2(Imm(0b11_10_11_10), vv, uu, v)
	VSHUFI64X2(Imm(0b11_01_10_00), u, u, u)
	VSHUFI64X2(Imm(0b11_01_10_00), v, v, v)

	VMOVDQU64(u, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8})

	ADDQ(Imm(16), i)

	Label("t_2_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_2_loop_body"))

	PERM_02461357, PERM_04152637 := ZMM(), ZMM()
	VMOVDQU64(NewDataAddr(NewStaticSymbol("PERM_02461357"), 0), PERM_02461357)
	VMOVDQU64(NewDataAddr(NewStaticSymbol("PERM_04152637"), 0), PERM_04152637)

	XORQ(i, i)
	JMP(LabelRef("t_1_loop_end"))
	Label("t_1_loop_body")

	VMOVDQU64(Mem{Base: tw, Index: wIdx, Scale: 8}, w)
	VMOVDQU64(Mem{Base: twS, Index: wIdx, Scale: 8}, wS)
	ADDQ(Imm(8), wIdx)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8}, v)

	VPERMQ(u, PERM_02461357, u)
	VPERMQ(v, PERM_02461357, v)

	VSHUFF64X2(Imm(0b01_00_01_00), v, u, uu)
	VSHUFF64X2(Imm(0b11_10_11_10), v, u, vv)

	if !isIFMA {
		FwdButterflyAVX512(uu, vv, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		FwdButterflyAVX512IFMA(uu, vv, w, wS, q, twoQ, zero, mask52)
	}

	VSHUFF64X2(Imm(0b01_00_01_00), vv, uu, u)
	VSHUFF64X2(Imm(0b11_10_11_10), vv, uu, v)

	VPERMQ(u, PERM_04152637, u)
	VPERMQ(v, PERM_04152637, v)

	VMOVDQU64(u, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8})

	ADDQ(Imm(16), i)

	Label("t_1_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_1_loop_body"))

	RET()
}

func InvNTTInPlacePow2UnrollAVX512(isIFMA bool) {
	if !isIFMA {
		TEXT("invNTTInPlacePow2UnrollAVX512", NOSPLIT, "func(coeffs, twInv, twInvS []uint64, q uint64)")
	} else {
		TEXT("invNTTInPlacePow2UnrollAVX512IFMA", NOSPLIT, "func(coeffs, twInv, twInvS []uint64, q uint64)")
	}
	Pragma("noescape")

	var maskLo, mask52, zero reg.VecVirtual
	if !isIFMA {
		maskLo = ZMM()
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)
	} else {
		mask52 = ZMM()
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_52"), 0), mask52)
		zero = ZMM()
		VPXORQ(zero, zero, zero)
	}

	coeffs := Load(Param("coeffs").Base(), GP64())
	twInv := Load(Param("twInv").Base(), GP64())
	twInvS := Load(Param("twInvS").Base(), GP64())
	N := Load(Param("coeffs").Len(), GP64())

	wIdx := GP64()

	q, twoQ := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPADDQ(q, q, twoQ)

	MOVQ(N, wIdx)
	SHRQ(Imm(1), wIdx)

	PERM_02461357, PERM_04152637 := ZMM(), ZMM()
	VMOVDQU64(NewDataAddr(NewStaticSymbol("PERM_02461357"), 0), PERM_02461357)
	VMOVDQU64(NewDataAddr(NewStaticSymbol("PERM_04152637"), 0), PERM_04152637)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("t_1_loop_end"))
	Label("t_1_loop_body")

	w, wS := ZMM(), ZMM()
	VMOVDQU64(Mem{Base: twInv, Index: wIdx, Scale: 8}, w)
	VMOVDQU64(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS)
	ADDQ(Imm(8), wIdx)

	var wSHi reg.VecVirtual
	if !isIFMA {
		wSHi = ZMM()
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	u, v := ZMM(), ZMM()
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8}, v)

	VPERMQ(u, PERM_02461357, u)
	VPERMQ(v, PERM_02461357, v)

	uu, vv := ZMM(), ZMM()
	VSHUFF64X2(Imm(0b01_00_01_00), v, u, uu)
	VSHUFF64X2(Imm(0b11_10_11_10), v, u, vv)

	if !isIFMA {
		InvButterflyAVX512(uu, vv, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		InvButterflyAVX512IFMA(uu, vv, w, wS, wSHi, q, twoQ, zero, mask52)
	}

	VSHUFF64X2(Imm(0b01_00_01_00), vv, uu, u)
	VSHUFF64X2(Imm(0b11_10_11_10), vv, uu, v)

	VPERMQ(u, PERM_04152637, u)
	VPERMQ(v, PERM_04152637, v)

	VMOVDQU64(u, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8})

	ADDQ(Imm(16), i)

	Label("t_1_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_1_loop_body"))

	PERM_00112233 := ZMM()
	VMOVDQU64(NewDataAddr(NewStaticSymbol("PERM_00112233"), 0), PERM_00112233)

	MOVQ(N, wIdx)
	SHRQ(Imm(2), wIdx)

	XORQ(i, i)
	JMP(LabelRef("t_2_loop_end"))
	Label("t_2_loop_body")

	VMOVDQU64(Mem{Base: twInv, Index: wIdx, Scale: 8}, w.AsY())
	VMOVDQU64(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS.AsY())
	ADDQ(Imm(4), wIdx)

	VPERMQ(w, PERM_00112233, w)
	VPERMQ(wS, PERM_00112233, wS)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8}, v)

	VSHUFI64X2(Imm(0b10_00_10_00), v, u, uu)
	VSHUFI64X2(Imm(0b11_01_11_01), v, u, vv)

	if !isIFMA {
		InvButterflyAVX512(uu, vv, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		InvButterflyAVX512IFMA(uu, vv, w, wS, wSHi, q, twoQ, zero, mask52)
	}

	VSHUFI64X2(Imm(0b01_00_01_00), vv, uu, u)
	VSHUFI64X2(Imm(0b11_10_11_10), vv, uu, v)
	VSHUFI64X2(Imm(0b11_01_10_00), u, u, u)
	VSHUFI64X2(Imm(0b11_01_10_00), v, v, v)

	VMOVDQU64(u, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8})

	ADDQ(Imm(16), i)

	Label("t_2_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_2_loop_body"))

	MOVQ(N, wIdx)
	SHRQ(Imm(3), wIdx)

	XORQ(i, i)
	JMP(LabelRef("t_4_loop_end"))
	Label("t_4_loop_body")

	w0, w0S := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w0)
	VPBROADCASTQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, w0S)
	INCQ(wIdx)

	w1, w1S := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w1)
	VPBROADCASTQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, w1S)
	INCQ(wIdx)

	VSHUFI64X2(Imm(0b10_10_00_00), w1, w0, w)
	VSHUFI64X2(Imm(0b10_10_00_00), w1S, w0S, wS)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8}, v)

	VSHUFI64X2(Imm(0b01_00_01_00), v, u, uu)
	VSHUFI64X2(Imm(0b11_10_11_10), v, u, vv)

	if !isIFMA {
		InvButterflyAVX512(uu, vv, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		InvButterflyAVX512IFMA(uu, vv, w, wS, wSHi, q, twoQ, zero, mask52)
	}

	VSHUFI64X2(Imm(0b01_00_01_00), vv, uu, u)
	VSHUFI64X2(Imm(0b11_10_11_10), vv, uu, v)

	VMOVDQU64(u, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: i, Disp: 8 * 8, Scale: 8})

	ADDQ(Imm(16), i)

	Label("t_4_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_4_loop_body"))

	t, m := GP64(), GP64()
	MOVQ(U64(8), t)
	MOVQ(N, m)
	SHRQ(Imm(4), m)
	JMP(LabelRef("m_loop_end"))
	Label("m_loop_body")

	MOVQ(m, wIdx)

	XORQ(i, i)
	JMP(LabelRef("i_loop_end"))
	Label("i_loop_body")

	j1, j2 := GP64(), GP64()
	MOVQ(i, j1)
	IMULQ(t, j1)
	SHLQ(Imm(1), j1)
	MOVQ(j1, j2)
	ADDQ(t, j2)

	VPBROADCASTQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	j, jt := GP64(), GP64()
	MOVQ(j1, j)
	MOVQ(j, jt)
	ADDQ(t, jt)
	JMP(LabelRef("j_loop_end"))
	Label("j_loop_body")

	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	if !isIFMA {
		InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		InvButterflyAVX512IFMA(u, v, w, wS, wSHi, q, twoQ, zero, mask52)
	}

	VMOVDQU64(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(8), j)
	ADDQ(Imm(8), jt)

	Label("j_loop_end")
	CMPQ(j, j2)
	JL(LabelRef("j_loop_body"))

	ADDQ(Imm(1), i)

	Label("i_loop_end")
	CMPQ(i, m)
	JL(LabelRef("i_loop_body"))

	SHLQ(Imm(1), t)
	SHRQ(Imm(1), m)

	Label("m_loop_end")
	CMPQ(m, Imm(2))
	JGE(LabelRef("m_loop_body"))

	NN := GP64()
	MOVQ(N, NN)
	SHRQ(Imm(1), NN)

	VPBROADCASTQ(Mem{Base: twInv, Disp: 8, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twInvS, Disp: 8, Scale: 8}, wS)

	if !isIFMA {
		VPSRLQ(Imm(32), wS, wSHi)
	} else {
		VPSRLQ(Imm(12), wS, wS)
	}

	XORQ(j, j)
	MOVQ(j, jt)
	ADDQ(t, jt)

	JMP(LabelRef("last_loop_end"))
	Label("last_loop_body")

	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	if !isIFMA {
		InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)
	} else {
		InvButterflyAVX512IFMA(u, v, w, wS, wSHi, q, twoQ, zero, mask52)
	}

	VMOVDQU64(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(8), j)
	ADDQ(Imm(8), jt)

	Label("last_loop_end")
	CMPQ(j, NN)
	JL(LabelRef("last_loop_body"))

	RET()
}
