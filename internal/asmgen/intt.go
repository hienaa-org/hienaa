package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
)

func InvNTTInPlacePow2UnrollAVX2() {
	TEXT("inttInPlacePow2UnrollAVX2", NOSPLIT, "func(coeffs, twInv, twInvS []uint64, q uint64)")
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	maskSign := YMM()
	VPSLLQ(Imm(63), allOne, maskSign)

	maskLo, maskHi := YMM(), YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)
	VPSLLQ(Imm(32), maskLo, maskHi)

	coeffs := Load(Param("coeffs").Base(), GP64())
	twInv := Load(Param("twInv").Base(), GP64())
	twInvS := Load(Param("twInvS").Base(), GP64())
	N := Load(Param("coeffs").Len(), GP64())

	wIdx := GP64()

	q := YMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)

	q64 := Load(Param("q"), GP64())
	twoQ64 := GP64()
	MOVQ(q64, twoQ64)
	ADDQ(q64, twoQ64)

	qSwap := YMM()
	ShuffleForLoAVX2(q, qSwap)

	w64, wS64 := GP64(), GP64()
	u64, v64 := GP64(), GP64()

	MOVQ(N, wIdx)
	SHRQ(Imm(1), wIdx)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("t_1_loop_end"))

	Label("t_1_loop_body")

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8})

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8})

	ADDQ(Imm(8), i)

	Label("t_1_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_1_loop_body"))

	MOVQ(N, wIdx)
	SHRQ(Imm(2), wIdx)

	XORQ(i, i)
	JMP(LabelRef("t_2_loop_end"))

	Label("t_2_loop_body")

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	ADDQ(Imm(4), i)

	Label("t_2_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_2_loop_body"))

	t, m := GP64(), GP64()
	MOVQ(U64(4), t)
	MOVQ(N, m)
	SHRQ(Imm(3), m)
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

	w, wS := YMM(), YMM()
	VPBROADCASTQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	wSwap, wSHi := YMM(), YMM()
	ShuffleForLoAVX2(w, wSwap)
	VPSRLQ(Imm(32), wS, wSHi)

	j, jt := GP64(), GP64()
	MOVQ(j1, j)
	MOVQ(j, jt)
	ADDQ(t, jt)
	JMP(LabelRef("j_loop_end"))
	Label("j_loop_body")

	u, v := YMM(), YMM()
	VMOVDQU(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	InvButterflyAVX2(u, v, w, wSwap, wS, wSHi, q, qSwap, maskLo, maskHi, maskSign, allOne)

	VMOVDQU(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(4), j)
	ADDQ(Imm(4), jt)

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

	ShuffleForLoAVX2(w, wSwap)
	VPSRLQ(Imm(32), wS, wSHi)

	XORQ(j, j)
	MOVQ(j, jt)
	ADDQ(t, jt)

	JMP(LabelRef("last_loop_end"))
	Label("last_loop_body")

	VMOVDQU(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	InvButterflyAVX2(u, v, w, wSwap, wS, wSHi, q, qSwap, maskLo, maskHi, maskSign, allOne)

	VMOVDQU(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(4), j)
	ADDQ(Imm(4), jt)

	Label("last_loop_end")
	CMPQ(j, NN)
	JL(LabelRef("last_loop_body"))

	RET()
}

func InvNTTInPlacePow2UnrollAVX512() {
	TEXT("inttInPlacePow2UnrollAVX512", NOSPLIT, "func(coeffs, twInv, twInvS []uint64, q uint64)")
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	maskLo256 := YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo256)

	coeffs := Load(Param("coeffs").Base(), GP64())
	twInv := Load(Param("twInv").Base(), GP64())
	twInvS := Load(Param("twInvS").Base(), GP64())
	N := Load(Param("coeffs").Len(), GP64())

	wIdx := GP64()

	q, twoQ := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPADDQ(q, q, twoQ)

	q256, twoQ256 := YMM(), YMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q256)
	VPADDQ(q256, q256, twoQ256)

	q64 := Load(Param("q"), GP64())
	twoQ64 := GP64()
	MOVQ(q64, twoQ64)
	ADDQ(q64, twoQ64)

	w64, wS64 := GP64(), GP64()
	u64, v64 := GP64(), GP64()

	MOVQ(N, wIdx)
	SHRQ(Imm(1), wIdx)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("t_1_loop_end"))

	Label("t_1_loop_body")

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8})

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8})

	ADDQ(Imm(8), i)

	Label("t_1_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_1_loop_body"))

	MOVQ(N, wIdx)
	SHRQ(Imm(2), wIdx)

	XORQ(i, i)
	JMP(LabelRef("t_2_loop_end"))

	Label("t_2_loop_body")

	MOVQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	InvButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	ADDQ(Imm(4), i)

	Label("t_2_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_2_loop_body"))

	MOVQ(N, wIdx)
	SHRQ(Imm(3), wIdx)

	XORQ(i, i)
	JMP(LabelRef("t_4_loop_end"))
	Label("t_4_loop_body")

	w256, wS256 := YMM(), YMM()
	VPBROADCASTQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w256)
	VPBROADCASTQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS256)
	INCQ(wIdx)

	wSHi256 := YMM()
	VPSRLQ(Imm(32), wS256, wSHi256)

	u256, v256 := YMM(), YMM()
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u256)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8}, v256)

	InvButterflyAVX512YMM(u256, v256, w256, wS256, wSHi256, q256, twoQ256, maskLo256)

	VMOVDQU64(u256, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v256, Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8})

	ADDQ(Imm(8), i)

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

	w, wS := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: twInv, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twInvS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	wSHi := ZMM()
	VPSRLQ(Imm(32), wS, wSHi)

	j, jt := GP64(), GP64()
	MOVQ(j1, j)
	MOVQ(j, jt)
	ADDQ(t, jt)
	JMP(LabelRef("j_loop_end"))
	Label("j_loop_body")

	u, v := ZMM(), ZMM()
	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)

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

	VPSRLQ(Imm(32), wS, wSHi)

	XORQ(j, j)
	MOVQ(j, jt)
	ADDQ(t, jt)

	JMP(LabelRef("last_loop_end"))
	Label("last_loop_body")

	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	InvButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)

	VMOVDQU64(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU64(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(8), j)
	ADDQ(Imm(8), jt)

	Label("last_loop_end")
	CMPQ(j, NN)
	JL(LabelRef("last_loop_body"))

	RET()
}
