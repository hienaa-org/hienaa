package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
)

func NTTConstants() {
	ConstData("MASK_LO", U64(1<<32-1))
}

func NTTInPlacePow2UnrollAVX2() {
	TEXT("nttInPlacePow2UnrollAVX2", NOSPLIT, "func(coeffs, tw, twS []uint64, q uint64)")
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	maskLo, maskHi := YMM(), YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)
	VPSLLQ(Imm(32), maskLo, maskHi)

	maskSign := YMM()
	VPSLLQ(Imm(63), allOne, maskSign)

	coeffs := Load(Param("coeffs").Base(), GP64())
	tw := Load(Param("tw").Base(), GP64())
	twS := Load(Param("twS").Base(), GP64())
	N := Load(Param("coeffs").Len(), GP64())

	wIdx := GP64()
	MOVQ(U64(1), wIdx)

	q, twoQ := YMM(), YMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPADDQ(q, q, twoQ)

	q64 := Load(Param("q"), GP64())
	twoQ64 := GP64()
	MOVQ(q64, twoQ64)
	ADDQ(q64, twoQ64)

	qSwap := YMM()
	ShuffleForLoAVX2(q, qSwap)

	w, wS := YMM(), YMM()
	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	wSwap, wSHi := YMM(), YMM()
	ShuffleForLoAVX2(w, wSwap)
	VPSRLQ(Imm(32), wS, wSHi)

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

	u, v := YMM(), YMM()
	VMOVDQU(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	ButterflyAVX2(u, v, w, wSwap, wS, wSHi, q, qSwap, twoQ, maskLo, maskHi, maskSign, allOne)

	VMOVDQU(u, Mem{Base: coeffs, Index: j, Scale: 8})
	VMOVDQU(v, Mem{Base: coeffs, Index: jt, Scale: 8})

	ADDQ(Imm(4), j)
	ADDQ(Imm(4), jt)

	Label("first_loop_end")
	CMPQ(j, NN)
	JL(LabelRef("first_loop_body"))

	MOVQ(N, NN)
	SHRQ(Imm(3), NN)

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

	ShuffleForLoAVX2(w, wSwap)
	VPSRLQ(Imm(32), wS, wSHi)

	MOVQ(j1, j)
	MOVQ(j2, jt)
	JMP(LabelRef("j_loop_end"))
	Label("j_loop_body")

	u, v = YMM(), YMM()
	VMOVDQU(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	ButterflyAVX2(u, v, w, wSwap, wS, wSHi, q, qSwap, twoQ, maskLo, maskHi, maskSign, allOne)

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

	SHLQ(Imm(1), m)

	Label("m_loop_end")
	CMPQ(m, NN)
	JLE(LabelRef("m_loop_body"))

	w64, wS64 := GP64(), GP64()
	u64, v64 := GP64(), GP64()

	XORQ(i, i)
	JMP(LabelRef("t_2_loop_end"))

	Label("t_2_loop_body")

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	ADDQ(Imm(4), i)

	Label("t_2_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_2_loop_body"))

	XORQ(i, i)
	JMP(LabelRef("t_1_loop_end"))

	Label("t_1_loop_body")

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8})

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8})

	ADDQ(Imm(8), i)

	Label("t_1_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_1_loop_body"))

	RET()
}

func NTTInPlacePow2UnrollAVX512() {
	TEXT("nttInPlacePow2UnrollAVX512", NOSPLIT, "func(coeffs, tw, twS []uint64, q uint64)")
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	maskLo256 := YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo256)

	coeffs := Load(Param("coeffs").Base(), GP64())
	tw := Load(Param("tw").Base(), GP64())
	twS := Load(Param("twS").Base(), GP64())
	N := Load(Param("coeffs").Len(), GP64())

	wIdx := GP64()
	MOVQ(U64(1), wIdx)

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

	w, wS := ZMM(), ZMM()
	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS)
	INCQ(wIdx)

	wSHi := ZMM()
	VPSRLQ(Imm(32), wS, wSHi)

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

	ButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)

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

	VPSRLQ(Imm(32), wS, wSHi)

	MOVQ(j1, j)
	MOVQ(j2, jt)
	JMP(LabelRef("j_loop_end"))
	Label("j_loop_body")

	u, v = ZMM(), ZMM()
	VMOVDQU64(Mem{Base: coeffs, Index: j, Scale: 8}, u)
	VMOVDQU64(Mem{Base: coeffs, Index: jt, Scale: 8}, v)

	ButterflyAVX512(u, v, w, wS, wSHi, q, twoQ, maskLo)

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

	w256, wS256 := YMM(), YMM()
	VPBROADCASTQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w256)
	VPBROADCASTQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS256)
	INCQ(wIdx)

	wSHi256 := YMM()
	VPSRLQ(Imm(32), wS256, wSHi256)

	u256, v256 := YMM(), YMM()
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u256)
	VMOVDQU64(Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8}, v256)

	ButterflyAVX512YMM(u256, v256, w256, wS256, wSHi256, q256, twoQ256, maskLo256)

	VMOVDQU64(u256, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	VMOVDQU64(v256, Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8})

	ADDQ(Imm(8), i)

	Label("t_4_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_4_loop_body"))

	w64, wS64 := GP64(), GP64()
	u64, v64 := GP64(), GP64()

	XORQ(i, i)
	JMP(LabelRef("t_2_loop_end"))

	Label("t_2_loop_body")

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	ADDQ(Imm(4), i)

	Label("t_2_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_2_loop_body"))

	XORQ(i, i)
	JMP(LabelRef("t_1_loop_end"))

	Label("t_1_loop_body")

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 0 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 1 * 8, Scale: 8})

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 2 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 3 * 8, Scale: 8})

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 4 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 5 * 8, Scale: 8})

	MOVQ(Mem{Base: tw, Index: wIdx, Scale: 8}, w64)
	MOVQ(Mem{Base: twS, Index: wIdx, Scale: 8}, wS64)
	INCQ(wIdx)

	MOVQ(Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8}, u64)
	MOVQ(Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8}, v64)

	ButterflyX86(u64, v64, w64, wS64, q64, twoQ64)

	MOVQ(u64, Mem{Base: coeffs, Index: i, Disp: 6 * 8, Scale: 8})
	MOVQ(v64, Mem{Base: coeffs, Index: i, Disp: 7 * 8, Scale: 8})

	ADDQ(Imm(8), i)

	Label("t_1_loop_end")
	CMPQ(i, N)
	JL(LabelRef("t_1_loop_body"))

	RET()
}
