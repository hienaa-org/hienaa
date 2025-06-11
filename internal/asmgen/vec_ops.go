package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func VecOpConstants() {
	ConstData("ONE", U64(1))
}

func VecAddToAVX2(isLazy bool) {
	if isLazy {
		TEXT("addLazyToAVX2", NOSPLIT, "func(v0, v1, vOut []uint64)")
	} else {
		TEXT("addToAVX2", NOSPLIT, "func(v0, v1 []uint64, q uint64, vOut []uint64)")
	}
	Pragma("noescape")

	var q reg.Register
	var qv, qvNegOne reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())

		qv = YMM()
		VPBROADCASTQ(NewParamAddr("q", 48), qv)

		one := YMM()
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("ONE"), 0), one)

		qvNegOne = YMM()
		VPSUBQ(one, qv, qvNegOne)
	}

	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	N := Load(Param("vOut").Len(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(2), M)
	SHLQ(Imm(2), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := YMM(), YMM()
	VMOVDQU(Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOVDQU(Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOut := YMM()
	VPADDQ(x1, x0, xOut)

	if !isLazy {
		subQ := YMM()
		GreaterThanAVX2(xOut, qvNegOne, subQ)
		VPAND(qv, subQ, subQ)
		VPSUBQ(subQ, xOut, xOut)
	}

	VMOVDQU(xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(4), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y0, y1 := GP64(), GP64()
	MOVQ(Mem{Base: v0, Index: i, Scale: 8}, y0)
	MOVQ(Mem{Base: v1, Index: i, Scale: 8}, y1)

	ADDQ(y1, y0)

	if !isLazy {
		ySubQ := GP64()
		MOVQ(y0, ySubQ)
		SUBQ(q, ySubQ)
		CMPQ(y0, q)
		CMOVQGT(ySubQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecSubToAVX2(isLazy bool) {
	if isLazy {
		TEXT("subLazyToAVX2", NOSPLIT, "func(v0, v1, vOut []uint64)")
	} else {
		TEXT("subToAVX2", NOSPLIT, "func(v0, v1 []uint64, q uint64, vOut []uint64)")
	}
	Pragma("noescape")

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())

		qv = YMM()
		VPBROADCASTQ(NewParamAddr("q", 48), qv)
	}

	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	N := Load(Param("vOut").Len(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(2), M)
	SHLQ(Imm(2), M)

	zero := YMM()
	VPXOR(zero, zero, zero)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := YMM(), YMM()
	VMOVDQU(Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOVDQU(Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOut := YMM()
	VPSUBQ(x1, x0, xOut)

	if !isLazy {
		subQ := YMM()
		LessThanAVX2(xOut, zero, subQ)
		VPAND(qv, subQ, subQ)
		VPADDQ(subQ, xOut, xOut)
	}

	VMOVDQU(xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(4), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y0, y1 := GP64(), GP64()
	MOVQ(Mem{Base: v0, Index: i, Scale: 8}, y0)
	MOVQ(Mem{Base: v1, Index: i, Scale: 8}, y1)

	SUBQ(y1, y0)

	if !isLazy {
		ySubQ := GP64()
		MOVQ(y0, ySubQ)
		ADDQ(q, ySubQ)
		CMPQ(y0, Imm(0))
		CMOVQLT(ySubQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecBMulToAVX2(isLazy bool) {
	if isLazy {
		TEXT("bMulLazyToAVX2", NOSPLIT, "func(v0, v1 []uint64, q, divHi, divLo uint64, vOut []uint64)")
	} else {
		TEXT("bMulToAVX2", NOSPLIT, "func(v0, v1 []uint64, q, divHi, divLo uint64, vOut []uint64)")
	}
	Pragma("noescape")

	// q := Load(Param("q"), GP64())
	// divHi := Load(Param("divHi"), GP64())
	// divLo := Load(Param("divLo"), GP64())

	q, divHi, divLo := YMM(), YMM(), YMM()
	VPBROADCASTQ(NewParamAddr("q", 48), q)
	VPBROADCASTQ(NewParamAddr("divHi", 48+8), divHi)
	VPBROADCASTQ(NewParamAddr("divLo", 48+16), divLo)

	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	N := Load(Param("vOut").Len(), GP64())

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	// x0, x1 := reg.RAX, GP64()
	// MOVQ(Mem{Base: v0, Index: i, Scale: 8}, x0)
	// MOVQ(Mem{Base: v1, Index: i, Scale: 8}, x1)

	// MULQ(x1)

	// xOutHi, xOutLo, xOut := GP64(), GP64(), GP64()
	// MOVQ(reg.RDX, xOutHi)
	// MOVQ(reg.RAX, xOutLo)

	// if isLazy {
	// 	BModLazy(xOutHi, xOutLo, q, divHi, divLo, xOut)
	// } else {
	// 	BMod(xOutHi, xOutLo, q, divHi, divLo, xOut)
	// }

	// MOVQ(xOut, Mem{Base: vOut, Index: i, Scale: 8})
	// ADDQ(Imm(1), i)

	x0, x1 := YMM(), YMM()
	VMOVDQU(Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOVDQU(Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOutHi, xOutLo, xOut := YMM(), YMM(), YMM()
	Mul64AVX2(x0, x1, xOutHi, xOutLo)
	if isLazy {
		BModLazyAVX2(xOutHi, xOutLo, q, divHi, divLo, xOut)
	} else {
		BModAVX2(xOutHi, xOutLo, q, divHi, divLo, xOut)
	}

	VMOVDQU(xOut, Mem{Base: vOut, Index: i, Scale: 8})
	ADDQ(Imm(4), i)

	Label("loop_end")
	CMPQ(i, N)
	JL(LabelRef("loop_body"))

	RET()
}
