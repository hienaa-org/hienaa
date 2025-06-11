package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func VecAddToAVX2(isLazy bool) {
	if isLazy {
		TEXT("addLazyToAVX2", NOSPLIT, "func(v0, v1, vOut []uint64)")
	} else {
		TEXT("addToAVX2", NOSPLIT, "func(v0, v1 []uint64, q uint64, vOut []uint64)")
	}
	Pragma("noescape")

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = YMM()
		VPBROADCASTQ(q, qv)
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
		leftOver := YMM()
		GreaterOrEqualThanAVX2(xOut, qv, leftOver)
		VPANDQ(qv, leftOver, leftOver)
		VPSUBQ(leftOver, xOut, xOut)
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
		leftOver := GP64()
		XORQ(leftOver, leftOver)
		CMPQ(y0, q)
		CMOVQGE(q, leftOver)
		SUBQ(leftOver, y0)
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
		leftOver := YMM()
		LessThanAVX2(xOut, zero, leftOver)
		VPANDQ(qv, leftOver, leftOver)
		VPADDQ(leftOver, xOut, xOut)
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
		leftOver := GP64()
		XORQ(leftOver, leftOver)
		CMPQ(y0, Imm(0))
		CMOVQLT(q, leftOver)
		ADDQ(leftOver, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}
