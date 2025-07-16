package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func AddVecToAVX2(isLazy bool) {
	if isLazy {
		TEXT("addLazyVecToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("addVecToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
	}
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = YMM()
		VPBROADCASTQ(NewParamAddr("q", 72), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())

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
		GreaterOrEqualThanAVX2(xOut, qv, allOne, subQ)
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
		subQ := GP64()
		MOVQ(y0, subQ)
		SUBQ(q, subQ)
		CMPQ(y0, q)
		CMOVQGE(subQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func AddVecToAVX512(isLazy bool) {
	if isLazy {
		TEXT("addLazyVecToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("addVecToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
	}
	Pragma("noescape")

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = ZMM()
		VPBROADCASTQ(NewParamAddr("q", 72), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(3), M)
	SHLQ(Imm(3), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := ZMM(), ZMM()
	VMOVDQU64(Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOVDQU64(Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOut := ZMM()
	VPADDQ(x1, x0, xOut)

	if !isLazy {
		subQ, subQMask := ZMM(), K()
		VPCMPQ(Imm(0o5), qv, xOut, subQMask)
		VMOVAPD_Z(qv, subQMask, subQ)
		VPSUBQ(subQ, xOut, xOut)
	}

	VMOVDQU64(xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(8), i)

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
		subQ := GP64()
		MOVQ(y0, subQ)
		SUBQ(q, subQ)
		CMPQ(y0, q)
		CMOVQGE(subQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func SubVecToAVX2(isLazy bool) {
	if isLazy {
		TEXT("subLazyVecToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("subVecToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
	}
	Pragma("noescape")

	zero := YMM()
	VPXOR(zero, zero, zero)

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = YMM()
		VPBROADCASTQ(NewParamAddr("q", 72), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())

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
		subQ := GP64()
		MOVQ(y0, subQ)
		ADDQ(q, subQ)
		CMPQ(y0, Imm(0))
		CMOVQLT(subQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func SubVecToAVX512(isLazy bool) {
	if isLazy {
		TEXT("subLazyVecToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("subVecToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
	}
	Pragma("noescape")

	zero := ZMM()
	VXORPD(zero, zero, zero)

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = ZMM()
		VPBROADCASTQ(NewParamAddr("q", 72), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(3), M)
	SHLQ(Imm(3), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := ZMM(), ZMM()
	VMOVDQU64(Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOVDQU64(Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOut := ZMM()
	VPSUBQ(x1, x0, xOut)

	if !isLazy {
		subQ, subQMask := ZMM(), K()
		VPCMPQ(Imm(0o1), zero, xOut, subQMask)
		VMOVAPD_Z(qv, subQMask, subQ)
		VPADDQ(subQ, xOut, xOut)
	}

	VMOVDQU64(xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(8), i)

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
		subQ := GP64()
		MOVQ(y0, subQ)
		ADDQ(q, subQ)
		CMPQ(y0, Imm(0))
		CMOVQLT(subQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}
