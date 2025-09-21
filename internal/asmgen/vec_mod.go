package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func AddVecToAVX2(isLazy bool) {
	if isLazy {
		TEXT("addLazyToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("addToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
	}
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	maskSign := YMM()
	VPSLLQ(Imm(63), allOne, maskSign)

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
		GreaterOrEqualThanAVX2(xOut, qv, maskSign, allOne, subQ)
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
		CMOVQCC(subQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func ScalarAddVecToAVX2(isLazy bool) {
	if isLazy {
		TEXT("scalarAddLazyToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	} else {
		TEXT("scalarAddToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
	}
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	maskSign := YMM()
	VPSLLQ(Imm(63), allOne, maskSign)

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = YMM()
		VPBROADCASTQ(NewParamAddr("q", 56), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	c64 := Load(Param("c"), GP64())
	c := YMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(2), M)
	SHLQ(Imm(2), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := YMM()
	VMOVDQU(Mem{Base: v, Index: i, Scale: 8}, x)

	xOut := YMM()
	VPADDQ(c, x, xOut)

	if !isLazy {
		subQ := YMM()
		GreaterOrEqualThanAVX2(xOut, qv, maskSign, allOne, subQ)
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

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	ADDQ(c64, y)

	if !isLazy {
		subQ := GP64()
		MOVQ(y, subQ)
		SUBQ(q, subQ)
		CMPQ(y, q)
		CMOVQCC(subQ, y)
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func AddVecToAVX512(isLazy bool) {
	if isLazy {
		TEXT("addLazyToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("addToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
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
		VPCMPUQ(Imm(0o5), qv, xOut, subQMask)
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
		CMOVQCC(subQ, y0)
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func ScalarAddVecToAVX512(isLazy bool) {
	if isLazy {
		TEXT("scalarAddLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	} else {
		TEXT("scalarAddToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
	}
	Pragma("noescape")

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = ZMM()
		VPBROADCASTQ(NewParamAddr("q", 56), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	c64 := Load(Param("c"), GP64())
	c := ZMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(3), M)
	SHLQ(Imm(3), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := ZMM()
	VMOVDQU64(Mem{Base: v, Index: i, Scale: 8}, x)

	xOut := ZMM()
	VPADDQ(c, x, xOut)

	if !isLazy {
		subQ, subQMask := ZMM(), K()
		VPCMPUQ(Imm(0o5), qv, xOut, subQMask)
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

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	ADDQ(c64, y)

	if !isLazy {
		subQ := GP64()
		MOVQ(y, subQ)
		SUBQ(q, subQ)
		CMPQ(y, q)
		CMOVQCC(subQ, y)
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func SubVecToAVX2(isLazy bool) {
	if isLazy {
		TEXT("subLazyToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("subToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
	}
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	maskSign := YMM()
	VPSLLQ(Imm(63), allOne, maskSign)

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
		GreaterOrEqualThanAVX2(xOut, qv, maskSign, allOne, subQ)
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

func ScalarSubVecToAVX2(isLazy bool) {
	if isLazy {
		TEXT("scalarSubLazyToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	} else {
		TEXT("scalarSubToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
	}
	Pragma("noescape")

	allOne := YMM()
	VPCMPEQQ(allOne, allOne, allOne)

	maskSign := YMM()
	VPSLLQ(Imm(63), allOne, maskSign)

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = YMM()
		VPBROADCASTQ(NewParamAddr("q", 56), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	c64 := Load(Param("c"), GP64())
	c := YMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(2), M)
	SHLQ(Imm(2), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := YMM()
	VMOVDQU(Mem{Base: v, Index: i, Scale: 8}, x)

	xOut := YMM()
	VPSUBQ(c, x, xOut)

	if !isLazy {
		subQ := YMM()
		GreaterOrEqualThanAVX2(xOut, qv, maskSign, allOne, subQ)
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

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	SUBQ(c64, y)

	if !isLazy {
		subQ := GP64()
		MOVQ(y, subQ)
		ADDQ(q, subQ)
		CMPQ(y, Imm(0))
		CMOVQLT(subQ, y)
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func SubVecToAVX512(isLazy bool) {
	if isLazy {
		TEXT("subLazyToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	} else {
		TEXT("subToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
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
	VPSUBQ(x1, x0, xOut)

	if !isLazy {
		subQ, subQMask := ZMM(), K()
		VPCMPUQ(Imm(0o5), qv, xOut, subQMask)
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

func ScalarSubVecToAVX512(isLazy bool) {
	if isLazy {
		TEXT("scalarSubLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	} else {
		TEXT("scalarSubToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
	}
	Pragma("noescape")

	var q reg.Register
	var qv reg.VecVirtual
	if !isLazy {
		q = Load(Param("q"), GP64())
		qv = ZMM()
		VPBROADCASTQ(NewParamAddr("q", 56), qv)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	c64 := Load(Param("c"), GP64())
	c := ZMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(3), M)
	SHLQ(Imm(3), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := ZMM()
	VMOVDQU64(Mem{Base: v, Index: i, Scale: 8}, x)

	xOut := ZMM()
	VPSUBQ(c, x, xOut)

	if !isLazy {
		subQ, subQMask := ZMM(), K()
		VPCMPUQ(Imm(0o5), qv, xOut, subQMask)
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

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	SUBQ(c64, y)

	if !isLazy {
		subQ := GP64()
		MOVQ(y, subQ)
		ADDQ(q, subQ)
		CMPQ(y, Imm(0))
		CMOVQLT(subQ, y)
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}
