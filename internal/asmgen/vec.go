package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func VecConstants() {
	ConstData("MASK_LO", U64(1<<32-1))
	ConstData("MASK_52", U64(1<<52-1))
}

func VecAddSubToAVX2(opType OpType, isWordOp bool) {
	switch opType {
	case OpAdd:
		if isWordOp {
			TEXT("addWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		} else {
			TEXT("addToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
		}
	case OpSub:
		if isWordOp {
			TEXT("subWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		} else {
			TEXT("subToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
		}
	}
	Pragma("noescape")

	maskSign := YMM()
	if !isWordOp {
		VPCMPEQQ(maskSign, maskSign, maskSign)
		VPSLLQ(Imm(63), maskSign, maskSign)
	}

	q64, q := GP64(), YMM()
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 72), q)
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
	switch opType {
	case OpAdd:
		VPADDQ(x1, x0, xOut)
		if !isWordOp {
			subQ := YMM()
			GreaterOrEqualThanAVX2(xOut, q, maskSign, subQ)
			VPAND(q, subQ, subQ)
			VPSUBQ(subQ, xOut, xOut)
		}
	case OpSub:
		VPSUBQ(x1, x0, xOut)
		if !isWordOp {
			subQ := YMM()
			GreaterOrEqualThanAVX2(xOut, q, maskSign, subQ)
			VPAND(q, subQ, subQ)
			VPADDQ(subQ, xOut, xOut)
		}
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

	switch opType {
	case OpAdd:
		ADDQ(y1, y0)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y0, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, y0)
			CMOVQLS(subQ, y0)
		}
	case OpSub:
		SUBQ(y1, y0)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y0, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, y0)
			CMOVQLS(subQ, y0)
		}
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecAddSubScalarToAVX2(opType OpType, isWordOp bool) {
	switch opType {
	case OpAdd:
		if isWordOp {
			TEXT("addScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		} else {
			TEXT("addScalarToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
		}
	case OpSub:
		if isWordOp {
			TEXT("subScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		} else {
			TEXT("subScalarToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
		}
	}
	Pragma("noescape")

	maskSign := YMM()
	if !isWordOp {
		VPCMPEQQ(maskSign, maskSign, maskSign)
		VPSLLQ(Imm(63), maskSign, maskSign)
	}

	q64, q := GP64(), YMM()
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 56), q)
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
	switch opType {
	case OpAdd:
		VPADDQ(c, x, xOut)
		if !isWordOp {
			subQ := YMM()
			GreaterOrEqualThanAVX2(xOut, q, maskSign, subQ)
			VPAND(q, subQ, subQ)
			VPSUBQ(subQ, xOut, xOut)
		}
	case OpSub:
		VPSUBQ(c, x, xOut)
		if !isWordOp {
			subQ := YMM()
			GreaterOrEqualThanAVX2(xOut, q, maskSign, subQ)
			VPAND(q, subQ, subQ)
			VPADDQ(subQ, xOut, xOut)
		}
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

	switch opType {
	case OpAdd:
		ADDQ(c64, y)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, y)
			CMOVQLS(subQ, y)
		}
	case OpSub:
		SUBQ(c64, y)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, y)
			CMOVQLS(subQ, y)
		}
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecAddSubToAVX512(opType OpType, isWordOp bool) {
	switch opType {
	case OpAdd:
		if isWordOp {
			TEXT("addWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		} else {
			TEXT("addToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
		}
	case OpSub:
		if isWordOp {
			TEXT("subWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		} else {
			TEXT("subToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64)")
		}
	}
	Pragma("noescape")

	q64, q := GP64(), ZMM()
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 72), q)
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
	switch opType {
	case OpAdd:
		VPADDQ(x1, x0, xOut)
		if !isWordOp {
			xSubQ := ZMM()
			VPSUBQ(q, xOut, xSubQ)
			VPMINUQ(xSubQ, xOut, xOut)
		}
	case OpSub:
		VPSUBQ(x1, x0, xOut)
		if !isWordOp {
			xAddQ := ZMM()
			VPADDQ(q, xOut, xAddQ)
			VPMINUQ(xAddQ, xOut, xOut)
		}
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

	switch opType {
	case OpAdd:
		ADDQ(y1, y0)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y0, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, y0)
			CMOVQLS(subQ, y0)
		}
	case OpSub:
		SUBQ(y1, y0)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y0, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, y0)
			CMOVQLS(subQ, y0)
		}
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecAddSubScalarToAVX512(opType OpType, isWordOp bool) {
	switch opType {
	case OpAdd:
		if isWordOp {
			TEXT("addScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		} else {
			TEXT("addScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
		}
	case OpSub:
		if isWordOp {
			TEXT("subScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		} else {
			TEXT("subScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64)")
		}
	}
	Pragma("noescape")

	q64, q := GP64(), ZMM()
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 56), q)
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
	switch opType {
	case OpAdd:
		VPADDQ(c, x, xOut)
		if !isWordOp {
			xSubQ := ZMM()
			VPSUBQ(q, xOut, xSubQ)
			VPMINUQ(xSubQ, xOut, xOut)
		}
	case OpSub:
		VPSUBQ(c, x, xOut)
		if !isWordOp {
			xAddQ := ZMM()
			VPADDQ(q, xOut, xAddQ)
			VPMINUQ(xAddQ, xOut, xOut)
		}
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

	switch opType {
	case OpAdd:
		ADDQ(c64, y)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, y)
			CMOVQLS(subQ, y)
		}
	case OpSub:
		SUBQ(c64, y)
		if !isWordOp {
			subQ := GP64()
			MOVQ(y, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, y)
			CMOVQLS(subQ, y)
		}
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecNegToAVX2(isWordOp bool) {
	if isWordOp {
		TEXT("negWordToAVX2", NOSPLIT, "func(vOut, v []uint64)")
	} else {
		TEXT("negToAVX2", NOSPLIT, "func(vOut, v []uint64, q uint64)")
	}
	Pragma("noescape")

	zero, one := YMM(), YMM()
	if isWordOp {
		VPXOR(zero, zero, zero)
	} else {
		VPCMPEQQ(one, one, one)
	}

	q64, q := GP64(), YMM()
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 48), q)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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
	if isWordOp {
		VPSUBQ(x, zero, xOut)
	} else {
		eq := YMM()
		VPSUBQ(x, q, x)
		EqualAVX2(x, q, eq)
		VPXOR(eq, one, eq)
		VPAND(x, eq, xOut)
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

	if isWordOp {
		NEGQ(y)
	} else {
		CMPQ(y, Imm(0))
		CMOVQEQ(q64, y)
		NEGQ(y)
		ADDQ(q64, y)
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecNegToAVX512(isWordOp bool) {
	if isWordOp {
		TEXT("negWordToAVX512", NOSPLIT, "func(vOut, v []uint64)")
	} else {
		TEXT("negToAVX512", NOSPLIT, "func(vOut, v []uint64, q uint64)")
	}
	Pragma("noescape")

	zero := ZMM()
	if isWordOp {
		VPXORQ(zero, zero, zero)
	}

	q64, q := GP64(), ZMM()
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 48), q)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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
	if isWordOp {
		VPSUBQ(x, zero, xOut)
	} else {
		VPSUBQ(x, q, x)
		eqMask := K()
		VPCMPQ(Imm(0o4), x, q, eqMask)
		VMOVAPD_Z(x, eqMask, xOut)
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

	if isWordOp {
		NEGQ(y)
	} else {
		CMPQ(y, Imm(0))
		CMOVQEQ(q64, y)
		NEGQ(y)
		ADDQ(q64, y)
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMFormToAVX512() {
	TEXT("mFormToAVX512", NOSPLIT, "func(vOut, v []uint64, q, divHi, divLo uint64)")
	Pragma("noescape")

	zero := ZMM()
	VPXORQ(zero, zero, zero)

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	divHi64 := Load(Param("divHi"), GP64())
	divLo64 := Load(Param("divLo"), GP64())

	q, divHi, divLo := ZMM(), ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 48), q)
	VPBROADCASTQ(NewParamAddr("divHi", 56), divHi)
	VPBROADCASTQ(NewParamAddr("divLo", 64), divLo)

	divLoHi := ZMM()
	VPSRLQ(Imm(32), divLo, divLoHi)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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

	xHi := ZMM()
	VPSRLQ(Imm(32), x, xHi)

	xM := ZMM()
	Mul64HiAVX512(x, xHi, divLo, divLoHi, maskLo, xM)
	VPMULLQ(x, divHi, x)
	VPADDQ(x, xM, xM)

	xOutM := ZMM()
	VPMULLQ(xM, q, xOutM)
	VPSUBQ(xOutM, zero, xOutM)

	xSubQ := ZMM()
	VPSUBQ(q, xOutM, xSubQ)
	VPMINUQ(xSubQ, xOutM, xOutM)

	VMOVDQU64(xOutM, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(8), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	yM, yDivHi := GP64(), GP64()
	MOVQ(divLo64, reg.RDX)
	MULXQ(y, yDivHi, yM)

	MOVQ(y, yDivHi)
	IMULQ(divHi64, yDivHi)
	ADDQ(yDivHi, yM)

	NEGQ(yM)
	IMULQ(q64, yM)

	subQ64 := GP64()
	MOVQ(yM, subQ64)
	SUBQ(q64, subQ64)
	CMPQ(q64, yM)
	CMOVQLS(subQ64, yM)

	MOVQ(yM, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecInvMFormToAVX512() {
	TEXT("invMFormToAVX512", NOSPLIT, "func(vOut, v []uint64, q, inv uint64)")
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	inv64 := Load(Param("inv"), GP64())

	q, inv := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 48), q)
	VPBROADCASTQ(NewParamAddr("inv", 56), inv)

	qHi := ZMM()
	VPSRLQ(Imm(32), q, qHi)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(3), M)
	SHLQ(Imm(3), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	xM := ZMM()
	VMOVDQU64(Mem{Base: v, Index: i, Scale: 8}, xM)

	xMInv, xMInvHi := ZMM(), ZMM()
	VPMULLQ(xM, inv, xMInv)
	VPSRLQ(Imm(32), xMInv, xMInvHi)

	xOut := ZMM()
	Mul64HiAVX512(xMInv, xMInvHi, q, qHi, maskLo, xOut)
	VPSUBQ(xOut, q, xOut)

	xSubQ := ZMM()
	VPSUBQ(q, xOut, xSubQ)
	VPMINUQ(xSubQ, xOut, xOut)

	VMOVDQU64(xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(8), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	yM := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, yM)

	y, yLo := GP64(), GP64()
	IMULQ(inv64, yM)
	MOVQ(q64, reg.RDX)
	MULXQ(yM, yLo, y)

	NEGQ(y)
	ADDQ(q64, y)

	subQ64 := GP64()
	MOVQ(y, subQ64)
	SUBQ(q64, subQ64)
	CMPQ(q64, y)
	CMOVQLS(subQ64, y)

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulScalarWordToAVX2(opType OpType) {
	switch opType {
	case OpPure:
		TEXT("mulScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	case OpAdd:
		TEXT("mulAddScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	case OpSub:
		TEXT("mulSubScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	}
	Pragma("noescape")

	maskHi := YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskHi)
	VPSLLQ(Imm(32), maskHi, maskHi)

	c64 := Load(Param("c"), GP64())
	c, cSwap := YMM(), YMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)
	ShuffleForLoAVX2(c, cSwap)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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

	xOut, xMul := YMM(), YMM()
	Mul64LoAVX2(x, c, cSwap, maskHi, xMul)

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
	case OpSub:
		VMOVDQU(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
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

	yOut := GP64()
	IMULQ(c64, y)

	switch opType {
	case OpPure:
		yOut = y
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulScalarWordToAVX512(opType OpType) {
	switch opType {
	case OpPure:
		TEXT("mulScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	case OpAdd:
		TEXT("mulAddScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	case OpSub:
		TEXT("mulSubScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
	}
	Pragma("noescape")

	c64 := Load(Param("c"), GP64())
	c := ZMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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

	xOut, xMul := ZMM(), ZMM()
	VPMULLQ(x, c, xMul)

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
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

	yOut := GP64()
	IMULQ(c64, y)

	switch opType {
	case OpPure:
		yOut = y
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMMulScalarToAVX512(opType OpType, isLazy bool) {
	if !isLazy {
		switch opType {
		case OpPure:
			TEXT("mMulScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q, inv uint64)")
		case OpAdd:
			TEXT("mMulAddScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q, inv uint64)")
		case OpSub:
			TEXT("mMulSubScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q, inv uint64)")
		}
	} else {
		switch opType {
		case OpPure:
			TEXT("mMulScalarLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q, inv uint64)")
		case OpAdd:
			TEXT("mMulAddScalarLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q, inv uint64)")
		case OpSub:
			TEXT("mMulSubScalarLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q, inv uint64)")
		}
	}
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	q, qHi := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 56), q)
	VPSRLQ(Imm(32), q, qHi)

	inv64 := Load(Param("inv"), GP64())
	inv := ZMM()
	VPBROADCASTQ(NewParamAddr("inv", 64), inv)

	c64 := Load(Param("c"), GP64())
	c, cHi := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)
	VPSRLQ(Imm(32), c, cHi)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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

	xOut, xMul := ZMM(), ZMM()

	xMulMHi, xMulMLo := ZMM(), ZMM()
	VPMULLQ(x, c, xMulMLo)
	VPMULLQ(inv, xMulMLo, xMulMLo)

	xHi := ZMM()
	VPSRLQ(Imm(32), x, xHi)
	Mul64HiAVX512(x, xHi, c, cHi, maskLo, xMulMHi)

	xMulMLoHi := ZMM()
	VPSRLQ(Imm(32), xMulMLo, xMulMLoHi)

	wHi := ZMM()
	Mul64HiAVX512(xMulMLo, xMulMLoHi, q, qHi, maskLo, wHi)

	VPSUBQ(wHi, xMulMHi, xMul)
	VPADDQ(q, xMul, xMul)

	if !isLazy {
		xSubQ := ZMM()
		VPSUBQ(q, xMul, xSubQ)
		VPMINUQ(xSubQ, xMul, xMul)
	}

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		if !isLazy {
			xSubQ := ZMM()
			VPSUBQ(q, xOut, xSubQ)
			VPMINUQ(xSubQ, xOut, xOut)
		}
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		if !isLazy {
			VPSUBQ(xMul, xOut, xOut)
			xAddQ := ZMM()
			VPADDQ(q, xOut, xAddQ)
			VPMINUQ(xAddQ, xOut, xOut)
		} else {
			VPADDQ(xMul, xOut, xOut)
		}
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

	yOut := GP64()

	yMulMHi, yMulMLo := GP64(), GP64()
	MOVQ(y, reg.RDX)
	MULXQ(c64, yMulMLo, yMulMHi)

	zHi := GP64()
	IMULQ(inv64, yMulMLo)
	MOVQ(yMulMLo, reg.RDX)
	MULXQ(q64, yMulMLo, zHi)

	MOVQ(q64, y)
	ADDQ(yMulMHi, y)
	SUBQ(zHi, y)

	if !isLazy {
		subQ := GP64()
		MOVQ(y, subQ)
		SUBQ(q64, subQ)
		CMPQ(q64, y)
		CMOVQLS(subQ, y)
	}

	switch opType {
	case OpPure:
		yOut = y
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y, yOut)
		if !isLazy {
			subQ := GP64()
			MOVQ(yOut, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		}
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		if !isLazy {
			SUBQ(y, yOut)
			subQ := GP64()
			MOVQ(yOut, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		} else {
			ADDQ(y, yOut)
		}
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecSMulScalarToAVX512(opType OpType, isLazy bool) {
	if !isLazy {
		switch opType {
		case OpPure:
			TEXT("sMulScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpAdd:
			TEXT("sMulAddScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpSub:
			TEXT("sMulSubScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		}
	} else {
		switch opType {
		case OpPure:
			TEXT("sMulScalarLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpAdd:
			TEXT("sMulAddScalarLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpSub:
			TEXT("sMulSubScalarLazyToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		}
	}
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	q, qHi := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 64), q)
	VPSRLQ(Imm(32), q, qHi)

	c64 := Load(Param("c"), GP64())
	cS64 := Load(Param("cS"), GP64())
	c, cS := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)
	VPBROADCASTQ(NewParamAddr("cS", 56), cS)

	cSHi := ZMM()
	VPSRLQ(Imm(32), cS, cSHi)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

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

	xOut, xMul := ZMM(), ZMM()

	xHi := ZMM()
	VPSRLQ(Imm(32), x, xHi)

	quo := ZMM()
	Mul64HiAVX512(x, xHi, cS, cSHi, maskLo, quo)
	VPMULLQ(quo, q, quo)

	VPMULLQ(x, c, xMul)
	VPSUBQ(quo, xMul, xMul)

	if !isLazy {
		xSubQ := ZMM()
		VPSUBQ(q, xMul, xSubQ)
		VPMINUQ(xSubQ, xMul, xMul)
	}

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		if !isLazy {
			xSubQ := ZMM()
			VPSUBQ(q, xOut, xSubQ)
			VPMINUQ(xSubQ, xOut, xOut)
		}
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		if !isLazy {
			VPSUBQ(xMul, xOut, xOut)
			xAddQ := ZMM()
			VPADDQ(q, xOut, xAddQ)
			VPMINUQ(xAddQ, xOut, xOut)
		} else {
			VPADDQ(xMul, xOut, xOut)
		}
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

	yOut := GP64()

	quo64 := GP64()
	MOVQ(cS64, reg.RDX)
	MULXQ(y, yOut, quo64)
	IMULQ(q64, quo64)

	IMULQ(c64, y)
	SUBQ(quo64, y)

	if !isLazy {
		subQ := GP64()
		MOVQ(y, subQ)
		SUBQ(q64, subQ)
		CMPQ(q64, y)
		CMOVQLS(subQ, y)
	}

	switch opType {
	case OpPure:
		yOut = y
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y, yOut)
		if !isLazy {
			subQ := GP64()
			MOVQ(yOut, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		}
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		if !isLazy {
			SUBQ(y, yOut)
			subQ := GP64()
			MOVQ(yOut, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		} else {
			ADDQ(y, yOut)
		}
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulWordToAVX2(opType OpType) {
	switch opType {
	case OpPure:
		TEXT("mulWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	case OpAdd:
		TEXT("mulAddWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	case OpSub:
		TEXT("mulSubWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	}
	Pragma("noescape")

	maskHi := YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskHi)
	VPSLLQ(Imm(32), maskHi, maskHi)

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

	x1Swap := YMM()
	ShuffleForLoAVX2(x1, x1Swap)

	xOut, xMul := YMM(), YMM()
	Mul64LoAVX2(x0, x1, x1Swap, maskHi, xMul)

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
	case OpSub:
		VMOVDQU(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
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

	yOut := GP64()
	IMULQ(y1, y0)

	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y0, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulWordToAVX512(opType OpType) {
	switch opType {
	case OpPure:
		TEXT("mulWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	case OpAdd:
		TEXT("mulAddWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	case OpSub:
		TEXT("mulSubWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
	}
	Pragma("noescape")

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

	xOut, xMul := ZMM(), ZMM()
	VPMULLQ(x0, x1, xMul)

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
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

	yOut := GP64()
	IMULQ(y1, y0)

	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y0, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMMulToAVX512(opType OpType, isLazy bool) {
	if !isLazy {
		switch opType {
		case OpPure:
			TEXT("mMulToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q, inv uint64)")
		case OpAdd:
			TEXT("mMulAddToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q, inv uint64)")
		case OpSub:
			TEXT("mMulSubToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q, inv uint64)")
		}
	} else {
		switch opType {
		case OpPure:
			TEXT("mMulLazyToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q, inv uint64)")
		case OpAdd:
			TEXT("mMulAddLazyToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q, inv uint64)")
		case OpSub:
			TEXT("mMulSubLazyToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q, inv uint64)")
		}
	}
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	inv64 := Load(Param("inv"), GP64())
	q, qHi := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPSRLQ(Imm(32), q, qHi)
	inv := ZMM()
	VPBROADCASTQ(NewParamAddr("inv", 80), inv)

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

	if opType == OpSub && isLazy {
		VPSUBQ(x0, q, x0)
	}

	xOut, xMul := ZMM(), ZMM()

	xMulMHi, xMulMLo := ZMM(), ZMM()
	VPMULLQ(x0, x1, xMulMLo)
	VPMULLQ(inv, xMulMLo, xMulMLo)

	x0Hi, x1Hi := ZMM(), ZMM()
	VPSRLQ(Imm(32), x0, x0Hi)
	VPSRLQ(Imm(32), x1, x1Hi)
	Mul64HiAVX512(x0, x0Hi, x1, x1Hi, maskLo, xMulMHi)

	xMulMLoHi := ZMM()
	VPSRLQ(Imm(32), xMulMLo, xMulMLoHi)

	wHi := ZMM()
	Mul64HiAVX512(xMulMLo, xMulMLoHi, q, qHi, maskLo, wHi)

	VPSUBQ(wHi, xMulMHi, xMul)
	VPADDQ(q, xMul, xMul)

	if !isLazy {
		xSubQ := ZMM()
		VPSUBQ(q, xMul, xSubQ)
		VPMINUQ(xSubQ, xMul, xMul)
	}

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		if !isLazy {
			xSubQ := ZMM()
			VPSUBQ(q, xOut, xSubQ)
			VPMINUQ(xSubQ, xOut, xOut)
		}
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		if !isLazy {
			VPSUBQ(xMul, xOut, xOut)
			xAddQ := ZMM()
			VPADDQ(q, xOut, xAddQ)
			VPMINUQ(xAddQ, xOut, xOut)
		} else {
			VPADDQ(xMul, xOut, xOut)
		}
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

	if opType == OpSub && isLazy {
		SUBQ(q64, y0)
		NEGQ(y0)
	}

	yOut := GP64()

	yMulMHi, yMulMLo := GP64(), GP64()
	MOVQ(y0, reg.RDX)
	MULXQ(y1, yMulMLo, yMulMHi)

	zHi := GP64()
	IMULQ(inv64, yMulMLo)
	MOVQ(yMulMLo, reg.RDX)
	MULXQ(q64, yMulMLo, zHi)

	MOVQ(q64, y0)
	ADDQ(yMulMHi, y0)
	SUBQ(zHi, y0)

	if !isLazy {
		subQ := GP64()
		MOVQ(y0, subQ)
		SUBQ(q64, subQ)
		CMPQ(q64, y0)
		CMOVQLS(subQ, y0)
	}

	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)
		if !isLazy {
			subQ := GP64()
			MOVQ(yOut, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		}
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		if !isLazy {
			SUBQ(y0, yOut)
			subQ := GP64()
			MOVQ(yOut, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		} else {
			ADDQ(y0, yOut)
		}
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecSMulToAVX512(opType OpType, isLazy bool) {
	if !isLazy {
		switch opType {
		case OpPure:
			TEXT("sMulToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpAdd:
			TEXT("sMulAddToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpSub:
			TEXT("sMulSubToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		}
	} else {
		switch opType {
		case OpPure:
			TEXT("sMulLazyToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpAdd:
			TEXT("sMulAddLazyToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpSub:
			TEXT("sMulSubLazyToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		}
	}
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	q, qHi := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 96), q)
	VPSRLQ(Imm(32), q, qHi)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())
	v1S := Load(Param("v1S").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(3), M)
	SHLQ(Imm(3), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1, x1S := ZMM(), ZMM(), ZMM()
	VMOVDQU64(Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOVDQU64(Mem{Base: v1, Index: i, Scale: 8}, x1)
	VMOVDQU64(Mem{Base: v1S, Index: i, Scale: 8}, x1S)

	if opType == OpSub && isLazy {
		VPSUBQ(x0, q, x0)
	}

	xOut, xMul := ZMM(), ZMM()

	x0Hi, x1SHi := ZMM(), ZMM()
	VPSRLQ(Imm(32), x0, x0Hi)
	VPSRLQ(Imm(32), x1S, x1SHi)

	quo := ZMM()
	Mul64HiAVX512(x0, x0Hi, x1S, x1SHi, maskLo, quo)
	VPMULLQ(quo, q, quo)

	VPMULLQ(x0, x1, xMul)
	VPSUBQ(quo, xMul, xMul)

	if !isLazy {
		xSubQ := ZMM()
		VPSUBQ(q, xMul, xSubQ)
		VPMINUQ(xSubQ, xMul, xMul)
	}

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		if !isLazy {
			xSubQ := ZMM()
			VPSUBQ(q, xOut, xSubQ)
			VPMINUQ(xSubQ, xOut, xOut)
		}
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		if !isLazy {
			VPSUBQ(xMul, xOut, xOut)
			xAddQ := ZMM()
			VPADDQ(q, xOut, xAddQ)
			VPMINUQ(xAddQ, xOut, xOut)
		} else {
			VPADDQ(xMul, xOut, xOut)
		}
	}

	VMOVDQU64(xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(8), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y0, y1, y1S := GP64(), GP64(), GP64()
	MOVQ(Mem{Base: v0, Index: i, Scale: 8}, y0)
	MOVQ(Mem{Base: v1, Index: i, Scale: 8}, y1)
	MOVQ(Mem{Base: v1S, Index: i, Scale: 8}, y1S)

	if opType == OpSub && isLazy {
		SUBQ(q64, y0)
		NEGQ(y0)
	}

	yOut := GP64()

	quo64 := GP64()
	MOVQ(y1S, reg.RDX)
	MULXQ(y0, y1S, quo64)
	IMULQ(q64, quo64)

	IMULQ(y1, y0)
	SUBQ(quo64, y0)

	if !isLazy {
		subQ := GP64()
		MOVQ(y0, subQ)
		SUBQ(q64, subQ)
		CMPQ(q64, y0)
		CMOVQLS(subQ, y0)
	}

	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)
		if !isLazy {
			subQ := GP64()
			MOVQ(yOut, subQ)
			SUBQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		}
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		if !isLazy {
			SUBQ(y0, yOut)
			subQ := GP64()
			MOVQ(yOut, subQ)
			ADDQ(q64, subQ)
			CMPQ(q64, yOut)
			CMOVQLS(subQ, yOut)
		} else {
			ADDQ(y0, yOut)
		}
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}
