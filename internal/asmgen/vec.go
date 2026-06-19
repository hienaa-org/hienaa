package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func VecConstants() {
	ConstData("MASK_LO", U64(1<<32-1))
	ConstData("MASK_52", U64(1<<52-1))
	ConstData("CVT_52", F64(1<<52))
}

func VecAddSubToAVX(avxType AVXType, opType OpType, isWordOp bool) {
	switch avxType {
	case TypeAVX2:
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
	case TypeAVX512:
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
	}
	Pragma("noescape")

	q64, q := GP64(), VMM(avxType)
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
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := VMM(avxType), VMM(avxType)
	VMOV(avxType, Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOV(avxType, Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOut := VMM(avxType)

	xOutQ := VMM(avxType)
	switch opType {
	case OpAdd:
		VPADDQ(x1, x0, xOut)
		if !isWordOp {
			switch avxType {
			case TypeAVX2:
				VPSUBQ(q, xOut, xOutQ)
				VBLENDVPD(xOutQ, xOut, xOutQ, xOut)
			case TypeAVX512, TypeAVX512IFMA:
				VPSUBQ(q, xOut, xOutQ)
				VPMINUQ(xOutQ, xOut, xOut)
			}
		}
	case OpSub:
		VPSUBQ(x1, x0, xOut)
		if !isWordOp {
			switch avxType {
			case TypeAVX2:
				VPADDQ(q, xOut, xOutQ)
				VBLENDVPD(xOut, xOutQ, xOut, xOut)
			case TypeAVX512, TypeAVX512IFMA:
				VPADDQ(q, xOut, xOutQ)
				VPMINUQ(xOutQ, xOut, xOut)
			}
		}
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y0, y1 := GP64(), GP64()
	MOVQ(Mem{Base: v0, Index: i, Scale: 8}, y0)
	MOVQ(Mem{Base: v1, Index: i, Scale: 8}, y1)

	y0Q := GP64()
	switch opType {
	case OpAdd:
		ADDQ(y1, y0)
		if !isWordOp {
			MOVQ(y0, y0Q)
			SUBQ(q64, y0Q)
			CMPQ(q64, y0)
			CMOVQLS(y0Q, y0)
		}
	case OpSub:
		SUBQ(y1, y0)
		if !isWordOp {
			MOVQ(y0, y0Q)
			ADDQ(q64, y0Q)
			CMPQ(q64, y0)
			CMOVQLS(y0Q, y0)
		}
	}

	MOVQ(y0, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecAddSubScalarToAVX(avxType AVXType, opType OpType, isWordOp bool) {
	switch avxType {
	case TypeAVX2:
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
	case TypeAVX512:
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
	}
	Pragma("noescape")

	q64, q := GP64(), VMM(avxType)
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 56), q)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	c64 := Load(Param("c"), GP64())
	c := VMM(avxType)
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := VMM(avxType)
	VMOV(avxType, Mem{Base: v, Index: i, Scale: 8}, x)

	xOut := VMM(avxType)

	xOutQ := VMM(avxType)
	switch opType {
	case OpAdd:
		VPADDQ(c, x, xOut)
		if !isWordOp {
			switch avxType {
			case TypeAVX2:
				VPSUBQ(q, xOut, xOutQ)
				VBLENDVPD(xOutQ, xOut, xOutQ, xOut)
			case TypeAVX512, TypeAVX512IFMA:
				VPSUBQ(q, xOut, xOutQ)
				VPMINUQ(xOutQ, xOut, xOut)
			}
		}
	case OpSub:
		VPSUBQ(c, x, xOut)
		if !isWordOp {
			switch avxType {
			case TypeAVX2:
				VPADDQ(q, xOut, xOutQ)
				VBLENDVPD(xOut, xOutQ, xOut, xOut)
			case TypeAVX512, TypeAVX512IFMA:
				VPADDQ(q, xOut, xOutQ)
				VPMINUQ(xOutQ, xOut, xOut)
			}
		}
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	yQ := GP64()
	switch opType {
	case OpAdd:
		ADDQ(c64, y)
		if !isWordOp {
			MOVQ(y, yQ)
			SUBQ(q64, yQ)
			CMPQ(q64, y)
			CMOVQLS(yQ, y)
		}
	case OpSub:
		SUBQ(c64, y)
		if !isWordOp {
			MOVQ(y, yQ)
			ADDQ(q64, yQ)
			CMPQ(q64, y)
			CMOVQLS(yQ, y)
		}
	}

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecNegToAVX(avxType AVXType, isWordOp bool) {
	switch avxType {
	case TypeAVX2:
		if isWordOp {
			TEXT("negWordToAVX2", NOSPLIT, "func(vOut, v []uint64)")
		} else {
			TEXT("negToAVX2", NOSPLIT, "func(vOut, v []uint64, q uint64)")
		}
	case TypeAVX512:
		if isWordOp {
			TEXT("negWordToAVX512", NOSPLIT, "func(vOut, v []uint64)")
		} else {
			TEXT("negToAVX512", NOSPLIT, "func(vOut, v []uint64, q uint64)")
		}
	}
	Pragma("noescape")

	zero := VMM(avxType)
	if isWordOp {
		VXORPD(zero, zero, zero)
	}

	one := VMM(avxType)
	if avxType == TypeAVX2 && !isWordOp {
		VPCMPEQQ(one, one, one)
	}

	q64, q := GP64(), VMM(avxType)
	if !isWordOp {
		Load(Param("q"), q64)
		VPBROADCASTQ(NewParamAddr("q", 48), q)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := VMM(avxType)
	VMOV(avxType, Mem{Base: v, Index: i, Scale: 8}, x)

	xOut := VMM(avxType)
	if isWordOp {
		VPSUBQ(x, zero, xOut)
	} else {
		switch avxType {
		case TypeAVX2:
			eq := VMM(avxType)
			VPSUBQ(x, q, x)
			VPCMPEQQ(x, q, eq)
			VPXOR(eq, one, eq)
			VPAND(x, eq, xOut)
		case TypeAVX512, TypeAVX512IFMA:
			VPSUBQ(x, q, x)
			eqMask := K()
			VPCMPQ(Imm(0o4), x, q, eqMask)
			VMOVAPD_Z(x, eqMask, xOut)
		}
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

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

func VecMulScalarWordToAVX(avxType AVXType, opType OpType) {
	switch avxType {
	case TypeAVX2:
		switch opType {
		case OpPure:
			TEXT("mulScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		case OpAdd:
			TEXT("mulAddScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		case OpSub:
			TEXT("mulSubScalarWordToAVX2", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		}
	case TypeAVX512:
		switch opType {
		case OpPure:
			TEXT("mulScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		case OpAdd:
			TEXT("mulAddScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		case OpSub:
			TEXT("mulSubScalarWordToAVX512", NOSPLIT, "func(vOut, v []uint64, c uint64)")
		}
	}

	Pragma("noescape")

	maskHi := VMM(avxType)
	if avxType == TypeAVX2 {
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskHi)
		VPSLLQ(Imm(32), maskHi, maskHi)
	}

	c64 := Load(Param("c"), GP64())
	c := VMM(avxType)
	VPBROADCASTQ(NewParamAddr("c", 48), c)

	cSwap := VMM(avxType)
	if avxType == TypeAVX2 {
		ShuffleForLoAVX2(c, cSwap)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := VMM(avxType)
	VMOV(avxType, Mem{Base: v, Index: i, Scale: 8}, x)

	xOut, xMul := VMM(avxType), VMM(avxType)
	switch avxType {
	case TypeAVX2:
		Mul64LoAVX2(x, c, cSwap, maskHi, xMul)
	case TypeAVX512:
		VPMULLQ(c, x, xMul)
	}

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
	case OpSub:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

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

func VecMulScalarToAVX(avxType AVXType, opType OpType) {
	switch avxType {
	case TypeAVX2:
		switch opType {
		case OpPure:
			TEXT("mulScalarToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64, qf, qfInv float64)")
		case OpAdd:
			TEXT("mulAddScalarToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64, qf, qfInv float64)")
		case OpSub:
			TEXT("mulSubScalarToAVX2", NOSPLIT, "func(vOut, v []uint64, c, q uint64, qf, qfInv float64)")
		}
	case TypeAVX512:
		switch opType {
		case OpPure:
			TEXT("mulScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64, qf, qfInv float64)")
		case OpAdd:
			TEXT("mulAddScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64, qf, qfInv float64)")
		case OpSub:
			TEXT("mulSubScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, q uint64, qf, qfInv float64)")
		}
	}
	Pragma("noescape")

	cvt52 := VMM(avxType)
	if avxType == TypeAVX2 {
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("CVT_52"), 0), cvt52)
	}

	q64 := Load(Param("q"), GP64())
	qf64 := Load(Param("qf"), XMM())
	qfInv64 := Load(Param("qfInv"), XMM())
	q, qf, qfInv := VMM(avxType), VMM(avxType), VMM(avxType)
	VPBROADCASTQ(NewParamAddr("q", 56), q)
	VPBROADCASTQ(NewParamAddr("qf", 64), qf)
	VPBROADCASTQ(NewParamAddr("qfInv", 72), qfInv)

	c64 := Load(Param("c"), GP64())
	cf64 := XMM()
	CVTSQ2SD(c64, cf64)

	c := VMM(avxType)
	VPBROADCASTQ(NewParamAddr("c", 48), c)
	switch avxType {
	case TypeAVX2:
		VORPD(cvt52, c, c)
		VSUBPD(cvt52, c, c)
	case TypeAVX512:
		VCVTUQQ2PD(c, c)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := VMM(avxType)
	VMOV(avxType, Mem{Base: v, Index: i, Scale: 8}, x)

	xOut, xMul := VMM(avxType), VMM(avxType)

	switch avxType {
	case TypeAVX2:
		VORPD(cvt52, x, x)
		VSUBPD(cvt52, x, x)
	case TypeAVX512:
		VCVTUQQ2PD(x, x)
	}

	VMULPD(x, c, xMul)

	VFMSUB213PD(xMul, c, x)
	lo := x

	quo := VMM(avxType)
	VMULPD(xMul, qfInv, quo)
	switch avxType {
	case TypeAVX2:
		VROUNDPD(Imm(1), quo, quo)
	case TypeAVX512:
		VRNDSCALEPD(Imm(1), quo, quo)
	}

	VFNMADD231PD(quo, qf, xMul)
	VADDPD(xMul, lo, xMul)
	VADDPD(xMul, qf, xMul)

	switch avxType {
	case TypeAVX2:
		VADDPD(cvt52, xMul, xMul)
		VXORPD(cvt52, xMul, xMul)
	case TypeAVX512:
		VCVTPD2UQQ(xMul, xMul)
	}

	xMulQ := VMM(avxType)
	switch avxType {
	case TypeAVX2:
		VPSUBQ(q, xMul, xMulQ)
		VBLENDVPD(xMulQ, xMul, xMulQ, xMul)
	case TypeAVX512:
		VPSUBQ(q, xMul, xMulQ)
		VPMINUQ(xMulQ, xMul, xMul)
	}

	xOutQ := VMM(avxType)
	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		switch avxType {
		case TypeAVX2:
			VPSUBQ(q, xOut, xOutQ)
			VBLENDVPD(xOutQ, xOut, xOutQ, xOut)
		case TypeAVX512:
			VPSUBQ(q, xOut, xOutQ)
			VPMINUQ(xOutQ, xOut, xOut)
		}
	case OpSub:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
		switch avxType {
		case TypeAVX2:
			VPADDQ(q, xOut, xOutQ)
			VBLENDVPD(xOut, xOutQ, xOut, xOut)
		case TypeAVX512:
			VPADDQ(q, xOut, xOutQ)
			VPMINUQ(xOutQ, xOut, xOut)
		}
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	yOut := GP64()

	yf := XMM()
	CVTSQ2SD(y, yf)

	yMulf := XMM()
	MOVSD(yf, yMulf)
	MULSD(cf64, yMulf)

	VFMSUB213SD(yMulf, cf64, yf)
	lo64 := yf

	quo64 := XMM()
	MOVSD(yMulf, quo64)
	MULSD(qfInv64, quo64)
	ROUNDSD(Imm(1), quo64, quo64)

	VFNMADD231SD(quo64, qf64, yMulf)
	ADDSD(lo64, yMulf)
	ADDSD(qf64, yMulf)

	CVTTSD2SQ(yMulf, y)

	yQ := GP64()
	MOVQ(y, yQ)
	SUBQ(q64, yQ)
	CMPQ(q64, y)
	CMOVQLS(yQ, y)

	yOutQ := GP64()
	switch opType {
	case OpPure:
		yOut = y
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y, yOut)
		MOVQ(yOut, yOutQ)
		SUBQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y, yOut)
		MOVQ(yOut, yOutQ)
		ADDQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecSMulScalarToAVX(avxType AVXType, opType OpType) {
	switch avxType {
	case TypeAVX512:
		switch opType {
		case OpPure:
			TEXT("sMulScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpAdd:
			TEXT("sMulAddScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpSub:
			TEXT("sMulSubScalarToAVX512", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		}
	case TypeAVX512IFMA:
		switch opType {
		case OpPure:
			TEXT("sMulScalarToAVX512IFMA", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpAdd:
			TEXT("sMulAddScalarToAVX512IFMA", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		case OpSub:
			TEXT("sMulSubScalarToAVX512IFMA", NOSPLIT, "func(vOut, v []uint64, c, cS, q uint64)")
		}
	}
	Pragma("noescape")

	maskLo, mask52, zero := ZMM(), ZMM(), ZMM()
	switch avxType {
	case TypeAVX512:
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)
	case TypeAVX512IFMA:
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_52"), 0), mask52)
		VPXORQ(zero, zero, zero)
	}

	q64 := Load(Param("q"), GP64())
	q := ZMM()
	VPBROADCASTQ(NewParamAddr("q", 64), q)

	c64 := Load(Param("c"), GP64())
	cS64 := Load(Param("cS"), GP64())
	c, cS := ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("c", 48), c)
	VPBROADCASTQ(NewParamAddr("cS", 56), cS)

	cSHi := ZMM()
	switch avxType {
	case TypeAVX512:
		VPSRLQ(Imm(32), cS, cSHi)
	case TypeAVX512IFMA:
		VPSRLQ(Imm(12), cS, cS)
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

	xOut, xMul := ZMM(), ZMM()

	quo := ZMM()
	switch avxType {
	case TypeAVX512:
		xHi := ZMM()
		VPSRLQ(Imm(32), x, xHi)

		Mul64HiAVX512(x, xHi, cS, cSHi, maskLo, quo)
		VPMULLQ(q, quo, quo)
		VPMULLQ(c, x, xMul)
		VPSUBQ(quo, xMul, xMul)

	case TypeAVX512IFMA:
		VPXORQ(quo, quo, quo)
		VPMADD52HUQ(cS, x, quo)

		VPXORQ(xMul, xMul, xMul)
		VPMADD52LUQ(q, quo, xMul)

		VPSUBQ(xMul, zero, xMul)
		VPMADD52LUQ(c, x, xMul)
		VPANDQ(xMul, mask52, xMul)
	}

	xMulQ := ZMM()
	VPSUBQ(q, xMul, xMulQ)
	VPMINUQ(xMulQ, xMul, xMul)

	xOutQ := ZMM()
	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		VPSUBQ(q, xOut, xOutQ)
		VPMINUQ(xOutQ, xOut, xOut)
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
		VPADDQ(q, xOut, xOutQ)
		VPMINUQ(xOutQ, xOut, xOut)
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

	yQ := GP64()
	MOVQ(y, yQ)
	SUBQ(q64, yQ)
	CMPQ(q64, y)
	CMOVQLS(yQ, y)

	yOutQ := GP64()
	switch opType {
	case OpPure:
		yOut = y
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y, yOut)
		MOVQ(yOut, yOutQ)
		SUBQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y, yOut)
		MOVQ(yOut, yOutQ)
		ADDQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulWordToAVX(avxType AVXType, opType OpType) {
	switch avxType {
	case TypeAVX2:
		switch opType {
		case OpPure:
			TEXT("mulWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		case OpAdd:
			TEXT("mulAddWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		case OpSub:
			TEXT("mulSubWordToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		}
	case TypeAVX512:
		switch opType {
		case OpPure:
			TEXT("mulWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		case OpAdd:
			TEXT("mulAddWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		case OpSub:
			TEXT("mulSubWordToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64)")
		}

	}
	Pragma("noescape")

	maskHi := VMM(avxType)
	if avxType == TypeAVX2 {
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskHi)
		VPSLLQ(Imm(32), maskHi, maskHi)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := VMM(avxType), VMM(avxType)
	VMOV(avxType, Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOV(avxType, Mem{Base: v1, Index: i, Scale: 8}, x1)

	x1Swap := VMM(avxType)
	if avxType == TypeAVX2 {
		ShuffleForLoAVX2(x1, x1Swap)
	}

	xOut, xMul := VMM(avxType), VMM(avxType)
	switch avxType {
	case TypeAVX2:
		Mul64LoAVX2(x0, x1, x1Swap, maskHi, xMul)
	case TypeAVX512:
		VPMULLQ(x1, x0, xMul)
	}

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
	case OpSub:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

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

func VecMulToAVX(avxType AVXType, opType OpType) {
	switch avxType {
	case TypeAVX2:
		switch opType {
		case OpPure:
			TEXT("mulToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
		case OpAdd:
			TEXT("mulAddToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
		case OpSub:
			TEXT("mulSubToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
		}
	case TypeAVX512:
		switch opType {
		case OpPure:
			TEXT("mulToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
		case OpAdd:
			TEXT("mulAddToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
		case OpSub:
			TEXT("mulSubToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
		}
	}
	Pragma("noescape")

	cvt52 := VMM(avxType)
	if avxType == TypeAVX2 {
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("CVT_52"), 0), cvt52)
	}

	q64 := Load(Param("q"), GP64())
	qf64 := Load(Param("qf"), XMM())
	qfInv64 := Load(Param("qfInv"), XMM())
	q, qf, qfInv := VMM(avxType), VMM(avxType), VMM(avxType)
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPBROADCASTQ(NewParamAddr("qf", 80), qf)
	VPBROADCASTQ(NewParamAddr("qfInv", 88), qfInv)

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v0 := Load(Param("v0").Base(), GP64())
	v1 := Load(Param("v1").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x0, x1 := VMM(avxType), VMM(avxType)
	VMOV(avxType, Mem{Base: v0, Index: i, Scale: 8}, x0)
	VMOV(avxType, Mem{Base: v1, Index: i, Scale: 8}, x1)

	xOut, xMul := VMM(avxType), VMM(avxType)

	switch avxType {
	case TypeAVX2:
		VORPD(cvt52, x0, x0)
		VSUBPD(cvt52, x0, x0)
		VORPD(cvt52, x1, x1)
		VSUBPD(cvt52, x1, x1)
	case TypeAVX512:
		VCVTUQQ2PD(x0, x0)
		VCVTUQQ2PD(x1, x1)
	}

	VMULPD(x0, x1, xMul)

	VFMSUB213PD(xMul, x1, x0)
	lo := x0

	quo := VMM(avxType)
	VMULPD(xMul, qfInv, quo)
	switch avxType {
	case TypeAVX2:
		VROUNDPD(Imm(1), quo, quo)
	case TypeAVX512:
		VRNDSCALEPD(Imm(1), quo, quo)
	}

	VFNMADD231PD(quo, qf, xMul)
	VADDPD(xMul, lo, xMul)
	VADDPD(xMul, qf, xMul)

	switch avxType {
	case TypeAVX2:
		VADDPD(cvt52, xMul, xMul)
		VXORPD(cvt52, xMul, xMul)
	case TypeAVX512:
		VCVTPD2UQQ(xMul, xMul)
	}

	xMulQ := VMM(avxType)
	switch avxType {
	case TypeAVX2:
		VPSUBQ(q, xMul, xMulQ)
		VBLENDVPD(xMulQ, xMul, xMulQ, xMul)
	case TypeAVX512:
		VPSUBQ(q, xMul, xMulQ)
		VPMINUQ(xMulQ, xMul, xMul)
	}

	xOutQ := VMM(avxType)
	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		switch avxType {
		case TypeAVX2:
			VPSUBQ(q, xOut, xOutQ)
			VBLENDVPD(xOutQ, xOut, xOutQ, xOut)
		case TypeAVX512:
			VPSUBQ(q, xOut, xOutQ)
			VPMINUQ(xOutQ, xOut, xOut)
		}
	case OpSub:
		VMOV(avxType, Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
		switch avxType {
		case TypeAVX2:
			VPADDQ(q, xOut, xOutQ)
			VBLENDVPD(xOut, xOutQ, xOut, xOut)
		case TypeAVX512:
			VPADDQ(q, xOut, xOutQ)
			VPMINUQ(xOutQ, xOut, xOut)
		}
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y0, y1 := GP64(), GP64()
	MOVQ(Mem{Base: v0, Index: i, Scale: 8}, y0)
	MOVQ(Mem{Base: v1, Index: i, Scale: 8}, y1)

	yOut := GP64()

	y0f, y1f := XMM(), XMM()
	CVTSQ2SD(y0, y0f)
	CVTSQ2SD(y1, y1f)

	yMulf := XMM()
	MOVSD(y0f, yMulf)
	MULSD(y1f, yMulf)

	VFMSUB213SD(yMulf, y1f, y0f)
	lo64 := y0f

	quo64 := XMM()
	MOVSD(yMulf, quo64)
	MULSD(qfInv64, quo64)
	ROUNDSD(Imm(1), quo64, quo64)

	VFNMADD231SD(quo64, qf64, yMulf)
	ADDSD(lo64, yMulf)
	ADDSD(qf64, yMulf)

	CVTTSD2SQ(yMulf, y0)

	y0Q := GP64()
	MOVQ(y0, y0Q)
	SUBQ(q64, y0Q)
	CMPQ(q64, y0)
	CMOVQLS(y0Q, y0)

	yOutQ := GP64()
	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)
		MOVQ(yOut, yOutQ)
		SUBQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y0, yOut)
		MOVQ(yOut, yOutQ)
		ADDQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulToAVX2(opType OpType) {
	switch opType {
	case OpPure:
		TEXT("mulToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
	case OpAdd:
		TEXT("mulAddToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
	case OpSub:
		TEXT("mulSubToAVX2", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
	}
	Pragma("noescape")

	cvt52 := YMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("CVT_52"), 0), cvt52)

	q64 := Load(Param("q"), GP64())
	qf64 := Load(Param("qf"), XMM())
	qfInv64 := Load(Param("qfInv"), XMM())
	q, qf, qfInv := YMM(), YMM(), YMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPBROADCASTQ(NewParamAddr("qf", 80), qf)
	VPBROADCASTQ(NewParamAddr("qfInv", 88), qfInv)

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

	xOut, xMul := YMM(), YMM()

	VORPD(cvt52, x0, x0)
	VSUBPD(cvt52, x0, x0)
	VORPD(cvt52, x1, x1)
	VSUBPD(cvt52, x1, x1)

	VMULPD(x0, x1, xMul)

	VFMSUB213PD(xMul, x1, x0)
	lo := x0

	quo := YMM()
	VMULPD(xMul, qfInv, quo)
	VROUNDPD(Imm(1), quo, quo)

	VFNMADD231PD(quo, qf, xMul)
	VADDPD(xMul, lo, xMul)
	VADDPD(xMul, qf, xMul)

	VADDPD(cvt52, xMul, xMul)
	VXORPD(cvt52, xMul, xMul)

	xSubQ := YMM()
	VPSUBQ(q, xMul, xSubQ)
	VBLENDVPD(xSubQ, xMul, xSubQ, xMul)

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		xSubQ := YMM()
		VPSUBQ(q, xOut, xSubQ)
		VBLENDVPD(xSubQ, xOut, xSubQ, xOut)
	case OpSub:
		VMOVDQU(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
		xAddQ := YMM()
		VPADDQ(q, xOut, xAddQ)
		VBLENDVPD(xOut, xAddQ, xOut, xOut)
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

	y0f, y1f := XMM(), XMM()
	CVTSQ2SD(y0, y0f)
	CVTSQ2SD(y1, y1f)

	yMulf := XMM()
	MOVSD(y0f, yMulf)
	MULSD(y1f, yMulf)

	VFMSUB213SD(yMulf, y1f, y0f)
	lo64 := y0f

	quo64 := XMM()
	MOVSD(yMulf, quo64)
	MULSD(qfInv64, quo64)
	ROUNDSD(Imm(1), quo64, quo64)

	VFNMADD231SD(quo64, qf64, yMulf)
	ADDSD(lo64, yMulf)
	ADDSD(qf64, yMulf)

	CVTTSD2SQ(yMulf, y0)

	subQ := GP64()
	MOVQ(y0, subQ)
	SUBQ(q64, subQ)
	CMPQ(q64, y0)
	CMOVQLS(subQ, y0)

	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)

		subQ := GP64()
		MOVQ(yOut, subQ)
		SUBQ(q64, subQ)
		CMPQ(q64, yOut)
		CMOVQLS(subQ, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y0, yOut)

		subQ := GP64()
		MOVQ(yOut, subQ)
		ADDQ(q64, subQ)
		CMPQ(q64, yOut)
		CMOVQLS(subQ, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecMulToAVX512(opType OpType) {
	switch opType {
	case OpPure:
		TEXT("mulToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
	case OpAdd:
		TEXT("mulAddToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
	case OpSub:
		TEXT("mulSubToAVX512", NOSPLIT, "func(vOut, v0, v1 []uint64, q uint64, qf, qfInv float64)")
	}
	Pragma("noescape")

	q64 := Load(Param("q"), GP64())
	qf64 := Load(Param("qf"), XMM())
	qfInv64 := Load(Param("qfInv"), XMM())
	q, qf, qfInv := ZMM(), ZMM(), ZMM()
	VPBROADCASTQ(NewParamAddr("q", 72), q)
	VPBROADCASTQ(NewParamAddr("qf", 80), qf)
	VPBROADCASTQ(NewParamAddr("qfInv", 88), qfInv)

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

	VCVTUQQ2PD(x0, x0)
	VCVTUQQ2PD(x1, x1)

	VMULPD(x0, x1, xMul)

	VFMSUB213PD(xMul, x1, x0)
	lo := x0

	quo := ZMM()
	VMULPD(xMul, qfInv, quo)
	VRNDSCALEPD(Imm(1), quo, quo)

	VFNMADD231PD(quo, qf, xMul)
	VADDPD(xMul, lo, xMul)
	VADDPD(xMul, qf, xMul)

	VCVTPD2UQQ(xMul, xMul)

	xSubQ := ZMM()
	VPSUBQ(q, xMul, xSubQ)
	VPMINUQ(xSubQ, xMul, xMul)

	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		xSubQ := ZMM()
		VPSUBQ(q, xOut, xSubQ)
		VPMINUQ(xSubQ, xOut, xOut)
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
		xAddQ := ZMM()
		VPADDQ(q, xOut, xAddQ)
		VPMINUQ(xAddQ, xOut, xOut)
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

	y0f, y1f := XMM(), XMM()
	CVTSQ2SD(y0, y0f)
	CVTSQ2SD(y1, y1f)

	yMulf := XMM()
	MOVSD(y0f, yMulf)
	MULSD(y1f, yMulf)

	VFMSUB213SD(yMulf, y1f, y0f)
	lo64 := y0f

	quo64 := XMM()
	MOVSD(yMulf, quo64)
	MULSD(qfInv64, quo64)
	ROUNDSD(Imm(1), quo64, quo64)

	VFNMADD231SD(quo64, qf64, yMulf)
	ADDSD(lo64, yMulf)
	ADDSD(qf64, yMulf)

	CVTTSD2SQ(yMulf, y0)

	subQ := GP64()
	MOVQ(y0, subQ)
	SUBQ(q64, subQ)
	CMPQ(q64, y0)
	CMOVQLS(subQ, y0)

	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)

		subQ := GP64()
		MOVQ(yOut, subQ)
		SUBQ(q64, subQ)
		CMPQ(q64, yOut)
		CMOVQLS(subQ, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y0, yOut)

		subQ := GP64()
		MOVQ(yOut, subQ)
		ADDQ(q64, subQ)
		CMPQ(q64, yOut)
		CMOVQLS(subQ, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecSMulToAVX(avxType AVXType, opType OpType) {
	switch avxType {
	case TypeAVX512:
		switch opType {
		case OpPure:
			TEXT("sMulToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpAdd:
			TEXT("sMulAddToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpSub:
			TEXT("sMulSubToAVX512", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		}
	case TypeAVX512IFMA:
		switch opType {
		case OpPure:
			TEXT("sMulToAVX512IFMA", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpAdd:
			TEXT("sMulAddToAVX512IFMA", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		case OpSub:
			TEXT("sMulSubToAVX512IFMA", NOSPLIT, "func(vOut, v0, v1, v1S []uint64, q uint64)")
		}
	}
	Pragma("noescape")

	maskLo, mask52, zero := ZMM(), ZMM(), ZMM()
	switch avxType {
	case TypeAVX512:
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)
	case TypeAVX512IFMA:
		VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_52"), 0), mask52)
		VPXORQ(zero, zero, zero)
	}

	q64 := Load(Param("q"), GP64())
	q := ZMM()
	VPBROADCASTQ(NewParamAddr("q", 96), q)

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

	xOut, xMul := ZMM(), ZMM()

	quo := ZMM()
	switch avxType {
	case TypeAVX512:
		x0Hi, x1SHi := ZMM(), ZMM()
		VPSRLQ(Imm(32), x0, x0Hi)
		VPSRLQ(Imm(32), x1S, x1SHi)

		Mul64HiAVX512(x0, x0Hi, x1S, x1SHi, maskLo, quo)
		VPMULLQ(q, quo, quo)

		VPMULLQ(x1, x0, xMul)
		VPSUBQ(quo, xMul, xMul)

	case TypeAVX512IFMA:
		VPSRLQ(Imm(12), x1S, x1S)

		VPXORQ(quo, quo, quo)
		VPMADD52HUQ(x1S, x0, quo)

		VPXORQ(xMul, xMul, xMul)
		VPMADD52LUQ(q, quo, xMul)

		VPSUBQ(xMul, zero, xMul)
		VPMADD52LUQ(x1, x0, xMul)
		VPANDQ(xMul, mask52, xMul)
	}

	xMulQ := ZMM()
	VPSUBQ(q, xMul, xMulQ)
	VPMINUQ(xMulQ, xMul, xMul)

	xOutQ := ZMM()
	switch opType {
	case OpPure:
		xOut = xMul
	case OpAdd:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPADDQ(xMul, xOut, xOut)
		VPSUBQ(q, xOut, xOutQ)
		VPMINUQ(xOutQ, xOut, xOut)
	case OpSub:
		VMOVDQU64(Mem{Base: vOut, Index: i, Scale: 8}, xOut)
		VPSUBQ(xMul, xOut, xOut)
		VPADDQ(q, xOut, xOutQ)
		VPMINUQ(xOutQ, xOut, xOut)
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

	yOut := GP64()

	quo64 := GP64()
	MOVQ(y1S, reg.RDX)
	MULXQ(y0, y1S, quo64)
	IMULQ(q64, quo64)

	IMULQ(y1, y0)
	SUBQ(quo64, y0)

	y0Q := GP64()
	MOVQ(y0, y0Q)
	SUBQ(q64, y0Q)
	CMPQ(q64, y0)
	CMOVQLS(y0Q, y0)

	yOutQ := GP64()
	switch opType {
	case OpPure:
		yOut = y0
	case OpAdd:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		ADDQ(y0, yOut)
		MOVQ(yOut, yOutQ)
		SUBQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	case OpSub:
		MOVQ(Mem{Base: vOut, Index: i, Scale: 8}, yOut)
		SUBQ(y0, yOut)
		MOVQ(yOut, yOutQ)
		ADDQ(q64, yOutQ)
		CMPQ(q64, yOut)
		CMOVQLS(yOutQ, yOut)
	}

	MOVQ(yOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecReduceToAVX512() {
	TEXT("reduceToAVX512", NOSPLIT, "func(vOut, v []uint64, q, divHi uint64)")
	Pragma("noescape")

	maskLo := ZMM()
	VPBROADCASTQ(NewDataAddr(NewStaticSymbol("MASK_LO"), 0), maskLo)

	q64 := Load(Param("q"), GP64())
	q := ZMM()
	VPBROADCASTQ(NewParamAddr("q", 48), q)

	divHi64 := Load(Param("divHi"), GP64())
	divHi := ZMM()
	VPBROADCASTQ(NewParamAddr("divHi", 56), divHi)

	divHiHi := ZMM()
	VPSRLQ(Imm(32), divHi, divHiHi)

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

	xHi := ZMM()
	VPSRLQ(Imm(32), x, xHi)

	quo := ZMM()
	Mul64HiAVX512(x, xHi, divHi, divHiHi, maskLo, quo)
	VPMULLQ(quo, q, quo)
	VPSUBQ(quo, x, xOut)

	xOutQ := ZMM()
	VPSUBQ(q, xOut, xOutQ)
	VPMINUQ(xOutQ, xOut, xOut)

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
	MOVQ(divHi64, reg.RDX)
	MULXQ(y, yOut, quo64)
	IMULQ(q64, quo64)

	SUBQ(quo64, y)

	yQ := GP64()
	MOVQ(y, yQ)
	SUBQ(q64, yQ)
	CMPQ(q64, y)
	CMOVQLS(yQ, y)

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}

func VecReduceFixedQToAVX(factor int, avxType AVXType) {
	switch avxType {
	case TypeAVX2:
		switch factor {
		case 2:
			TEXT("reduce2QToAVX2", NOSPLIT, "func(vOut, v []uint64, q uint64)")
		case 4:
			TEXT("reduce4QToAVX2", NOSPLIT, "func(vOut, v []uint64, q uint64)")
		default:
			panic("factor should be 2 or 4")
		}
	case TypeAVX512:
		switch factor {
		case 2:
			TEXT("reduce2QToAVX512", NOSPLIT, "func(vOut, v []uint64, q uint64)")
		case 4:
			TEXT("reduce4QToAVX512", NOSPLIT, "func(vOut, v []uint64, q uint64)")
		default:
			panic("factor should be 2 or 4")
		}
	}

	Pragma("noescape")

	q64 := Load(Param("q"), GP64())
	q := VMM(avxType)
	VPBROADCASTQ(NewParamAddr("q", 48), q)

	twoQ64 := GP64()
	twoQ := VMM(avxType)
	if factor == 4 {
		MOVQ(q64, twoQ64)
		ADDQ(q64, twoQ64)
		VPADDQ(q, q, twoQ)
	}

	N := Load(Param("vOut").Len(), GP64())
	vOut := Load(Param("vOut").Base(), GP64())
	v := Load(Param("v").Base(), GP64())

	M := GP64()
	MOVQ(N, M)
	SHRQ(Imm(LogWidth(avxType)), M)
	SHLQ(Imm(LogWidth(avxType)), M)

	i := GP64()
	XORQ(i, i)
	JMP(LabelRef("loop_end"))
	Label("loop_body")

	x := VMM(avxType)
	VMOV(avxType, Mem{Base: v, Index: i, Scale: 8}, x)

	xOut, xOutQ := VMM(avxType), VMM(avxType)
	switch avxType {
	case TypeAVX2:
		switch factor {
		case 2:
			VPSUBQ(q, x, xOutQ)
			VBLENDVPD(xOutQ, x, xOutQ, xOut)
		case 4:
			VPSUBQ(twoQ, x, xOutQ)
			VBLENDVPD(xOutQ, x, xOutQ, xOut)
			VPSUBQ(q, xOut, xOutQ)
			VBLENDVPD(xOutQ, xOut, xOutQ, xOut)
		}

	case TypeAVX512:
		switch factor {
		case 2:
			VPSUBQ(q, x, xOutQ)
			VPMINUQ(xOutQ, x, xOut)
		case 4:
			VPSUBQ(twoQ, x, xOutQ)
			VPMINUQ(xOutQ, x, xOut)
			VPSUBQ(q, xOut, xOutQ)
			VPMINUQ(xOutQ, xOut, xOut)
		}
	}

	VMOV(avxType, xOut, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1<<LogWidth(avxType)), i)

	Label("loop_end")
	CMPQ(i, M)
	JL(LabelRef("loop_body"))

	JMP(LabelRef("leftover_loop_end"))
	Label("leftover_loop_body")

	y := GP64()
	MOVQ(Mem{Base: v, Index: i, Scale: 8}, y)

	yQ := GP64()
	if factor == 4 {
		MOVQ(y, yQ)
		SUBQ(twoQ64, yQ)
		CMPQ(twoQ64, y)
		CMOVQLS(yQ, y)
	}
	MOVQ(y, yQ)
	SUBQ(q64, yQ)
	CMPQ(q64, y)
	CMOVQLS(yQ, y)

	MOVQ(y, Mem{Base: vOut, Index: i, Scale: 8})

	ADDQ(Imm(1), i)

	Label("leftover_loop_end")
	CMPQ(i, N)
	JL(LabelRef("leftover_loop_body"))

	RET()
}
