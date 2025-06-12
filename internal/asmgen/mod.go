package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func BMod64(x, q, divHi, xOut reg.Register) {
	MOVQ(x, reg.RDX)
	MULXQ(divHi, reg.RDX, xOut) // xOut = tmp

	IMULQ(q, reg.RDX)
	MOVQ(x, xOut)
	SUBQ(reg.RDX, xOut)

	xOutSubQ := GP64()
	MOVQ(xOut, xOutSubQ)
	SUBQ(q, xOutSubQ)
	CMPQ(xOut, q)
	CMOVQGE(xOutSubQ, xOut)
}

func BMod128(xHi, xLo, q, divHi, divLo, xOut reg.Register) {
	BMod128Lazy(xHi, xLo, q, divHi, divLo, xOut)

	xOutSubQ := GP64()
	MOVQ(xOut, xOutSubQ)
	SUBQ(q, xOutSubQ)
	CMPQ(xOut, q)
	CMOVQGE(xOutSubQ, xOut)
}

func BMod128Lazy(xHi, xLo, q, divHi, divLo, xOut reg.Register) {
	// quo := xHi * divHi
	quo := GP64()
	MOVQ(xHi, quo)
	IMULQ(divHi, quo)

	// quoMid0, quoMid0Lo := xLo * divHi
	// quo += quoMid0
	quoMid0, quoMid0Lo := GP64(), GP64()
	MOVQ(xLo, reg.RDX)
	MULXQ(divHi, quoMid0Lo, quoMid0)
	ADDQ(quoMid0, quo)

	// quoLo, _ := xLo * divLo
	quoLo := GP64()
	MULXQ(divLo, xOut, quoLo) // xOut = tmp

	// quoMid1, quoMid1Lo := xHi * divLo
	// quo += quoMid1
	quoMid1, quoMid1Lo := GP64(), GP64()
	MOVQ(xHi, reg.RDX)
	MULXQ(divLo, quoMid1Lo, quoMid1)
	ADDQ(quoMid1, quo)

	// quoMidSum = quoMid0Lo + quoLo
	// quo += carry
	ADDQ(quoMid0Lo, quoLo)
	ADCQ(Imm(0), quo)

	// quoMidSum = quoMid1Lo + quoMidSum
	// quo += carry
	ADDQ(quoMid1Lo, quoLo)
	ADCQ(Imm(0), quo)

	IMULQ(q, quo)
	MOVQ(xLo, xOut)
	SUBQ(quo, xOut)
}

func MMul(x0M, x1M, q, inv, xOutM reg.Register) {
	MMulLazy(x0M, x1M, q, inv, xOutM)

	xOutSubQ := GP64()
	MOVQ(xOutM, xOutSubQ)
	SUBQ(q, xOutSubQ)
	CMPQ(xOutM, q)
	CMOVQGE(xOutSubQ, xOutM)
}

func MMulLazy(x0M, x1M, q, inv, xOutM reg.Register) {
	MOVQ(x0M, reg.RDX)
	MULXQ(x1M, reg.RDX, xOutM)

	tmp := GP64()
	IMULQ(inv, reg.RDX)
	MULXQ(q, tmp, reg.RDX)

	SUBQ(reg.RDX, xOutM)
	ADDQ(q, xOutM)
}

func BMod128AVX2(xHi, xLo, q, divHi, divLo, xOut reg.VecVirtual) {
	BMod128LazyAVX2(xHi, xLo, q, divHi, divLo, xOut)

	subQ := YMM()
	GreaterOrEqualThanAVX2(xOut, q, subQ)
	VPAND(q, subQ, subQ)
	VPSUBQ(subQ, xOut, xOut)
}

func BMod128LazyAVX2(xHi, xLo, q, divHi, divLo, xOut reg.VecVirtual) {
	// quo := xHi * divHi
	quo := YMM()
	Mul64LoAVX2(xHi, divHi, quo)

	// quoMid0, quoMid0Lo := xLo * divHi
	// quo += quoMid0
	quoMid0, quoMid0Lo := YMM(), YMM()
	Mul64AVX2(xLo, divHi, quoMid0, quoMid0Lo)
	VPADDQ(quoMid0, quo, quo)

	// quoMid1, quoMid1Lo := xHi * divLo
	// quo += quoMid1
	quoMid1, quoMid1Lo := YMM(), YMM()
	Mul64AVX2(xHi, divLo, quoMid1, quoMid1Lo)
	VPADDQ(quoMid1, quo, quo)

	// quoLo, _ := xLo * divLo
	quoLo := YMM()
	Mul64HiAVX2(xLo, divLo, quoLo)

	// quoMidSum = quoMid0Lo + quoLo
	// quo += carry
	quoMidSum, quoMidCarry := YMM(), YMM()
	Add64AVX2(quoMid0Lo, quoLo, quoMidSum, quoMidCarry)
	VPADDQ(quoMidCarry, quo, quo)

	// quoMidSum = quoMid1Lo + quoMidSum
	// quo += carry
	Add64AVX2(quoMid1Lo, quoMidSum, quoMidSum, quoMidCarry)
	VPADDQ(quoMidCarry, quo, quo)

	Mul64LoAVX2(q, quo, quo)
	VPSUBQ(quo, xLo, xOut)
}

func MMulAVX2(x0M, x1M, q, inv, xOutM reg.VecVirtual) {
	MMulLazyAVX2(x0M, x1M, q, inv, xOutM)

	subQ := YMM()
	GreaterOrEqualThanAVX2(xOutM, q, subQ)
	VPAND(q, subQ, subQ)
	VPSUBQ(subQ, xOutM, xOutM)
}

func MMulLazyAVX2(x0M, x1M, q, inv, xOutM reg.VecVirtual) {
	xOutMHi, xOutMLo := YMM(), YMM()
	Mul64AVX2(x0M, x1M, xOutMLo, xOutMHi)

	Mul64LoAVX2(xOutMLo, inv, xOutMLo)
	Mul64HiAVX2(xOutMLo, q, xOutMLo)

	VPSUBQ(xOutMLo, xOutMHi, xOutM)
	VPADDQ(q, xOutMHi, xOutMHi)
}
