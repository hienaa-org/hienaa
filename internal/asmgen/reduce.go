package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func BMod(xHi, xLo, q, divHi, divLo, xOut reg.Register) {
	BModLazy(xHi, xLo, q, divHi, divLo, xOut)

	xOutSubQ := GP64()
	MOVQ(xOut, xOutSubQ)
	SUBQ(q, xOutSubQ)
	CMPQ(xOut, q)
	CMOVQGE(xOutSubQ, xOut)
}

func BModLazy(xHi, xLo, q, divHi, divLo, xOut reg.Register) {
	// quo := xHi * divHi
	quo := GP64()
	MOVQ(xHi, quo)
	IMULQ(divHi, quo)

	// quoMid0, quoMid0Lo := xLo * divHi
	// quo += quoMid0
	MOVQ(xLo, reg.RAX)
	MULQ(divHi)
	ADDQ(reg.RDX, quo)

	quoMid0Lo := GP64()
	MOVQ(reg.RAX, quoMid0Lo)

	// quoMid1, quoMid1Lo := xHi * divLo
	// quo += quoMid1
	MOVQ(xHi, reg.RAX)
	MULQ(divLo)
	ADDQ(reg.RDX, quo)

	quoMid1Lo := GP64()
	MOVQ(reg.RAX, quoMid1Lo)

	// quoLo, _ := xLo * divLo
	// RDX = quoLo
	MOVQ(xLo, reg.RAX)
	MULQ(divLo)

	// quoMidSum = quoMid0Lo + quoLo
	// quo += carry
	ADDQ(quoMid0Lo, reg.RDX)
	ADCQ(Imm(0), quo)

	// quoMidSum = quoMid1Lo + quoMidSum
	// quo += carry
	ADDQ(quoMid1Lo, reg.RDX)
	ADCQ(Imm(0), quo)

	IMULQ(q, quo)
	MOVQ(xLo, xOut)
	SUBQ(quo, xOut)
}

func MMul(xHi, xLo, q, inv, xOut reg.Register) {
	MMulLazy(xHi, xLo, q, inv, xOut)

	xOutSubQ := GP64()
	MOVQ(xOut, xOutSubQ)
	SUBQ(q, xOutSubQ)
	CMPQ(xOut, q)
	CMOVQGE(xOutSubQ, xOut)
}

func MMulLazy(xHi, xLo, q, inv, xOut reg.Register) {
	xOutMHi := GP64()
	MOVQ(xLo, reg.RAX)
	MULQ(xHi)
	MOVQ(reg.RDX, xOutMHi)

	// RAX = xOutMLo
	MULQ(inv)
	MULQ(q)

	// RDX = wHi
	SUBQ(reg.RDX, xOutMHi)
	ADDQ(q, xOutMHi)

	MOVQ(xOutMHi, xOut)
}

func BModAVX2(xHi, xLo, q, divHi, divLo, xOut reg.VecVirtual) {
	BModLazyAVX2(xHi, xLo, q, divHi, divLo, xOut)

	subQ := YMM()
	GreaterOrEqualThanAVX2(xOut, q, subQ)
	VPAND(q, subQ, subQ)
	VPSUBQ(subQ, xOut, xOut)
}

func BModLazyAVX2(xHi, xLo, q, divHi, divLo, xOut reg.VecVirtual) {
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
