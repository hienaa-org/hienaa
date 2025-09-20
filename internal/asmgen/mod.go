package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func BMod128(xHi, xLo, q, divHi, divLo, xOut reg.Register) {
	BMod128Lazy(xHi, xLo, q, divHi, divLo, xOut)

	xOutSubQ := GP64()
	MOVQ(xOut, xOutSubQ)
	SUBQ(q, xOutSubQ)
	CMPQ(xOut, q)
	CMOVQCC(xOutSubQ, xOut)
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
