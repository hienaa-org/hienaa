package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func BMod(xHi, xLo, q, divHi, divLo, rem reg.Register) {
	BModLazy(xHi, xLo, q, divHi, divLo, rem)

	leftOver := GP64()
	XORQ(leftOver, leftOver)

	CMPQ(rem, q)
	CMOVQGE(q, leftOver)
	SUBQ(leftOver, rem)
}

func BModLazy(xHi, xLo, q, divHi, divLo, rem reg.Register) {
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
	MOVQ(xLo, rem)
	SUBQ(quo, rem)
}
