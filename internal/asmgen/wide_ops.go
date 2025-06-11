package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func Add64AVX2(x0, x1, xOut, carry reg.VecVirtual) {
	xAnd, xOr := YMM(), YMM()
	VPAND(x0, x1, xAnd)
	VPOR(x0, x1, xOr)

	VPADDQ(x0, x1, xOut)

	VPANDN(xOr, xOut, carry)
	VPOR(xAnd, carry, carry)
	VPSRLQ(Imm(63), carry, carry)
}

func Mul64LoAVX2(x0, x1, xOut reg.VecVirtual) {
	x0Hi, x0Lo := YMM(), x0
	x1Lo := x1

	VPSHUFD(Imm(0b10110001), x0, x0Hi)

	xOutMid0, xOutMid1 := YMM(), YMM()
	VPMULLD(x0Hi, x1Lo, xOutMid0)
	VPSRLQ(Imm(32), xOutMid0, xOutMid1)

	VPADDQ(xOutMid1, xOutMid0, xOutMid0)
	VPSLLQ(Imm(32), xOutMid0, xOutMid0)

	VPMULUDQ(x0Lo, x1Lo, xOut)
	VPADDQ(xOutMid0, xOut, xOut)
}

func Mul64HiAVX2(x0, x1, xOut reg.VecVirtual) {
	x0Hi, x0Lo := YMM(), x0
	x1Hi, x1Lo := YMM(), x1

	VPSHUFD(Imm(0b10110001), x0, x0Hi)
	VPSHUFD(Imm(0b10110001), x1, x1Hi)

	xOutMid0, xOutMid1, xOutLo := YMM(), YMM(), YMM()
	VPMULUDQ(x0Lo, x1Lo, xOutLo)
	VPMULUDQ(x0Hi, x1Lo, xOutMid0)
	VPMULUDQ(x0Lo, x1Hi, xOutMid1)
	VPMULUDQ(x0Hi, x1Hi, xOut)

	xOutMidHi, xOutMidLo := YMM(), YMM()
	VPSRLQ(Imm(32), xOutLo, xOutMidLo)
	VPADDQ(xOutMid0, xOutMidLo, xOutMidLo)

	VPSRLQ(Imm(32), xOutMidLo, xOutMidHi)
	VPSLLQ(Imm(32), xOutMidLo, xOutMidLo)
	VPSRLQ(Imm(32), xOutMidLo, xOutMidLo)

	VPADDQ(xOutMid1, xOutMidLo, xOutMidLo)
	VPSRLQ(Imm(32), xOutMidLo, xOutMidLo)

	VPADDQ(xOutMidHi, xOut, xOut)
	VPADDQ(xOutMidLo, xOut, xOut)
}

func Mul64AVX2(x0, x1, xOutHi, xOutLo reg.VecVirtual) {
	x0Hi, x0Lo := YMM(), x0
	x1Hi, x1Lo := YMM(), x1

	VPSHUFD(Imm(0b10110001), x0, x0Hi)
	VPSHUFD(Imm(0b10110001), x1, x1Hi)

	xOutMid0, xOutMid1 := YMM(), YMM()
	VPMULUDQ(x0Lo, x1Lo, xOutLo)
	VPMULUDQ(x0Hi, x1Lo, xOutMid0)
	VPMULUDQ(x0Lo, x1Hi, xOutMid1)
	VPMULUDQ(x0Hi, x1Hi, xOutHi)

	xOutMid := YMM()
	VPADDQ(xOutMid0, xOutMid1, xOutMid)
	VPSLLQ(Imm(32), xOutMid, xOutMid)

	xOutMidHi, xOutMidLo := YMM(), YMM()
	VPSRLQ(Imm(32), xOutLo, xOutMidLo)
	VPADDQ(xOutMid0, xOutMidLo, xOutMidLo)

	VPSRLQ(Imm(32), xOutMidLo, xOutMidHi)
	VPSLLQ(Imm(32), xOutMidLo, xOutMidLo)
	VPSRLQ(Imm(32), xOutMidLo, xOutMidLo)

	VPADDQ(xOutMid1, xOutMidLo, xOutMidLo)
	VPSRLQ(Imm(32), xOutMidLo, xOutMidLo)

	VPADDQ(xOutMid, xOutLo, xOutLo)

	VPADDQ(xOutMidHi, xOutHi, xOutHi)
	VPADDQ(xOutMidLo, xOutHi, xOutHi)
}
