package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func Add64AVX2(x0, x1, xOut, carry reg.VecVirtual) {
	// Port of [bits.Add64]

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

	x0Hix1Lo, x0Lox1Hi := YMM(), YMM()
	VPMULLD(x0Hi, x1Lo, x0Hix1Lo)
	VPSRLQ(Imm(32), x0Hix1Lo, x0Lox1Hi)

	xOutMid := YMM()
	VPADDQ(x0Hix1Lo, x0Lox1Hi, xOutMid)
	VPSLLQ(Imm(32), xOutMid, xOutMid)

	VPMULUDQ(x0Lo, x1Lo, xOut)
	VPADDQ(xOutMid, xOut, xOut)
}

func Mul64HiAVX2(x0, x1, xOut reg.VecVirtual) {
	// Port of [bits.Mul64]

	x0Hi, x0Lo := YMM(), x0
	x1Hi, x1Lo := YMM(), x1

	VPSHUFD(Imm(0b10110001), x0, x0Hi)
	VPSHUFD(Imm(0b10110001), x1, x1Hi)

	w0 := YMM()
	VPMULUDQ(x0Lo, x1Lo, w0)
	VPSRLQ(Imm(32), w0, w0)

	t := YMM()
	VPMULUDQ(x0Hi, x1Lo, t)
	VPADDQ(w0, t, t)

	w1 := YMM()
	VPSLLQ(Imm(32), t, w1)
	VPSRLQ(Imm(32), w1, w1)

	w2 := YMM()
	VPSRLQ(Imm(32), t, w2)

	VPMULUDQ(x0Lo, x1Hi, w0)
	VPADDQ(w0, w1, w1)
	VPSRLQ(Imm(32), w1, w1)

	VPMULUDQ(x0Hi, x1Hi, xOut)
	VPADDQ(w2, xOut, xOut)
	VPADDQ(w1, xOut, xOut)
}

func Mul64AVX2(x0, x1, xOutHi, xOutLo reg.VecVirtual) {
	x0Hi, x0Lo := YMM(), x0
	x1Hi, x1Lo := YMM(), x1

	VPSHUFD(Imm(0b10110001), x0, x0Hi)
	VPSHUFD(Imm(0b10110001), x1, x1Hi)

	w0 := YMM()
	VPMULUDQ(x0Lo, x1Lo, xOutLo)
	VPSRLQ(Imm(32), xOutLo, w0)

	t := YMM()
	VPMULUDQ(x0Hi, x1Lo, t)

	xOutCross := YMM()
	VPSLLQ(Imm(32), t, xOutCross)
	VPADDQ(xOutCross, xOutLo, xOutLo)

	VPADDQ(w0, t, t)

	w1 := YMM()
	VPSLLQ(Imm(32), t, w1)
	VPSRLQ(Imm(32), w1, w1)

	w2 := YMM()
	VPSRLQ(Imm(32), t, w2)

	VPMULUDQ(x0Lo, x1Hi, w0)
	VPSLLQ(Imm(32), w0, xOutCross)
	VPADDQ(xOutCross, xOutLo, xOutLo)

	VPADDQ(w0, w1, w1)
	VPSRLQ(Imm(32), w1, w1)

	VPMULUDQ(x0Hi, x1Hi, xOutHi)
	VPADDQ(w2, xOutHi, xOutHi)
	VPADDQ(w1, xOutHi, xOutHi)
}
