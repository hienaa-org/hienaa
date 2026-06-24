package main

import (
	. "github.com/mmcloughlin/avo/build"
	. "github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

func ShuffleForLoAVX2(x, xSwap reg.VecVirtual) {
	VPSHUFD(Imm(0b10_11_00_01), x, xSwap)
}

func Mul64LoAVX2(x0, x1, x1Swap, maskHi, xOut reg.VecVirtual) {
	xLoLo := YMM()
	VPMULUDQ(x1, x0, xLoLo)

	xMid0, xMid1 := YMM(), YMM()
	VPMULLD(x1Swap, x0, xMid0)
	VPSLLQ(Imm(32), xMid0, xMid1)
	VPAND(maskHi, xMid0, xOut)

	VPADDQ(xMid1, xOut, xOut)
	VPADDQ(xLoLo, xOut, xOut)
}

func Mul64HiAVX2(x0, x0Hi, x1, x1Hi, maskLo, xOut reg.VecVirtual) {
	xLoLo, xLoHi, xHiLo := YMM(), YMM(), YMM()
	VPMULUDQ(x1, x0, xLoLo)
	VPMULUDQ(x1Hi, x0, xLoHi)
	VPMULUDQ(x1, x0Hi, xHiLo)
	VPMULUDQ(x1Hi, x0Hi, xOut)

	VPSRLQ(Imm(32), xLoLo, xLoLo)

	xMidHi, xMidLo := YMM(), YMM()
	VPADDQ(xLoLo, xLoHi, xMidHi)
	VPAND(maskLo, xMidHi, xMidLo)
	VPSRLQ(Imm(32), xMidHi, xMidHi)
	VPADDQ(xOut, xMidHi, xOut)

	VPADDQ(xMidLo, xHiLo, xMidHi)
	VPSRLQ(Imm(32), xMidHi, xMidHi)
	VPADDQ(xOut, xMidHi, xOut)
}

func Mul64HiAVX512(x0, x0Hi, x1, x1Hi, maskLo, xOut reg.VecVirtual) {
	xLoLo, xLoHi, xHiLo := ZMM(), ZMM(), ZMM()
	VPMULUDQ(x1, x0, xLoLo)
	VPMULUDQ(x1Hi, x0, xLoHi)
	VPMULUDQ(x1, x0Hi, xHiLo)
	VPMULUDQ(x1Hi, x0Hi, xOut)

	VPSRLQ(Imm(32), xLoLo, xLoLo)

	xMidHi, xMidLo := ZMM(), ZMM()
	VPADDQ(xLoLo, xLoHi, xMidHi)
	VPANDQ(maskLo, xMidHi, xMidLo)
	VPSRLQ(Imm(32), xMidHi, xMidHi)
	VPADDQ(xOut, xMidHi, xOut)

	VPADDQ(xMidLo, xHiLo, xMidHi)
	VPSRLQ(Imm(32), xMidHi, xMidHi)
	VPADDQ(xOut, xMidHi, xOut)
}

func Mul64HiApproxAVX512(x0, x0Hi, x1, x1Hi, maskLo, xOut reg.VecVirtual) {
	xLoHi, xHiLo := ZMM(), ZMM()
	VPMULUDQ(x1Hi, x0, xLoHi)
	VPMULUDQ(x1, x0Hi, xHiLo)
	VPMULUDQ(x1Hi, x0Hi, xOut)

	VPSRLQ(Imm(32), xLoHi, xLoHi)
	VPSRLQ(Imm(32), xHiLo, xHiLo)
	VPADDQ(xOut, xLoHi, xOut)
	VPADDQ(xOut, xHiLo, xOut)
}
