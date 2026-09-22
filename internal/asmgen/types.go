package main

import (
	. "github.com/mmcloughlin/avo/build"
	"github.com/mmcloughlin/avo/operand"
	"github.com/mmcloughlin/avo/reg"
)

type OpType int

const (
	OpPure OpType = iota
	OpAdd
	OpSub
)

type AVXType int

const (
	TypeAVX2 AVXType = iota
	TypeAVX512
)

func LogWidth(avxType AVXType) uint64 {
	switch avxType {
	case TypeAVX2:
		return 2
	case TypeAVX512:
		return 3
	}
	panic("avxType unavailable")
}

func VMM(avxType AVXType) reg.VecVirtual {
	switch avxType {
	case TypeAVX2:
		return YMM()
	case TypeAVX512:
		return ZMM()
	}
	panic("avxType unavailable")
}

func VMOV(avxType AVXType, mxy, mxy1 operand.Op) {
	switch avxType {
	case TypeAVX2:
		VMOVDQU(mxy, mxy1)
		return
	case TypeAVX512:
		VMOVDQU64(mxy, mxy1)
		return
	}
	panic("avxType unavailable")
}
