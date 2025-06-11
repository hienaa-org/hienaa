package main

import (
	. "github.com/mmcloughlin/avo/build"
	"github.com/mmcloughlin/avo/reg"
)

func EqualAVX2(x, y, cmpOut reg.VecVirtual) {
	VPCMPEQQ(x, y, cmpOut)
}

func LessThanAVX2(x, y, cmpOut reg.VecVirtual) {
	VPCMPGTQ(x, y, cmpOut)
}

func LessOrEqualThanAVX2(x, y, cmpOut reg.VecVirtual) {
	ltOut, eqOut := YMM(), YMM()
	LessThanAVX2(x, y, ltOut)
	EqualAVX2(x, y, eqOut)
	VPORQ(ltOut, eqOut, cmpOut)
}

func GreaterThanAVX2(x, y, cmpOut reg.VecVirtual) {
	VPCMPGTQ(y, x, cmpOut)
}

func GreaterOrEqualThanAVX2(x, y, cmpOut reg.VecVirtual) {
	gtOut, eqOut := YMM(), YMM()
	GreaterThanAVX2(x, y, gtOut)
	EqualAVX2(x, y, eqOut)
	VPORQ(gtOut, eqOut, cmpOut)
}
