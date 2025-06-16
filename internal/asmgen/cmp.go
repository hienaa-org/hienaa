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

func LessOrEqualThanAVX2(x, y, allOne, cmpOut reg.VecVirtual) {
	GreaterThanAVX2(x, y, cmpOut)
	VPXOR(allOne, cmpOut, cmpOut)
}

func GreaterThanAVX2(x, y, cmpOut reg.VecVirtual) {
	VPCMPGTQ(y, x, cmpOut)
}

func GreaterOrEqualThanAVX2(x, y, allOne, cmpOut reg.VecVirtual) {
	LessThanAVX2(x, y, cmpOut)
	VPXOR(allOne, cmpOut, cmpOut)
}
