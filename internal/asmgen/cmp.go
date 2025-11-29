package main

import (
	. "github.com/mmcloughlin/avo/build"
	"github.com/mmcloughlin/avo/reg"
)

func EqualAVX2(x, y, cmpOut reg.VecVirtual) {
	VPCMPEQQ(x, y, cmpOut)
}

func LessThanAVX2(x, y, maskSign, cmpOut reg.VecVirtual) {
	xSigned, ySigned := YMM(), YMM()
	VPXOR(maskSign, x, xSigned)
	VPXOR(maskSign, y, ySigned)
	VPCMPGTQ(xSigned, ySigned, cmpOut)
}

func LessOrEqualThanAVX2(x, y, maskSign, cmpOut reg.VecVirtual) {
	GreaterThanAVX2(y, x, maskSign, cmpOut)
}

func GreaterThanAVX2(x, y, maskSign, cmpOut reg.VecVirtual) {
	xSigned, ySigned := YMM(), YMM()
	VPXOR(maskSign, x, xSigned)
	VPXOR(maskSign, y, ySigned)
	VPCMPGTQ(ySigned, xSigned, cmpOut)
}

func GreaterOrEqualThanAVX2(x, y, maskSign, cmpOut reg.VecVirtual) {
	LessThanAVX2(y, x, maskSign, cmpOut)
}
