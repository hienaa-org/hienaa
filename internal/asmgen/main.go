//go:generate go run . -vec -out ../../math/vec/asm_ops_amd64.s -stubs ../../math/vec/asm_ops_stub_amd64.go -pkg=vec
package main

import (
	"flag"

	. "github.com/mmcloughlin/avo/build"
	"github.com/mmcloughlin/avo/buildtags"
)

type OpType int

const (
	Mul OpType = iota
	MulAdd
	MulSub
)

var (
	vec = flag.Bool("vec", false, "asm_vec_amd64.s")
)

func main() {
	flag.Parse()

	Constraint(buildtags.Term("amd64"))
	Constraint(buildtags.Not("purego"))

	if *vec {
		VecOpConstants()

		VecAddToAVX2(false)
		VecAddToAVX2(true)

		VecSubToAVX2(false)
		VecSubToAVX2(true)

		VecBMulToX86(false, Mul)
		VecBMulToX86(false, MulAdd)
		VecBMulToX86(false, MulSub)

		VecBMulToX86(true, Mul)
		VecBMulToX86(true, MulAdd)
		VecBMulToX86(true, MulSub)

		VecMMulToX86(false, Mul)
		VecMMulToX86(false, MulAdd)
		VecMMulToX86(false, MulSub)

		VecMMulToX86(true, Mul)
		VecMMulToX86(true, MulAdd)
		VecMMulToX86(true, MulSub)
	}

	Generate()
}
