//go:generate go run . -vec -out ../../math/mod/asm_vec_ops_amd64.s -stubs ../../math/mod/asm_vec_ops_stub_amd64.go -pkg=mod
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
	vec = flag.Bool("vec", false, "asm_vec_ops_amd64.s")
)

func main() {
	flag.Parse()

	Constraint(buildtags.Term("amd64"))
	Constraint(buildtags.Not("purego"))

	if *vec {
		VecOpConstants()

		AddVecToAVX2(false)
		AddVecToAVX2(true)

		SubVecToAVX2(false)
		SubVecToAVX2(true)

		ScalarMulVecToX86(false, Mul)
		ScalarMulVecToX86(false, MulAdd)
		ScalarMulVecToX86(false, MulSub)

		ScalarMulVecToX86(true, Mul)
		ScalarMulVecToX86(true, MulAdd)
		ScalarMulVecToX86(true, MulSub)

		MulVecToX86(false, Mul)
		MulVecToX86(false, MulAdd)
		MulVecToX86(false, MulSub)

		MulVecToX86(true, Mul)
		MulVecToX86(true, MulAdd)
		MulVecToX86(true, MulSub)
	}

	Generate()
}
