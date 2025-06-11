//go:generate go run . -vec -out ../../math/vec/asm_ops_amd64.s -stubs ../../math/vec/asm_ops_stub_amd64.go -pkg=vec
package main

import (
	"flag"

	. "github.com/mmcloughlin/avo/build"
	"github.com/mmcloughlin/avo/buildtags"
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

		VecBMulToAVX2(false)
		VecBMulToAVX2(true)
	}

	Generate()
}
