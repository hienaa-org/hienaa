//go:generate go run . -vec -out ../../math/vec/asm_mod_amd64.s -stubs ../../math/vec/asm_mod_stub_amd64.go -pkg=vec
//go:generate go run . -ntt -out ../../math/dft/asm_ntt_pow2.s -stubs ../../math/dft/asm_ntt_pow2_stub_amd64.go -pkg=dft
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
	vec = flag.Bool("vec", false, "vec/asm_mod_amd64.s")
	ntt = flag.Bool("ntt", false, "rns/asm_ntt_pow2.s")
)

func main() {
	flag.Parse()

	Constraint(buildtags.Term("amd64"))
	Constraint(buildtags.Not("purego"))

	if *vec {
		AddVecToAVX2(false)
		AddVecToAVX2(true)

		AddVecToAVX512(false)
		AddVecToAVX512(true)

		ScalarAddVecToAVX2(false)
		ScalarAddVecToAVX2(true)

		ScalarAddVecToAVX512(false)
		ScalarAddVecToAVX512(true)

		SubVecToAVX2(false)
		SubVecToAVX2(true)

		SubVecToAVX512(false)
		SubVecToAVX512(true)

		ScalarSubVecToAVX2(false)
		ScalarSubVecToAVX2(true)

		ScalarSubVecToAVX512(false)
		ScalarSubVecToAVX512(true)
	}

	if *ntt {
		NTTConstants()

		NTTInPlacePow2UnrollAVX2()
		NTTInPlacePow2UnrollAVX512()

		InvNTTInPlacePow2UnrollAVX2()
		InvNTTInPlacePow2UnrollAVX512()
	}

	Generate()
}
