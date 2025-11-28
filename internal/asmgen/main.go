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
	vec = flag.Bool("vec", false, "asm_mod_amd64.s")
	ntt = flag.Bool("ntt", false, "asm_ntt_pow2.s")
)

func main() {
	flag.Parse()

	Constraint(buildtags.Term("amd64"))
	Constraint(buildtags.Not("purego"))

	if *vec {
		VecConstants()

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

		MulWordVecToAVX512(OpMul)
		MulWordVecToAVX512(OpMulAdd)
		MulWordVecToAVX512(OpMulSub)

		MMulVecToAVX512(OpMul, false)
		MMulVecToAVX512(OpMulAdd, false)
		MMulVecToAVX512(OpMulSub, false)

		MMulVecToAVX512(OpMul, true)
		MMulVecToAVX512(OpMulAdd, true)
		MMulVecToAVX512(OpMulSub, true)
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
