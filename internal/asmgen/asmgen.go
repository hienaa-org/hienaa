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
	OpPure OpType = iota
	OpAdd
	OpSub
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

		AddSubVecToAVX2(OpAdd, false)
		AddSubVecToAVX2(OpAdd, true)
		AddSubVecToAVX2(OpSub, false)
		AddSubVecToAVX2(OpSub, true)

		AddSubVecToAVX512(OpAdd, false)
		AddSubVecToAVX512(OpAdd, true)
		AddSubVecToAVX512(OpSub, false)
		AddSubVecToAVX512(OpSub, true)

		ScalarAddSubVecToAVX2(OpAdd, false)
		ScalarAddSubVecToAVX2(OpAdd, true)
		ScalarAddSubVecToAVX2(OpSub, false)
		ScalarAddSubVecToAVX2(OpSub, true)

		ScalarAddSubVecToAVX512(OpAdd, false)
		ScalarAddSubVecToAVX512(OpAdd, true)
		ScalarAddSubVecToAVX512(OpSub, false)
		ScalarAddSubVecToAVX512(OpSub, true)

		NegVecToAVX2(false)
		NegVecToAVX2(true)

		NegVecToAVX512(false)
		NegVecToAVX512(true)

		MFormVecToAVX512()
		InvMFormVecToAVX512()

		ScalarMulWordVecToAVX2(OpPure)
		ScalarMulWordVecToAVX2(OpAdd)
		ScalarMulWordVecToAVX2(OpSub)

		ScalarMulWordVecToAVX512(OpPure)
		ScalarMulWordVecToAVX512(OpAdd)
		ScalarMulWordVecToAVX512(OpSub)

		ScalarMulVecToAVX512(OpPure, false)
		ScalarMulVecToAVX512(OpAdd, false)
		ScalarMulVecToAVX512(OpSub, false)

		ScalarMulVecToAVX512(OpPure, true)
		ScalarMulVecToAVX512(OpAdd, true)
		ScalarMulVecToAVX512(OpSub, true)

		ScalarMMulVecToAVX512(OpPure, false)
		ScalarMMulVecToAVX512(OpAdd, false)
		ScalarMMulVecToAVX512(OpSub, false)

		ScalarMMulVecToAVX512(OpPure, true)
		ScalarMMulVecToAVX512(OpAdd, true)
		ScalarMMulVecToAVX512(OpSub, true)

		MulWordVecToAVX2(OpPure)
		MulWordVecToAVX2(OpAdd)
		MulWordVecToAVX2(OpSub)

		MulWordVecToAVX512(OpPure)
		MulWordVecToAVX512(OpAdd)
		MulWordVecToAVX512(OpSub)

		MMulVecToAVX512(OpPure, false)
		MMulVecToAVX512(OpAdd, false)
		MMulVecToAVX512(OpSub, false)

		MMulVecToAVX512(OpPure, true)
		MMulVecToAVX512(OpAdd, true)
		MMulVecToAVX512(OpSub, true)

		SMulVecToAVX512(OpPure, false)
		SMulVecToAVX512(OpAdd, false)
		SMulVecToAVX512(OpSub, false)

		SMulVecToAVX512(OpPure, true)
		SMulVecToAVX512(OpAdd, true)
		SMulVecToAVX512(OpSub, true)
	}

	if *ntt {
		NTTConstants()

		FwdNTTInPlacePow2UnrollAVX2()
		FwdNTTInPlacePow2UnrollAVX512()

		InvNTTInPlacePow2UnrollAVX2()
		InvNTTInPlacePow2UnrollAVX512()
	}

	Generate()
}
