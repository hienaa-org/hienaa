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

		VecAddSubToAVX2(OpAdd, false)
		VecAddSubToAVX2(OpAdd, true)
		VecAddSubToAVX2(OpSub, false)
		VecAddSubToAVX2(OpSub, true)

		VecAddSubToAVX512(OpAdd, false)
		VecAddSubToAVX512(OpAdd, true)
		VecAddSubToAVX512(OpSub, false)
		VecAddSubToAVX512(OpSub, true)

		VecAddSubScalarToAVX2(OpAdd, false)
		VecAddSubScalarToAVX2(OpAdd, true)
		VecAddSubScalarToAVX2(OpSub, false)
		VecAddSubScalarToAVX2(OpSub, true)

		VecAddSubScalarToAVX512(OpAdd, false)
		VecAddSubScalarToAVX512(OpAdd, true)
		VecAddSubScalarToAVX512(OpSub, false)
		VecAddSubScalarToAVX512(OpSub, true)

		VecNegToAVX2(false)
		VecNegToAVX2(true)

		VecNegToAVX512(false)
		VecNegToAVX512(true)

		VecMFormToAVX512()
		VecInvMFormToAVX512()

		VecMulScalarWordToAVX2(OpPure)
		VecMulScalarWordToAVX2(OpAdd)
		VecMulScalarWordToAVX2(OpSub)

		VecMulScalarWordToAVX512(OpPure)
		VecMulScalarWordToAVX512(OpAdd)
		VecMulScalarWordToAVX512(OpSub)

		VecSMulScalarToAVX512(OpPure, false)
		VecSMulScalarToAVX512(OpAdd, false)
		VecSMulScalarToAVX512(OpSub, false)

		VecSMulScalarToAVX512(OpPure, true)
		VecSMulScalarToAVX512(OpAdd, true)
		VecSMulScalarToAVX512(OpSub, true)

		VecMulWordToAVX2(OpPure)
		VecMulWordToAVX2(OpAdd)
		VecMulWordToAVX2(OpSub)

		VecMulWordToAVX512(OpPure)
		VecMulWordToAVX512(OpAdd)
		VecMulWordToAVX512(OpSub)

		VecSMulToAVX512(OpPure, false)
		VecSMulToAVX512(OpAdd, false)
		VecSMulToAVX512(OpSub, false)

		VecSMulToAVX512(OpPure, true)
		VecSMulToAVX512(OpAdd, true)
		VecSMulToAVX512(OpSub, true)

		VecReduceToAVX512()
	}

	if *ntt {
		NTTConstants()

		FwdNTTInPlacePow2UnrollAVX512()
		InvNTTInPlacePow2UnrollAVX512()
	}

	Generate()
}
