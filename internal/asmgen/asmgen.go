//go:generate go run . -vec -out ../../math/vec/asm_mod_amd64.s -stubs ../../math/vec/asm_mod_stub_amd64.go -pkg=vec
package main

// go:generate go run . -ntt -out ../../math/dft/asm_ntt_pow2.s -stubs ../../math/dft/asm_ntt_pow2_stub_amd64.go -pkg=dft

import (
	"flag"

	. "github.com/mmcloughlin/avo/build"
	"github.com/mmcloughlin/avo/buildtags"
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

		VecAddSubToAVX(TypeAVX2, OpAdd, false)
		VecAddSubToAVX(TypeAVX2, OpAdd, true)
		VecAddSubToAVX(TypeAVX2, OpSub, false)
		VecAddSubToAVX(TypeAVX2, OpSub, true)

		VecAddSubToAVX(TypeAVX512, OpAdd, false)
		VecAddSubToAVX(TypeAVX512, OpAdd, true)
		VecAddSubToAVX(TypeAVX512, OpSub, false)
		VecAddSubToAVX(TypeAVX512, OpSub, true)

		VecAddSubScalarToAVX(TypeAVX2, OpAdd, false, false)
		VecAddSubScalarToAVX(TypeAVX2, OpAdd, true, false)
		VecAddSubScalarToAVX(TypeAVX2, OpSub, false, false)
		VecAddSubScalarToAVX(TypeAVX2, OpSub, true, false)
		VecAddSubScalarToAVX(TypeAVX2, OpSub, false, true)
		VecAddSubScalarToAVX(TypeAVX2, OpSub, true, true)

		VecAddSubScalarToAVX(TypeAVX512, OpAdd, false, false)
		VecAddSubScalarToAVX(TypeAVX512, OpAdd, true, false)
		VecAddSubScalarToAVX(TypeAVX512, OpSub, false, false)
		VecAddSubScalarToAVX(TypeAVX512, OpSub, true, false)
		VecAddSubScalarToAVX(TypeAVX512, OpSub, false, true)
		VecAddSubScalarToAVX(TypeAVX512, OpSub, true, true)

		VecNegToAVX(TypeAVX2, false)
		VecNegToAVX(TypeAVX2, true)

		VecNegToAVX(TypeAVX512, true)
		VecNegToAVX(TypeAVX512, false)

		VecMulWordToAVX(TypeAVX2, OpPure)
		VecMulWordToAVX(TypeAVX2, OpAdd)
		VecMulWordToAVX(TypeAVX2, OpSub)

		VecMulWordToAVX(TypeAVX512, OpPure)
		VecMulWordToAVX(TypeAVX512, OpAdd)
		VecMulWordToAVX(TypeAVX512, OpSub)

		VecMulScalarWordToAVX(TypeAVX2, OpPure)
		VecMulScalarWordToAVX(TypeAVX2, OpAdd)
		VecMulScalarWordToAVX(TypeAVX2, OpSub)

		VecMulScalarWordToAVX(TypeAVX512, OpPure)
		VecMulScalarWordToAVX(TypeAVX512, OpAdd)
		VecMulScalarWordToAVX(TypeAVX512, OpSub)

		VecMFormToAVX512()
		VecInvMFormToAVX512()

		VecMMulToAVX512(OpPure, false)
		VecMMulToAVX512(OpAdd, false)
		VecMMulToAVX512(OpSub, false)

		VecMMulToAVX512(OpPure, true)
		VecMMulToAVX512(OpAdd, true)
		VecMMulToAVX512(OpSub, true)

		VecMMulScalarToAVX512(OpPure, false)
		VecMMulScalarToAVX512(OpAdd, false)
		VecMMulScalarToAVX512(OpSub, false)

		VecMMulScalarToAVX512(OpPure, true)
		VecMMulScalarToAVX512(OpAdd, true)
		VecMMulScalarToAVX512(OpSub, true)

		VecSMulToAVX512(OpPure, false)
		VecSMulToAVX512(OpAdd, false)
		VecSMulToAVX512(OpSub, false)

		VecSMulToAVX512(OpPure, true)
		VecSMulToAVX512(OpAdd, true)
		VecSMulToAVX512(OpSub, true)

		VecSMulScalarToAVX512(OpPure, false)
		VecSMulScalarToAVX512(OpAdd, false)
		VecSMulScalarToAVX512(OpSub, false)

		VecSMulScalarToAVX512(OpPure, true)
		VecSMulScalarToAVX512(OpAdd, true)
		VecSMulScalarToAVX512(OpSub, true)

		VecReduceToAVX512()

		VecReduceFixedQToAVX(2, TypeAVX2)
		VecReduceFixedQToAVX(2, TypeAVX512)

		VecReduceFixedQToAVX(4, TypeAVX2)
		VecReduceFixedQToAVX(4, TypeAVX512)
	}

	if *ntt {
		// NTTConstants()

		// FwdNTTInPlacePow2StrideUnrollAVX512()
		// FwdNTTInPlacePow2UnrollAVX512()

		// InvNTTInPlacePow2StrideUnrollAVX512()
		// InvNTTInPlacePow2UnrollAVX512()
	}

	Generate()
}
