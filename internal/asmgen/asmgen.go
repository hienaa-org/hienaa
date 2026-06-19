//go:generate go run . -vec -out ../../math/vec/asm_mod_amd64.s -stubs ../../math/vec/asm_mod_stub_amd64.go -pkg=vec
//go:generate go run . -ntt -out ../../math/dft/asm_ntt_pow2.s -stubs ../../math/dft/asm_ntt_pow2_stub_amd64.go -pkg=dft
package main

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

		VecAddSubScalarToAVX(TypeAVX2, OpAdd, false)
		VecAddSubScalarToAVX(TypeAVX2, OpAdd, true)
		VecAddSubScalarToAVX(TypeAVX2, OpSub, false)
		VecAddSubScalarToAVX(TypeAVX2, OpSub, true)

		VecAddSubScalarToAVX(TypeAVX512, OpAdd, false)
		VecAddSubScalarToAVX(TypeAVX512, OpAdd, true)
		VecAddSubScalarToAVX(TypeAVX512, OpSub, false)
		VecAddSubScalarToAVX(TypeAVX512, OpSub, true)

		VecNegToAVX(TypeAVX2, false)
		VecNegToAVX(TypeAVX2, true)

		VecNegToAVX(TypeAVX512, true)
		VecNegToAVX(TypeAVX512, false)

		VecMulScalarWordToAVX(TypeAVX2, OpPure)
		VecMulScalarWordToAVX(TypeAVX2, OpAdd)
		VecMulScalarWordToAVX(TypeAVX2, OpSub)

		VecMulScalarWordToAVX(TypeAVX512, OpPure)
		VecMulScalarWordToAVX(TypeAVX512, OpAdd)
		VecMulScalarWordToAVX(TypeAVX512, OpSub)

		VecMulScalarToAVX(TypeAVX2, OpPure)
		VecMulScalarToAVX(TypeAVX2, OpAdd)
		VecMulScalarToAVX(TypeAVX2, OpSub)

		VecMulScalarToAVX(TypeAVX512, OpPure)
		VecMulScalarToAVX(TypeAVX512, OpAdd)
		VecMulScalarToAVX(TypeAVX512, OpSub)

		VecSMulScalarToAVX(TypeAVX512, OpPure)
		VecSMulScalarToAVX(TypeAVX512, OpAdd)
		VecSMulScalarToAVX(TypeAVX512, OpSub)

		VecSMulScalarToAVX(TypeAVX512IFMA, OpPure)
		VecSMulScalarToAVX(TypeAVX512IFMA, OpAdd)
		VecSMulScalarToAVX(TypeAVX512IFMA, OpSub)

		VecMulWordToAVX(TypeAVX2, OpPure)
		VecMulWordToAVX(TypeAVX2, OpAdd)
		VecMulWordToAVX(TypeAVX2, OpSub)

		VecMulWordToAVX(TypeAVX512, OpPure)
		VecMulWordToAVX(TypeAVX512, OpAdd)
		VecMulWordToAVX(TypeAVX512, OpSub)

		VecMulToAVX(TypeAVX2, OpPure)
		VecMulToAVX(TypeAVX2, OpAdd)
		VecMulToAVX(TypeAVX2, OpSub)

		VecMulToAVX(TypeAVX512, OpPure)
		VecMulToAVX(TypeAVX512, OpAdd)
		VecMulToAVX(TypeAVX512, OpSub)

		VecSMulToAVX(TypeAVX512, OpPure)
		VecSMulToAVX(TypeAVX512, OpAdd)
		VecSMulToAVX(TypeAVX512, OpSub)

		VecSMulToAVX(TypeAVX512IFMA, OpPure)
		VecSMulToAVX(TypeAVX512IFMA, OpAdd)
		VecSMulToAVX(TypeAVX512IFMA, OpSub)

		VecReduceToAVX512()

		VecReduceFixedQToAVX(2, TypeAVX2)
		VecReduceFixedQToAVX(2, TypeAVX512)

		VecReduceFixedQToAVX(4, TypeAVX2)
		VecReduceFixedQToAVX(4, TypeAVX512)
	}

	if *ntt {
		NTTConstants()

		FwdNTTInPlacePow2UnrollAVX512()
		InvNTTInPlacePow2UnrollAVX512()
	}

	Generate()
}
