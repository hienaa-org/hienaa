package matrix

import (
	"slices"

	"github.com/hienaa-org/hienaa/fhe/rlwe"
)

// TODO: Implement the matrix multiplication algorithm.

// BSGSParams is the parameters for the BSGS matrix multiplication algorithm.
type BSGSParams struct {
	babyStep  []int
	giantStep []int
}

// PlainMatrix is a plaintext matrix.
type PlainMatrix struct {
	diag       map[int]*rlwe.Element
	bsgsParams BSGSParams
}

// NewPlainMatrix creates a new plaintext matrix.
func NewPlainMatrix(params rlwe.Parameters, bsgsParams BSGSParams) *PlainMatrix {
	return &PlainMatrix{
		diag:       make(map[int]*rlwe.Element),
		bsgsParams: bsgsParams,
	}
}

// RequiredRotationIndex returns the required rotation indices for the BSGS matrix multiplication algorithm.
func RequiredRotationIndex(mat *PlainMatrix) [][]int {
	dimLen := len(mat.bsgsParams.babyStep)

	res := make([][]int, 0)
	for dim, gs := range mat.bsgsParams.giantStep {
		for i := 0; i < gs; i++ {
			needGs := false
			for j := 0; j < mat.bsgsParams.babyStep[dim]; j++ {
				if mat.diag[i*gs+j] == nil {
					needGs = true

					rotIdx := make([]int, dimLen)
					rotIdx[dim] = j

					if j != 0 && !slices.ContainsFunc(res, func(x []int) bool {
						return slices.Equal(x, rotIdx)
					}) {
						res = append(res, rotIdx)
					}
				}

				if i != 0 && !needGs {
					rotIdx := make([]int, dimLen)
					rotIdx[dim] = i * gs
					res = append(res, rotIdx)
				}
			}
		}
	}

	slices.SortFunc(res, func(a, b []int) int {
		return slices.Compare(a, b)
	})

	return res
}
