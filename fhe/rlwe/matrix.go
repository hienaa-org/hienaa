package rlwe

import (
	"slices"
)

// TODO: Is this the right place for matrix?

// BSGSParams is the parameters for the BSGS matrix multiplication algorithm.
type BSGSParams struct {
	babyStep  []int
	giantStep []int
}

// PlainMatrix is a plaintext matrix.
type PlainMatrix struct {
	// map for the diagonal encodings.
	diag map[int]map[int]*Element
	// parameters for the BSGS matrix multiplication algorithm.
	bsgsParams BSGSParams
}

// NewPlainMatrix creates a new [PlainMatrix]
func NewPlainMatrix(params Parameters, bsgsParams BSGSParams) *PlainMatrix {
	return &PlainMatrix{
		diag:       make(map[int]map[int]*Element),
		bsgsParams: bsgsParams,
	}
}

// GetDiag returns the diagonal element for the given giant step and baby step.
func (mat *PlainMatrix) GetDiag(gs, bs int) *Element {
	return mat.diag[gs][bs]
}

// SetDiag sets the diagonal element for the given giant step and baby step.
func (mat *PlainMatrix) SetDiag(gs, bs int, val *Element) {
	// Sanity check.
	if !(slices.Contains(mat.bsgsParams.babyStep, bs) && slices.Contains(mat.bsgsParams.giantStep, gs)) {
		panic("invalid giant step or baby step")
	}

	if mat.diag[gs] == nil {
		mat.diag[gs] = make(map[int]*Element)
	}
	mat.diag[gs][bs] = val
}

// RequiredAutIndex returns the required automorphism indices for the BSGS matrix multiplication algorithm.
func RequiredAutIndex(mat *PlainMatrix) []int {
	res := make([]int, 0)

	for gs, bsMap := range mat.diag {
		if bsMap != nil {
			res = append(res, gs)

			for bs, val := range bsMap {
				if bs != 1 && val != nil && !slices.Contains(res, bs) {
					res = append(res, bs)
				}
			}
		}
	}
	slices.Sort(res)

	return res
}

// MulPlainMatrix computes ctOut = mat * ct.
func (op *Operator) MulPlainMatrix(mat *PlainMatrix, ct *Ciphertext, atk map[int]*AutomorphismKey) *Ciphertext {
	cOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), true)
	op.MulPlainMatrixTo(cOut, mat, ct, atk)
	return cOut
}

// TODO: Implement double-hoisting BSGS matrix multiplication algorithm.
// MulPlainMatrixTo computes ctOut = mat * ct using Halevi-Shoup BSGS matrix multiplication algorithm.
func (op *Operator) MulPlainMatrixTo(cOut *Ciphertext, mat *PlainMatrix, ct *Ciphertext, atk map[int]*AutomorphismKey) {
	if ct.AuxModLen() > 0 {
		panic("input ciphertext must not have auxiliary modulus")
	}

	ctLen := ct.BaseModLen()

	babyStep := make(map[int]*Ciphertext, len(mat.bsgsParams.babyStep))
	for _, bs := range mat.bsgsParams.babyStep {
		check := false
		for _, bsMap := range mat.diag {
			if bsMap[bs] != nil {
				check = true
				break
			}
		}

		if check {
			babyStep[bs] = op.ctPool.Get()
			defer op.ctPool.Put(babyStep[bs])
			babyStep[bs] = babyStep[bs].WithModLen(ctLen, 0)

			op.AutTo(babyStep[bs], ct, atk[bs], true)
		}
	}

	bsAcc := op.ctPool.Get()
	defer op.ctPool.Put(bsAcc)
	bsAcc = bsAcc.WithModLen(ctLen, 0)

	for gs, bsMap := range mat.diag {
		bsAcc.Clear()
		check := false

		for bs, diag := range bsMap {
			if diag != nil {
				check = true
				diag = diag.WithModLen(ctLen, 0)
				op.MulAddElementTo(bsAcc, babyStep[bs], diag)
			}
		}

		if check {
			if gs == 1 {
				cOut.CopyFrom(bsAcc)
			} else {
				op.AutTo(bsAcc, bsAcc, atk[gs], true)
				op.AddTo(cOut, cOut, bsAcc)
			}
		}
	}
}
