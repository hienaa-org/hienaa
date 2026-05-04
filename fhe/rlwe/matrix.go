package rlwe

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/vec"
)

// BSGSParams is the parameters for the BSGS matrix multiplication algorithm.
type BSGSParams struct {
	BabyStep  []int
	GiantStep []int
}

// PlainMatrix is a plaintext matrix.
type PlainMatrix struct {
	// map for the diagonal encodings.
	diag map[int]map[int]*Element
	// parameters for the BSGS matrix multiplication algorithm.
	bsgsParams BSGSParams
}

// NewPlainMatrix creates a new [PlainMatrix]
func NewPlainMatrix(diag map[int]map[int]*Element, bsgsParams BSGSParams) *PlainMatrix {
	return &PlainMatrix{
		diag:       diag,
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
	if !(slices.Contains(mat.bsgsParams.BabyStep, bs) && slices.Contains(mat.bsgsParams.GiantStep, gs)) {
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
			if gs != 1 && !slices.Contains(res, gs) {
				res = append(res, gs)
			}

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

	// Decompose the mask.
	ctLen := ct.BaseModLen()
	dcmp := op.dcmpPool.Get()
	defer op.dcmpPool.Put(dcmp)
	dcmp = dcmp.Slice(vec.Range(0, op.dcmp.DecomposeLen(ctLen))...)
	dcmp = dcmp.WithModLen(ctLen, op.dcmp.AuxModLen(ctLen))

	bufNTT := op.pPool.Get()
	defer op.pPool.Put(bufNTT)
	bufNTT = bufNTT.WithModLen(ctLen, 0)
	if !ct.Mask.IsNTT() {
		bufNTT.CopyFrom(ct.Mask)
	} else {
		op.plainOp.InvNTTTo(bufNTT, ct.Mask)
	}
	op.dcmp.DecomposeTo(dcmp, bufNTT, true)

	// Compute the baby steps.
	babyStep := make(map[int]*Ciphertext, len(mat.bsgsParams.BabyStep))
	for _, bs := range mat.bsgsParams.BabyStep {
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

			if bs == 1 {
				babyStep[bs].CopyFrom(ct)
			} else {
				op.HoistedAutTo(babyStep[bs], dcmp, ct, atk[bs], true)
			}
		}
	}

	// Compute the giant steps.
	bsAcc := op.ctPool.Get()
	defer op.ctPool.Put(bsAcc)
	bsAcc = bsAcc.WithModLen(ctLen, 0)
	for gs, bsMap := range mat.diag {
		bsAcc.Clear()

		for bs, diag := range bsMap {
			diag = diag.WithModLen(ctLen, 0)
			op.MulAddElementTo(bsAcc, babyStep[bs], diag)
		}

		if gs == 1 {
			cOut.CopyFrom(bsAcc)
		} else {
			op.AutTo(bsAcc, bsAcc, atk[gs], true)
			op.AddTo(cOut, cOut, bsAcc)
		}
	}
}
