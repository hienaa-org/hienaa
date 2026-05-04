package rlwe

import (
	"slices"

	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/num"
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

// Diag returns the diagonal element for the given giant step and baby step.
func (mat *PlainMatrix) Diag() map[int]map[int]*Element {
	return mat.diag
}

// BSGSParams returns the BSGS parameters.
func (mat *PlainMatrix) BSGSParams() BSGSParams {
	return mat.bsgsParams
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
func (op *Operator) MulPlainMatrix(mat *PlainMatrix, ct *Ciphertext, atk map[int]*AutomorphismKey, isNTT bool) *Ciphertext {
	cOut := NewCiphertextCustom(ct.Rank(), ct.BaseModLen(), ct.AuxModLen(), true)
	op.MulPlainMatrixTo(cOut, mat, ct, atk, isNTT)
	return cOut
}

// TODO: Implement double-hoisting BSGS matrix multiplication algorithm.
// MulPlainMatrixTo computes ctOut = mat * ct using Halevi-Shoup BSGS matrix multiplication algorithm.
func (op *Operator) MulPlainMatrixTo(cOut *Ciphertext, mat *PlainMatrix, ct *Ciphertext, atk map[int]*AutomorphismKey, isNTT bool) {
	if ct.AuxModLen() > 0 {
		panic("input ciphertext must not have auxiliary modulus")
	}

	// Define parameters.
	ctLen := ct.BaseModLen()
	auxLen := op.dcmp.AuxModLen(ctLen)

	// Compute constants.
	auxMod := op.pPool.Get(crt.TypeScalar)
	defer op.pPool.Put(auxMod)
	auxMod = auxMod.WithModLen(ctLen, 0)
	for i := 0; i < ctLen; i++ {
		auxMod.Value.Coeffs[i][0] = 1
		for j := 0; j < auxLen; j++ {
			auxMod.Value.Coeffs[i][0] = num.Mul(auxMod.Value.Coeffs[i][0], op.Params.auxMod[j].Value(), op.Params.baseMod[i])
		}
	}

	// Decompose the mask.
	dcmp := op.dcmpPool.Get()
	defer op.dcmpPool.Put(dcmp)
	dcmp = dcmp.Slice(vec.Range(0, op.dcmp.DecomposeLen(ctLen))...)
	dcmp = dcmp.WithModLen(ctLen, auxLen)

	bufNTT := op.pPool.Get(crt.TypePoly)
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
			babyStep[bs] = babyStep[bs].WithModLen(ctLen, auxLen)

			if bs == 1 {
				bsBase := babyStep[bs].WithModLen(ctLen, 0)
				bsAux := babyStep[bs].WithModLen(0, auxLen)
				if !ct.IsNTT() {
					op.FwdNTTTo(bsBase, ct)
				} else {
					bsBase.CopyFrom(ct)
				}

				op.MulElementTo(bsBase, bsBase, auxMod)
				bsAux.Clear()
			} else {
				op.HoistedAutLazyTo(babyStep[bs], dcmp, ct, atk[bs], true)
			}
		}
	}

	// Compute the giant steps.
	bsAcc := op.ctPool.Get()
	defer op.ctPool.Put(bsAcc)
	bsAcc = bsAcc.WithModLen(ctLen, auxLen)
	bsBase := bsAcc.WithModLen(ctLen, 0)

	gsAcc := op.ctPool.Get()
	defer op.ctPool.Put(gsAcc)
	gsAcc = gsAcc.WithModLen(ctLen, auxLen)

	cOut.Clear()
	cOut.Body.Value.IsNTT = true
	cOut.Mask.Value.IsNTT = true
	for gs, bsMap := range mat.diag {
		bsAcc.Clear()

		for bs, diag := range bsMap {
			diag = diag.WithModLen(ctLen, auxLen)
			op.MulAddElementTo(bsAcc, babyStep[bs], diag)
		}

		if gs != 1 {
			if op.Params.HasAuxModulus() {
				op.DivByAuxModulusTo(bsBase, bsAcc, true)
			}
			op.AutLazyTo(bsAcc, bsBase, atk[gs], true)
		}
		op.AddTo(gsAcc, gsAcc, bsAcc)
	}

	if op.Params.HasAuxModulus() {
		op.DivByAuxModulusTo(cOut, gsAcc, isNTT)
	} else if !isNTT {
		op.InvNTTTo(cOut, gsAcc)
	} else {
		cOut.CopyFrom(gsAcc)
	}
}
