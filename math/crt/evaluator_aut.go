package crt

// polyAutEvaluator is the evaluator for automorphism.
type polyAutEvaluator interface {
	// Aut returns aut_idx(p).
	// Panics when automorphism is invalid.
	// Notable cases include:
	//
	//	- In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
	//	- In any other rings, automorphism is not supported and it always panics.
	Aut(p *Poly, idx uint64) *Poly
	// AutTo computes pOut = aut_idx(p).
	// Panics when automorphism is invalid.
	// Notable cases include:
	//
	//	- In cyclotomic/autfixed rings, it panics when idx is not coprime with the cyclotomic order.
	//	- In any other rings, automorphism is not supported and it always panics.
	AutTo(pOut, p *Poly, idx uint64)
	// safeCopy returns a thread-safe copy.
	safeCopy() polyAutEvaluator
}

// polyAutEvaluatorBuffer is a buffer for [polyAutEvaluator].
type polyAutEvaluatorBuffer struct {
	p []uint64
}

// newPolyAutEvaluatorBuffer creates a new [polyAutEvaluatorBuffer].
func newPolyAutEvaluatorBuffer(rank int) polyAutEvaluatorBuffer {
	return polyAutEvaluatorBuffer{
		p: make([]uint64, rank),
	}
}

// polyAutEvaluatorPanic always panics.
// Used for [dft.Cyclic] rings.
type polyAutEvaluatorPanic struct{}

func (e *polyAutEvaluatorPanic) Aut(p *Poly, idx uint64) *Poly {
	panic("Aut: automorphism not supported in this ring")
}

func (e *polyAutEvaluatorPanic) AutTo(pOut, p *Poly, idx uint64) {
	panic("AutTo: automorphism not supported in this ring")
}

func (e *polyAutEvaluatorPanic) safeCopy() polyAutEvaluator {
	return &polyAutEvaluatorPanic{}
}
