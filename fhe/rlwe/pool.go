package rlwe

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/crt"
)

// ElementPool is a thin wrapper around [sync.Pool] with [rlwe.Element].
type ElementPool struct {
	pPool *sync.Pool
	sPool *sync.Pool
}

// NewElementPool creates a new [ElementPool].
func NewElementPool(params Parameters, hasAux, isNTT bool) *ElementPool {
	return &ElementPool{
		pPool: &sync.Pool{
			New: func() any {
				return NewPoly(params, hasAux, isNTT)
			},
		},
		sPool: &sync.Pool{
			New: func() any {
				return NewScalar(params, hasAux)
			},
		},
	}
}

// Get returns an element from the pool.
func (p *ElementPool) Get(eType crt.ElementType) *Element {
	switch eType {
	case crt.TypeScalar:
		return p.sPool.Get().(*Element)
	case crt.TypePoly:
		return p.pPool.Get().(*Element)
	}
	return nil
}

// Put returns an element to the pool.
func (p *ElementPool) Put(e *Element) {
	switch e.Type() {
	case crt.TypeScalar:
		p.sPool.Put(e)
	case crt.TypePoly:
		p.pPool.Put(e)
	}
}
