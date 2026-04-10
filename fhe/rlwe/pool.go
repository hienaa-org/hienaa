package rlwe

import (
	"github.com/hienaa-org/hienaa/internal/pool"
	"github.com/hienaa-org/hienaa/math/crt"
)

// ElementPool is a thin wrapper around [pool.Pool] with [rlwe.Element].
type ElementPool struct {
	pPool *pool.Pool[*Element]
	sPool *pool.Pool[*Element]
}

// NewElementPool creates a new [ElementPool].
func NewElementPool(params Parameters, hasAux, isNTT bool) *ElementPool {
	return &ElementPool{
		pPool: pool.NewPool(func() *Element {
			return NewPoly(params, hasAux, isNTT)
		}),
		sPool: pool.NewPool(func() *Element {
			return NewScalar(params, hasAux)
		}),
	}
}

// Get returns an element from the pool.
func (p *ElementPool) Get(eType crt.ElementType) *Element {
	switch eType {
	case crt.TypeScalar:
		return p.sPool.Get()
	case crt.TypePoly:
		return p.pPool.Get()
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
