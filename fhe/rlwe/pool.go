package rlwe

import (
	"sync"

	"github.com/hienaa-org/hienaa/math/crt"
)

// ElementPool is a thin wrapper around [sync.Pool] with [rlwe.Element].
type ElementPool struct {
	pool *sync.Pool
}

// NewElementPool creates a new [ElementPool].
func NewElementPool(params Parameters, hasAux, isNTT bool, eType crt.ElementType) *ElementPool {
	var pool *sync.Pool

	switch eType {
	case crt.TypeScalar:
		pool = &sync.Pool{
			New: func() any {
				return NewScalar(params, hasAux)
			},
		}
	case crt.TypePoly:
		pool = &sync.Pool{
			New: func() any {
				return NewPoly(params, hasAux, isNTT)
			},
		}
	}

	return &ElementPool{pool: pool}
}

// Get returns an element from the pool.
func (p *ElementPool) Get() *Element {
	return p.pool.Get().(*Element)
}

// Put returns an element to the pool.
func (p *ElementPool) Put(e *Element) {
	p.pool.Put(e)
}
