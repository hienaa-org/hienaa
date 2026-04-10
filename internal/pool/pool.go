package pool

import "sync"

// Pool is a thin wrapper around [sync.Pool] using Generics.
type Pool[T any] struct {
	pool sync.Pool
}

// NewPool creates a new [Pool].
func NewPool[T any](New func() T) *Pool[T] {
	return &Pool[T]{
		sync.Pool{
			New: func() any {
				return any(New())
			},
		},
	}
}

// Put adds x to the pool.
func (p *Pool[T]) Put(x T) {
	p.pool.Put(x)
}

// Get selects an arbitrary item from the Pool, removes it from the Pool, and returns it to the caller.
// Get may choose to ignore the pool and treat it as empty.
// Callers should not assume any relation between values passed to Pool.Put and the values returned by Get.
func (p *Pool[T]) Get() T {
	return p.pool.Get().(T)
}
