package num_test

import (
	crand "crypto/rand"
	"math/rand"
	"testing"

	"github.com/hienaa-org/hienaa/math/mod"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = rand.New(rand.NewSource(0))
)

func TestIsPrime(t *testing.T) {
	composite := uint64(rSrc.Uint32()) * uint64(rSrc.Uint32())
	assert.False(t, num.IsPrime(composite))

	prime, err := crand.Prime(rSrc, num.Log2(mod.MaxModulus)-1)
	assert.NoError(t, err)
	assert.True(t, num.IsPrime(prime.Uint64()))
}

func TestFactor(t *testing.T) {
	x := rSrc.Uint64()
	factors := num.Factor(x)

	xComp := uint64(1)
	for f, e := range factors {
		for i := uint64(0); i < e; i++ {
			xComp *= f
		}
	}

	assert.Equal(t, x, xComp)
}
