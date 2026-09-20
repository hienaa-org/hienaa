package num_test

import (
	crand "crypto/rand"
	"math"
	"testing"

	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

func TestIsPrime(t *testing.T) {
	composite := rSrc.Uint64N(1<<(num.MaxModulusBits/2)) * rSrc.Uint64N(1<<(num.MaxModulusBits/2))
	assert.False(t, num.IsPrime(composite))

	prime, err := crand.Prime(nil, num.MaxModulusBits-1)
	assert.NoError(t, err)
	assert.True(t, num.IsPrime(prime.Uint64()))
}

func TestFactor(t *testing.T) {
	x := rSrc.Uint64N(num.MaxModulus)
	primes, exps := num.Factor(x)

	xComp := uint64(1)
	for i := range primes {
		for j := uint64(0); j < exps[i]; j++ {
			xComp *= primes[i]
		}
	}

	assert.Equal(t, x, xComp)
}

func TestDivRound(t *testing.T) {
	assert.Equal(t, int64(0), num.DivRound[int64](1, math.MinInt64))
	assert.Equal(t, int64(1), num.DivRound[int64](2, 3))
	assert.Equal(t, int64(-1), num.DivRound[int64](-2, 3))
}
