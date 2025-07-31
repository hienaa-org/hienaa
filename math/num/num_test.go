package num_test

import (
	crand "crypto/rand"
	"testing"

	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/stretchr/testify/assert"
)

var (
	rSrc = csprng.NewUniformSamplerWithSeed(nil)
)

func TestIsPrime(t *testing.T) {
	composite := rSrc.SampleN(1<<(num.MaxModulusBits/2)) * rSrc.SampleN(1<<(num.MaxModulusBits/2))
	assert.False(t, num.IsPrime(composite))

	prime, err := crand.Prime(rSrc, num.MaxModulusBits-1)
	assert.NoError(t, err)
	assert.True(t, num.IsPrime(prime.Uint64()))
}

func TestFactor(t *testing.T) {
	x := rSrc.SampleN(num.MaxModulus)
	factors := num.Factor(x)

	xComp := uint64(1)
	for f, e := range factors {
		for i := uint64(0); i < e; i++ {
			xComp *= f
		}
	}

	assert.Equal(t, x, xComp)
}
