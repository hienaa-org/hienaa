package csprng

import (
	"crypto/aes"
	"crypto/cipher"
	"crypto/rand"
	"crypto/sha512"
	"math"
	"unsafe"

	"github.com/hienaa-org/hienaa/math/num"
)

// bufSize is the default buffer size of [UniformSampler].
const bufSize = 8192

// UniformSampler samples values from uniform distribution.
// This uses AES-CTR as a underlying prng.
type UniformSampler struct {
	prng cipher.Stream

	buf [bufSize]byte
	ptr uint64
}

// NewUniformSampler creates a new [UniformSampler].
//
// Panics when read from crypto/rand or AES initialization fails.
func NewUniformSampler() *UniformSampler {
	var seed [32]byte
	if _, err := rand.Read(seed[:]); err != nil {
		panic(err)
	}

	return NewUniformSamplerWithSeed(seed[:])
}

// NewUniformSamplerWithSeed creates a new [UniformSampler], with user supplied seed.
//
// Panics when AES initialization fails.
func NewUniformSamplerWithSeed(seed []byte) *UniformSampler {
	r := sha512.Sum384(seed)

	block, err := aes.NewCipher(r[:32])
	if err != nil {
		panic(err)
	}

	prng := cipher.NewCTR(block, r[32:])

	return &UniformSampler{
		prng: prng,

		buf: [bufSize]byte{},
		ptr: bufSize,
	}
}

// Read implements the [io.Reader] interface.
// It always succeeds, with n == len(p) and err == nil.
func (s *UniformSampler) Read(p []byte) (n int, err error) {
	clear(p)
	s.prng.XORKeyStream(p, p)
	return len(p), nil
}

// Sample uniformly samples a random value.
func (s *UniformSampler) Sample[T num.Integer]() T {
	bytes := uint64(unsafe.Sizeof(T(0)))

	if s.ptr+bytes > bufSize {
		s.prng.XORKeyStream(s.buf[:], s.buf[:])
		s.ptr = 0
	}

	var r uint64
	switch bytes {
	case 1:
		r |= uint64(s.buf[s.ptr+0]) << 0
	case 2:
		r |= uint64(s.buf[s.ptr+0]) << 0
		r |= uint64(s.buf[s.ptr+1]) << 8
	case 4:
		r |= uint64(s.buf[s.ptr+0]) << 0
		r |= uint64(s.buf[s.ptr+1]) << 8
		r |= uint64(s.buf[s.ptr+2]) << 16
		r |= uint64(s.buf[s.ptr+3]) << 24
	case 8:
		r |= uint64(s.buf[s.ptr+0]) << 0
		r |= uint64(s.buf[s.ptr+1]) << 8
		r |= uint64(s.buf[s.ptr+2]) << 16
		r |= uint64(s.buf[s.ptr+3]) << 24
		r |= uint64(s.buf[s.ptr+4]) << 32
		r |= uint64(s.buf[s.ptr+5]) << 40
		r |= uint64(s.buf[s.ptr+6]) << 48
		r |= uint64(s.buf[s.ptr+7]) << 56
	}
	s.ptr += bytes

	return T(r)
}

// SampleN uniformly samples a random value in [start, end).
//
// Panics if start >= end.
func (s *UniformSampler) SampleN[T num.Integer](start, end T) T {
	if start >= end {
		panic("range invalid")
	}

	n := uint64(end) - uint64(start)
	bound := -n % n
	for {
		r := s.Sample[uint64]()
		if r >= bound {
			return T(uint64(start) + r%n)
		}
	}
}

// SampleFloat samples a random float64 value in [0, 1).
func (s *UniformSampler) SampleFloat() float64 {
	r := s.Sample[uint64]() % (1 << floatPrec)
	rf := math.Float64frombits(r | ((1023 + floatPrec) << floatPrec))
	return (rf / (1 << floatPrec)) - 1
}
