package bfv

// Packer packs a vector of uint64 into [*Plaintext].
type Packer struct {
}

// PackLen returns the length of packable vector.
func (p *Packer) PackLen() int {
	return 0
}

// SafeCopy returns a thread-safe copy.
func (p *Packer) SafeCopy() *Packer {
	return p
}

// Pack packs v.
// TODO: Panics when...
func (p *Packer) Pack(v []uint64) *Plaintext {
	panic("unimplemented")
}

// PackTo packs v to ptOut.
func (p *Packer) PackTo(ptOut *Plaintext, v []uint64) {
	panic("unimplemented")
}

// UnPack unpacks pt.
func (p *Packer) UnPack(pt *Plaintext) []uint64 {
	panic("unimplemented")
}

// UnPack unpacks pt to vOut.
func (p *Packer) UnPackTo(vOut []uint64, pt *Plaintext) {
	panic("unimplemented")
}
