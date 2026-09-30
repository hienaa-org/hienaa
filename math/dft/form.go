package dft

// Form represents the form of a polynomial.
type Form byte

const (
	// FormCoeff represents the standard, coefficient form.
	FormCoeff Form = iota
	// FormNTT represents the NTT form.
	FormNTT
)
