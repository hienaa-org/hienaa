package dft

// Form represents the form of a polynomial.
type Form bool

const (
	// FormCoeff represents the standard, coefficient form.
	FormCoeff Form = false
	// FormNTT represents the NTT form.
	FormNTT Form = true
)
