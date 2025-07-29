package rns

type Reducer struct {
	reducer reducer
}

type reducer interface {
	reduceTo(p, pOut *Poly)
}

// func newReducer(deg int, modpoly *Poly, modulus *mod.Modulus) *Reducer {

// }

// reducerNTTBuffer is a buffer for [cyclotomicReducerNTTModulus].
type reducerNTTBuffer struct {
	// pIn is a buffer for the input polynomial.
	pIn [][]uint64
	// pQuo is a buffer for the quotient polynomial.
	pQuo [][]uint64
	// pRem is a buffer for the remainder polynomial.
	pRem [][]uint64
}

func newReducerNTTBuffer(bufLen, sizeIn, sizeQuo, sizeRem int) reducerNTTBuffer {
	pIn := make([][]uint64, bufLen)
	pQuo := make([][]uint64, bufLen)
	pRem := make([][]uint64, bufLen)

	for i := range pIn {
		pIn[i] = make([]uint64, sizeIn)
		pQuo[i] = make([]uint64, sizeQuo)
		pRem[i] = make([]uint64, sizeRem)
	}

	return reducerNTTBuffer{
		pIn:  pIn,
		pQuo: pQuo,
		pRem: pRem,
	}
}
