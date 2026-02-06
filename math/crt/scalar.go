package crt

// // Scalar represents a scalar in CRT basis.
// type Scalar []uint64

// // NewScalar creates a new [Scalar].
// func NewScalar[T num.Integer | *big.Int](x T, mod []*num.Modulus) Scalar {
// 	r := make(Scalar, len(mod))

// 	var z T
// 	switch any(z).(type) {
// 	case *big.Int:
// 		u := any(x).(*big.Int)
// 		q, t := new(big.Int), new(big.Int)
// 		for i := range r {
// 			q.SetUint64(mod[i].Value())
// 			r[i] = t.Mod(u, q).Uint64()
// 		}

// 	case int8:
// 		u := any(x).(int8)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case int16:
// 		u := any(x).(int16)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case int32:
// 		u := any(x).(int32)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case int64:
// 		u := any(x).(int64)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case int:
// 		u := any(x).(int)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}

// 	case uint8:
// 		u := any(x).(uint8)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case uint16:
// 		u := any(x).(uint16)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case uint32:
// 		u := any(x).(uint32)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case uint64:
// 		u := any(x).(uint64)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	case uint:
// 		u := any(x).(uint)
// 		for i := range r {
// 			r[i] = num.Reduce(u, mod[i])
// 		}
// 	}

// 	return r
// }

// // isTernaryScalarToOperable panics if inputs are not consistent.
// func mustScalarToOperable(modLen int, x ...Scalar) {
// 	for i := range x {
// 		if len(x[i]) != modLen {
// 			panic("input(s) not consistent")
// 		}
// 	}
// }

// // AddScalar returns x0 + x1.
// func AddScalar(x0, x1 Scalar, mod []*num.Modulus) Scalar {
// 	xOut := make(Scalar, len(mod))
// 	AddScalarTo(xOut, x0, x1, mod)
// 	return xOut
// }

// // AddScalarTo computes xOut = x0 + x1.
// func AddScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
// 	mustScalarToOperable(len(mod), xOut, x0, x1)

// 	for i := range mod {
// 		xOut[i] = num.Add(x0[i], x1[i], mod[i])
// 	}
// }

// // SubScalar returns x0 - x1.
// func SubScalar(x0, x1 Scalar, mod []*num.Modulus) Scalar {
// 	xOut := make(Scalar, len(mod))
// 	SubScalarTo(xOut, x0, x1, mod)
// 	return xOut
// }

// // SubScalarTo computes xOut = x0 - x1.
// func SubScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
// 	mustScalarToOperable(len(mod), xOut, x0, x1)

// 	for i := range mod {
// 		xOut[i] = num.Sub(x0[i], x1[i], mod[i])
// 	}
// }

// // NegScalar returns -x.
// func NegScalar(x Scalar, mod []*num.Modulus) Scalar {
// 	xOut := make(Scalar, len(mod))
// 	NegScalarTo(xOut, x, mod)
// 	return xOut
// }

// // NegScalarTo computes xOut = -x.
// func NegScalarTo(xOut, x Scalar, mod []*num.Modulus) {
// 	mustScalarToOperable(len(mod), xOut, x)

// 	for i := range mod {
// 		xOut[i] = num.Neg(x[i], mod[i])
// 	}
// }

// // MulScalar returns x0 * x1.
// func MulScalar(x0, x1 Scalar, mod []*num.Modulus) Scalar {
// 	xOut := make(Scalar, len(mod))
// 	MulScalarTo(xOut, x0, x1, mod)
// 	return xOut
// }

// // MulScalarTo computes xOut = x0 * x1.
// func MulScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
// 	mustScalarToOperable(len(mod), xOut, x0, x1)

// 	for i := range mod {
// 		xOut[i] = num.Mul(x0[i], x1[i], mod[i])
// 	}
// }

// // MulAddScalarTo computes xOut += x0 * x1.
// func MulAddScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
// 	mustScalarToOperable(len(mod), xOut, x0, x1)

// 	for i := range mod {
// 		xOut[i] = num.Add(xOut[i], num.Mul(x0[i], x1[i], mod[i]), mod[i])
// 	}
// }

// // MulSubScalarTo computes xOut -= x0 * x1.
// func MulSubScalarTo(xOut, x0, x1 Scalar, mod []*num.Modulus) {
// 	mustScalarToOperable(len(mod), xOut, x0, x1)

// 	for i := range mod {
// 		xOut[i] = num.Sub(xOut[i], num.Mul(x0[i], x1[i], mod[i]), mod[i])
// 	}
// }

// // AsBigScalar returns s as *[big.Int].
// func AsBigScalar(s Scalar, mod []*num.Modulus) *big.Int {
// 	if len(mod) != len(s) {
// 		panic("input(s) not consistent")
// 	}

// 	modBig := make([]*big.Int, len(mod))
// 	modProd := big.NewInt(1)
// 	for i := range mod {
// 		modBig[i] = new(big.Int).SetUint64(mod[i].Value())
// 		modProd.Mul(modProd, modBig[i])
// 	}
// 	modProdHalf := new(big.Int).Rsh(modProd, 1)

// 	gadget := make([]*big.Int, len(mod))
// 	for i := range modBig {
// 		qStar := new(big.Int).Div(modProd, modBig[i])
// 		qStarInv := new(big.Int).ModInverse(qStar, modBig[i])
// 		gadget[i] = new(big.Int).Mul(qStar, qStarInv)
// 		gadget[i].Mod(gadget[i], modProd)
// 	}

// 	sBig := big.NewInt(0)
// 	for i := range mod {
// 		c := new(big.Int).SetUint64(s[i])
// 		c.Mul(c, gadget[i])
// 		sBig.Add(sBig, c)
// 	}
// 	sBig.Mod(sBig, modProd)
// 	if sBig.Cmp(modProdHalf) > 0 {
// 		sBig.Sub(sBig, modProd)
// 	}

// 	return sBig
// }
