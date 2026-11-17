package regev

import (
	"math/big"

	"github.com/hienaa-org/hienaa/fhe/internal/gnum"
	"github.com/hienaa-org/hienaa/fhe/internal/gr"
	"github.com/hienaa-org/hienaa/math/crt"
	"github.com/hienaa-org/hienaa/math/csprng"
	"github.com/hienaa-org/hienaa/math/dft"
	"github.com/hienaa-org/hienaa/math/num"
	"github.com/hienaa-org/hienaa/math/vec"
)

// gaussIntNthRoot returns a primitive n-th root of unity of the Gaussian integers modulo q.
// The modulus q must be a power of a prime congruent to 3 modulo 4.
//
// Panics when n <= 0, q is not of this form, or no such root exists.
func gaussIntNthRoot(n int, q *num.Modulus) gnum.GaussianInt {
	if n <= 0 {
		panic("n must be positive")
	}

	qPrimes, _ := num.Factor(q.Value())
	if !(len(qPrimes) == 1 && qPrimes[0]%4 == 3) {
		panic("q must be power of prime that is 3 modulo 4")
	}

	if n == 1 {
		return gnum.GaussianInt{Real: 1}
	}

	ord := new(big.Int).SetUint64(qPrimes[0])
	ord.Mul(ord, ord)
	ord.Sub(ord, big.NewInt(1))
	ord.Mul(ord, new(big.Int).SetUint64(q.Value()/qPrimes[0]))

	rootOrd, rootOrdMod := new(big.Int).DivMod(ord, big.NewInt(int64(n)), new(big.Int))
	if rootOrdMod.Sign() != 0 {
		panic("no nth root of unity exists")
	}

	primes, _ := num.Factor(n)
	checkOrd := make([]*big.Int, len(primes))
	for i := range checkOrd {
		checkOrd[i] = new(big.Int).Div(ord, big.NewInt(int64(primes[i])))
	}

	rSrc := csprng.NewUniformSamplerWithSeed(nil)
	for {
		genEl := gnum.GaussianInt{
			Real: rSrc.SampleN(0, q.Value()),
			Imag: rSrc.SampleN(0, q.Value()),
		}

		if genEl.Real%qPrimes[0] == 0 && genEl.Imag%qPrimes[0] == 0 {
			continue
		}

		ok := true
		for i := range checkOrd {
			genElPow := gnum.Exp(genEl, checkOrd[i], q)
			if genElPow.Real == 1 && genElPow.Imag == 0 {
				ok = false
				break
			}
		}

		if ok {
			root := gnum.Exp(genEl, rootOrd, q)
			if !(root.Real == 1 && root.Imag == 0) {
				return root
			}
		}
	}
}

// galoisRingNthRoot returns a primitive n-th root of unity of the Galois ring.
// n must divide p^d - 1, where p is the underlying prime and d is the rank of r.
//
// Panics when n <= 0 or n does not divide p^d - 1.
func galoisRingNthRoot(r *gr.GaloisRing, n int) *gr.Element {
	if n <= 0 {
		panic("n must be positive")
	}

	if n == 1 {
		return r.NewElementFrom(1)
	}

	qPrimes, _ := num.Factor(r.Mod())
	ord := new(big.Int).SetUint64(qPrimes[0])
	ord.Exp(ord, big.NewInt(int64(r.Rank())), nil)
	ord.Sub(ord, big.NewInt(1))

	rootOrd, rootOrdMod := new(big.Int).DivMod(ord, big.NewInt(int64(n)), new(big.Int))
	if rootOrdMod.Sign() != 0 {
		panic("n must divide p^d - 1")
	}

	pExp := new(big.Int).SetUint64(r.Mod() / qPrimes[0])
	ord.Mul(ord, pExp)
	rootOrd.Mul(rootOrd, pExp)

	primes, _ := num.Factor(n)
	checkOrd := make([]*big.Int, len(primes))
	for i := range checkOrd {
		checkOrd[i] = new(big.Int).Div(ord, big.NewInt(int64(primes[i])))
	}

	genEl := r.NewElement()
	root := r.NewElement()
	one := r.NewElementFrom(1)
	genElPow := r.NewElement()
	genCoeffs := genEl.Coeffs()

	rSrc := csprng.NewUniformSamplerWithSeed(nil)
	for {
		for i := 0; i < r.Rank(); i++ {
			genCoeffs[i] = rSrc.SampleN(0, r.Mod())
		}

		if !r.IsUnit(genEl) {
			continue
		}

		ok := true
		for i := range checkOrd {
			r.ExpTo(genElPow, genEl, checkOrd[i])
			if genElPow.IsEqual(one) {
				ok = false
				break
			}
		}

		if ok {
			r.ExpTo(root, genEl, rootOrd)
			if !root.IsEqual(one) {
				return root
			}
		}
	}
}

// resolutionOfUnity returns the resolution of unity of a prime cyclotomic ring, modulo prime^exp.
// Panics if cyclotomic index and prime is not prime.
func resolutionOfUnity(params dft.RingParameters, prime, exp uint64) []uint64 {
	if !num.IsPrime(params.CycloIdx()) {
		panic("cycloIdx must be prime")
	} else if !num.IsPrime(prime) {
		panic("modulus must be prime")
	}

	ord := num.Order(prime, uint64(params.CycloIdx()))
	rank := int(num.Totient(uint64(params.CycloIdx())) / ord)
	r := gr.NewGaloisRing(num.Exp(prime, exp, nil), int(ord))
	root := galoisRingNthRoot(r, params.CycloIdx())

	// Compute the first factor of the cyclotomic polynomial.
	factorPoly := make([]*gr.Element, int(ord)+1)
	tmpPoly := make([]*gr.Element, int(ord)+1)
	for i := 0; i <= int(ord); i++ {
		if i == 0 {
			factorPoly[i] = r.NewElementFrom(1)
		} else {
			factorPoly[i] = r.NewElement()
		}
		tmpPoly[i] = r.NewElement()
	}

	rootPow := r.NewElement()
	for i := 0; i < int(ord); i++ {
		r.ExpTo(rootPow, root, num.Exp(prime, uint64(i), params.CycloIdxMod()))

		// tmpPoly = factor * rootPow
		for j := 0; j <= i; j++ {
			r.MulTo(tmpPoly[j], factorPoly[j], rootPow)
		}

		// factor = factor * X
		for j := i + 1; j > 0; j-- {
			factorPoly[j].CopyFrom(factorPoly[j-1])
		}
		factorPoly[0].Clear()

		// factor = factor - tmpPoly
		for j := 0; j <= i; j++ {
			r.SubTo(factorPoly[j], factorPoly[j], tmpPoly[j])
		}
	}

	factor := make([]uint64, int(ord)+1)
	factorSigned := make([]int64, int(ord)+1)
	for i := 0; i <= int(ord); i++ {
		factor[i] = uint64(factorPoly[i].Coeffs()[0])
		factorSigned[i] = int64(factorPoly[i].Coeffs()[0])
	}

	// remPoly = (cycloPoly/factor) % factor.
	modulus := num.NewModulus(num.Exp(prime, exp, nil))
	factorReducer := crt.NewLongDivReducer(len(params.ModPoly()), []*num.Modulus{modulus}, factorSigned)
	cycloPolyMod := crt.NewElement(len(params.ModPoly()), 1, crt.TypePoly, dft.FormCoeff)
	vec.ReduceTo(cycloPolyMod.Coeffs[0], params.ModPoly(), modulus)
	quoPoly := factorReducer.Quotient(cycloPolyMod)
	remPoly := factorReducer.Reduce(quoPoly)

	// remEl = (1/remPoly) % factor.
	rFactor := gr.NewGaloisRingCustom(modulus.Value(), factorSigned)
	remEl := rFactor.NewElement()
	copy(remEl.Coeffs(), remPoly.Coeffs[0])
	rFactor.InvTo(remEl, remEl)

	// Compute the resolution of unity.
	ntt := crt.NewOperator(dft.NewCyclicParameters(num.NextProdPower(params.CycloIdx()+1, []int{2})), []*num.Modulus{modulus})
	resolPoly := ntt.NewPoly(dft.FormCoeff)
	remInvPoly := ntt.NewPoly(dft.FormCoeff)

	// resolPoly = quoPoly * remEl.
	copy(resolPoly.Coeffs[0][:len(quoPoly.Coeffs[0])], quoPoly.Coeffs[0])
	copy(remInvPoly.Coeffs[0][:len(remEl.Coeffs())], remEl.Coeffs())

	ntt.FwdNTTTo(resolPoly, resolPoly)
	ntt.FwdNTTTo(remInvPoly, remInvPoly)
	ntt.MulTo(resolPoly, resolPoly, remInvPoly)
	ntt.InvNTTTo(resolPoly, resolPoly)

	// Extract the meaningful values.
	resol := make([]uint64, rank)
	gen := num.Generators(uint64(params.CycloIdx()))[0]
	idx := uint64(1)
	for i := 0; i < rank; i++ {
		resol[i] = resolPoly.Coeffs[0][idx] + modulus.Value() - resolPoly.Coeffs[0][0]
		if resol[i] >= modulus.Value() {
			resol[i] -= modulus.Value()
		}
		idx = num.Mul(idx, gen, params.CycloIdxMod())
	}

	return resol
}
