package rns

var (
	// cyclicNTTFactors are the factors of the degree for native cyclic NTT.
	cyclicNTTFactors = []uint64{2, 3, 5}
	// bluesteinFactos are the factors of the ambient degree used in Bluestein NTT.
	ambientDegreeFactors = []uint64{2, 3}
)
