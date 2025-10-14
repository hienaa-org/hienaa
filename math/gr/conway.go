package gr

import (
	"bufio"
	"bytes"
	"compress/gzip"
	_ "embed"
	"strconv"
)

//go:embed conway.gz
var conwayData []byte

// findConway finds the Conway polynomial for given prime p and rank r.
func findConway(p uint64, r int) []int64 {
	gz, err := gzip.NewReader(bytes.NewReader(conwayData))
	if err != nil {
		panic(err)
	}
	defer gz.Close()

	sc := bufio.NewScanner(gz)
	for sc.Scan() {
		data := bytes.Split(sc.Bytes(), []byte(","))
		prime, err := strconv.ParseUint(string(data[0]), 10, 64)
		if err != nil {
			panic(err)
		}
		if prime != p {
			continue
		}

		rank, err := strconv.Atoi(string(data[1]))
		if err != nil {
			panic(err)
		}
		if rank != r {
			continue
		}

		poly := make([]int64, 0)
		for _, c := range data[2:] {
			ci, err := strconv.ParseInt(string(c), 10, 64)
			if err != nil {
				panic(err)
			}
			poly = append(poly, ci)
		}
		return poly
	}

	panic("findConway: conway polynomial not found")
}
