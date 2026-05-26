package crt

// isEqualType checks if given elements are of the same type.
func isEqualType(e0, e1 *Element, t ElementType) bool {
	return e0.Type() == t && e1.Type() == t
}

// orderByType returns e0, e1 as the order of scalar and poly.
// Assumes that one of e0, e1 is scalar and the other is poly.
func orderByType(e0, e1 *Element) (c *Element, p *Element) {
	if e0.Type() == TypeScalar {
		return e0, e1
	}
	return e1, e0
}
