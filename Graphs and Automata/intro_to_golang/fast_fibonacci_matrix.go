package main
import "fmt"
import "math/big"
			   
func fib(n int) *big.Int {
        mul := func(tl1, tr1, bl1, br1, tl2, tr2, bl2, br2 *big.Int) (a, b, c, d *big.Int) {
		a, b, c, d = new(big.Int), new(big.Int), new(big.Int), new(big.Int)
		a1, a2 := new(big.Int), new(big.Int)
		a.Add(a1.Mul(tl1, tl2), a2.Mul(tr1, bl2)) // a = (topleft1 * topleft2) + (topright1 * botleft2)
		b.Add(a1.Mul(tl1, tr2), a2.Mul(tr1, br2))
		c.Add(a1.Mul(bl1, tl2), a2.Mul(br1, bl2))
		d.Add(a1.Mul(bl1, tr2), a2.Mul(br1, br2))
		return
	}
	//   |A B| * |E F|
	//   |C D|   |G H|
	a, b, c, d := big.NewInt(0), big.NewInt(1), big.NewInt(1), big.NewInt(1)
	e, f, g, h := big.NewInt(1), big.NewInt(0), big.NewInt(0), big.NewInt(1)
	for n != 0 {
		if (n & 1) != 0 { e, f, g, h = mul(e, f, g, h, a, b, c, d) }
		a, b, c, d = mul(a, b, c, d, a, b, c, d)
		n >>= 1
	}
	return f  // f = topright2
}

func main() {
	var n int
	fmt.Scanf("%d", &n)
	res := new(big.Int)
	fmt.Pritf("%d", fib(n))
}
//
// fib(1'000'000) has 208,988 decimals

