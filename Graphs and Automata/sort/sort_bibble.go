package main

import "fmt"
import "strconv"

func main() {
        var bubble func(n int, a []int)
	bubble = func(n int, a []int) {
		var bound, i, t int = 0, 0, n - 1

		var combo func(a, b int) int
		combo = func(a, b int) int {
			abt := strconv.Itoa(a) + strconv.Itoa(b)
			bat := strconv.Itoa(b) + strconv.Itoa(a)
			ab, _ := strconv.Atoi(abt)
			ba, _ := strconv.Atoi(bat)
			if ab > ba {
				return 1
			} else {
				return 0
			}
			return 5
		}
		for t > 0 {
			bound, t, i = t, 0, 0
			for i < bound {
				if 0 != combo(a[i+1], a[i]) {
					a[i+1], a[i] = a[i], a[i+1]
					t = i
				}
				i++
			}
		}
		return
	}

	var n, x, i int
	fmt.Scanf("%d", &n)
	i = 0
	a := make([]int, 0, 0)

	for i < n {
		fmt.Scanf("%d", &x)
		a = append(a, x)
		i++
	}

	bubble(n, a)

	for _, x := range a {
		fmt.Printf("%d", x)
	}
}

