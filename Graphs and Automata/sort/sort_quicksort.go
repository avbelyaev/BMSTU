package main

import "fmt"

func partition(less func(i, j int) bool, swap func(i, j int), low, high int) int {
        i, j := low, low

	for j < high {
		if less(j, high) {
			swap(i, j)
			i++
		}
		j++
	}
	swap(i, high)
	return i
}

func qsort(n int, less func(i, j int) bool, swap func(i, j int)) {
	var qsortrec func(less func(i, j int) bool, swap func(i, j int), low int, high int)
	qsortrec = func(less func(i, j int) bool, swap func(i, j int), low int, high int) {
		if low < high {
			q := partition(less, swap, low, high)
			qsortrec(less, swap, low, q - 1)
			qsortrec(less, swap, q + 1, high)
		}
	}
	qsortrec(less, swap, 0, n - 1)
}

func main() {
	var i, n, b int
	fmt.Scanf("%d", &n)
	a := make([]int, 0, n)

	for i < n {
		fmt.Scanf("%d", &b)
		a = append(a, b)
		i++
	}

	qsort(len(a), func(i, j int) bool { return a[i] < a[j] }, func(i, j int) { a[i], a[j] = a[j], a[i] },)

	for _, b := range a {
		fmt.Printf("%d ", b)
	}
}

