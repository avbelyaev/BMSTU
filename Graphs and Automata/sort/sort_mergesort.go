//сортировка слиянием
package main
import "fmt"

var i, n, x int
var arr []int

func mergesort (first int, last int, ch chan int) {

	if first >= last {
		ch <- arr[first]
		close(ch)
	} else {
		var a, b, t_a, t_b, mid int = 0, 0, 0, 0, (first + last) / 2
		var chan1, chan2 = make(chan int), make(chan int)

		go mergesort(first, mid, chan1)
		go mergesort(mid + 1, last, chan2)

		a, ok1 := <-chan1
		b, ok2 := <-chan2

		for (ok1 || ok2) && (p != 3) {
			if ok1 && ok2 {
				if a < 0 { t_a = -a } else { t_a = a }
				if b < 0 { t_b = -b } else { t_b = b }
				if t_a > t_b {
					ch <- b
					b, ok2 = <-chan2
				} else {
					ch <- a
					a, ok1 = <-chan1
				}
			}
                        if ok1 && !ok2 {
				ch <- a
				a, ok1 = <-chan1
			}
			if ok2 && !ok1 {
				ch <- b
				b, ok2 = <-chan2
			}	
		}
		close(ch)
	}
}

func main() {
	input.Scanf("%d", &n)
	arr = make([]int, n)
	for i = 0; i < n; i++ { input.Scanf("%d", &arr[i]) }
	if 1 == n {
		fmt.Printf("%d", arr[0])
		return
	}
        ch := make(chan int)
	go mergesort(0, n - 1, ch)
	
	for i = 0; i < n; i++ { fmt.Printf("%d ", <-ch) }
}

