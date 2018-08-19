package main
import "fmt"
import "github.com/skorobogatov/input"

var n, m, q_start, time int
var h1 []int
var h2 []int

func DFS(q int, delta [][]int) { if 0 > h1[q] { h1[q], time, h2 = time, time+1, append(h2, q); for i := 0; i < len(delta[q]); i++ { DFS(delta[q][i], delta) } } }

func Init_canon_auto(delta [][]int, phi [][]string) (delta1 [][]int, phi1 [][]string) {
	delta1, phi1 = make([][]int, time), make([][]string, time)
	for i := 0; i < time; i++ { delta1[i], phi1[i] = make([]int, m), make([]string, m) }
	for i := 0; i < time; i++ {
		x, curr1, curr2 := h2[i], -1, "x"
		for j := 0; j < len(delta[x]); j++ {
			curr1 = delta[x][j]
			delta1[i][j] = h1[curr1]
		}
		for j := 0; j < len(phi[x]); j++ {
			curr2 = phi[x][j]
			phi1[i][j] = curr2
		}
	}
	return
}

func Visualization(delta1 [][]int, phi1 [][]string) {
	fmt.Printf("%d\n%d\n0\n", time, m)
	for i := 0; i < time; i++ { for j := 0; j < m; j++ { fmt.Printf("%d ", delta1[i][j]) }; fmt.Printf("\n") }
	for i := 0; i < time; i++ { for j := 0; j < m; j++ { fmt.Printf("%s ", phi1[i][j]) }; fmt.Printf("\n") }
}

func main() {
	input.Scanf("%d%d%d", &n, &m, &q_start)
	h1, h2, time = make([]int, n), make([]int, 0, n), 0
	delta, phi := make([][]int, n), make([][]string, n)
	for i := 0; i < n; i++ {
		delta[i], h1[i] = make([]int, m), -1
		for j := 0; j < m; j++ { input.Scanf("%d", &delta[i][j]) }
	}
	for i := 0; i < n; i++ {
		phi[i] = make([]string, m)
		for j := 0; j < m; j++ { input.Scanf("%s", &phi[i][j]) }
	}

	DFS(q_start, delta)

	m = len(delta[q_start])
	new_delta, new_phi := Init_canon_auto(delta, phi)

	Visualization(new_delta, new_phi)
}

