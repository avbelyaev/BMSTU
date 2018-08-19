package main
import "fmt"

var k, m, n, q0, k1, m1, n1, q01 int

func new_delta(n, m int) [][]int {
        a := make([][]int, n)
	for i := 0; i < n; i++ { a[i] = make([]int, m) }
	return a
}

func new_str_matrix(n, m int) [][]string {
	a := make([][]string, n)
	for i := 0; i < n; i++ { a[i] = make([]string, m) }
	for i := 0; i < n; i++ { for j := 0; j < m; j++ { a[i][j] = "" } }
	return a
}

func Split1_Recognizer(F []int) (m int, pi []int) {
	pi, m = make([]int, 0), 2
	for i := 0; i < n; i++ { pi = append(pi, 0) }
	for q, _ := range F { if F[q] == 1 { pi[q] = 1 } else { pi[q] = 0 } }
	return
}

func Split(delta [][]int, pi []int) (int, []int) {
	var Union func(int, int, int)
	m_classes, pi1 := len(delta), make([]int, 0, len(delta))
	
	Union = func(x, y, z int) {
		for z < len(pi1) { if pi1[z] == x { pi1[z] = y }; z++ }
	}

	for i, _ := range delta { pi1 = append(pi1, i) }
	for i, _ := range delta {
		for j, _ := range delta {
			if pi[i] == pi[j] && pi1[i] != pi1[j] {
				eq := true
				for k := 0; k < len(delta[q0]); k++ {
					w1, w2 := delta[i][k], delta[j][k]
					if pi[w1] != pi[w2] {
						eq = false
						break
					}
				}
				if eq { Union(pi1[i], pi1[j], 0); m_classes-- }
			}
		}
	}
	return m_classes, pi1
}

func AufenkampHohn_Recognizer(delta [][]int, F []int) ([][]int, []int) {
        m, pi := Split1_Recognizer(F)
	m1, pi1 := -1, make([]int, 0, n)
	for {
		if m1, pi1 = Split(delta, pi); m == m1 { break }
		m, pi = m1, pi1
	}
	m1 = 0
	pi2, h3 := new_delta(n, n), make(map[int]int)

	var Exists func(int) bool
	Exists = func(v int) bool {
		_, ok := h3[v]
		if ok { return true }
		return false
	}

	for i := 0; i < n; i++ {
		if !Exists(pi1[i]) {
			h3[pi1[i]] = m1
			for j := 0; j < n; j++ { if pi1[i] == pi1[j] { pi2[j][0] = h3[pi1[i]] } }
			m1++
		}
	}

	delta1, F1 := new_delta(m, k), make([]int, m)
	for i := 0; i < n; i++ {
		for j := 0; j < k; j++ { delta1[pi2[i][0]][j] = pi2[delta[i][j]][0] }
		F1[pi2[i][0]] = F[i]
	}

	k1, n1, q01 = k, m, pi2[q0][0]
	return delta1, F1
}

func Vis_Recognizer(alphabet []string, delta1 [][]int, F1 []int) {
	matr := new_str_matrix(n1, n1)
	fmt.Printf("digraph {\n\trankdir = LR\n\tdummy [label = \"\", shape = none]\n")
	for i := 0; i < n1; i++ {
		fmt.Printf("\t%d [shape = ", i)
		if F1[i] == 1 { fmt.Printf("double") }
		fmt.Printf("circle]\n")
	}
	fmt.Printf("\tdummy -> %d\n", q01)

	for i := 0; i < n1; i++ { for j, x := range delta1[i] { matr[i][x] += alphabet[j]+"," } }
	for i := 0; i < n1; i++ { for j := 0; j < n1; j++ {
		if matr[i][j] != "" { 
                        fmt.Printf("\t%d -> %d [label = \"%s\"]\n", i, j, matr[i][j]) } } }
	fmt.Printf("}\n")
}

func Read() ([]string, [][]int, []int) {
	fmt.Scanf("%d\n", &k)
	alphabet := make([]string, k)
	for i := 0; i < k; i++ { fmt.Scanf("%s", &alphabet[i]) }

	fmt.Scanf("\n%d\n", &n)
	delta := new_delta(n, k)
	for i := 0; i < n; i++ { for j := 0; j < k; j++ { fmt.Scanf("%d", &delta[i][j]) }; fmt.Scanf("\n") }
	F := make([]int, n)

	for i := 0; i < n; i++ { fmt.Scanf("%d", &F[i]) }
	fmt.Scanf("\n%d", &q0)
	return alphabet, delta, F
}

func main() {
	alphabet, delta, F := Read()

	delta1, F1 := new_delta(n, k), make([]int, 0)
	delta1, F1 = AufenkampHohn_Recognizer(delta, F)

	Vis_Recognizer(alphabet, delta1, F1)
}

