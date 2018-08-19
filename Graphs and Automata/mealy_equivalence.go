package main

import "os"
import "fmt"
import "github.com/skorobogatov/input"

var time int
var h1, h2 []int

func Split(delta [][]int, phi [][]string, pi []int, q0 int) (int, []int) {
	m_classes := len(delta)
	pi1 := make([]int, 0, m_classes)
        var Union func(int, int, int)
	Union = func(x, y, z int) {
		for z < len(pi1) {
			if pi1[z] == x { pi1[z] = y }
			z++
		}
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
				if eq {
					Union(pi1[i], pi1[j], 0)
					m_classes--
				}
			}
		}
	}
	return m_classes, pi1
}

func DFS(q int, delta [][]int, phi [][]string) {
	if 0 > h1[q] {
		h1[q], time, h2 = time, time+1, append(h2, q)
		for i := 0; i < len(delta[q]); i++ { DFS(delta[q][i], delta, phi) }
	}
}

func Init_canon_auto(delta [][]int, phi [][]string, q0 int) (delta1 [][]int, phi1 [][]string, yes bool) {
	delta1, phi1, yes = make([][]int, time), make([][]string, time), false
	for i := 0; i < time; i++ { delta1[i], phi1[i] = make([]int, len(delta[0])), make([]string, len(delta[q0])) }
	for i := 0; i < time; i++ {
		yes = true
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

func AufenkampHohn(delta [][]int, phi [][]string, q0 int) ([][]int, [][]string) {
//Split1 {
	m := len(delta)
	pi := make([]int, 0, m)
	for i, _ := range delta { pi = append(pi, i) }
	for i, _ := range delta {
		for j, _ := range delta {
			if pi[i] != pi[j] {
				eq := true
				for k := 0; k < len(delta[0]); k++ {
					if phi[i][k] != phi[j][k] {
						eq = false
						break
					}
				}
				if eq {
					for q := 0; q < len(pi); q++ { if pi[q] == pi[i] { pi[q] = pi[j] } }
					m--
				}
			}
		}
	}
// } Split1
	for {
		m1, pi1 := Split(delta, phi, pi, q0)
		if m == m1 { break }
		m, pi = m1, pi1
	}
	h3 := make(map[int]int, len(pi))
	
	var Exists func(int) bool
	Exists = func(v int) bool { _, ok := h3[v]; if ok { return true }; return false }

	for i, x := range pi {
		if Exists(x) { y := h3[x]; pi[i] = y } else { pi[i], h3[x] = len(h3), len(h3) }
	}

	delta1, phi1, q01 := make([][]int, m), make([][]string, m), 0
	for i := 0; i < m; i++ { delta1[i], phi1[i] = make([]int, len(delta[q0])), make([]string, len(delta[q0]) )}
	
	for i := 0; i < len(delta); i++ {
		q := pi[i]
		if i == q0 { q01 = q }
		for j := 0; j < len(delta[q0]); j++ { delta1[q][j], phi1[q][j] = pi[delta[i][j]], phi[i][j] }
	}

	DFS(q01, delta1, phi1)
	
	new_delta, new_phi, ok := Init_canon_auto(delta1, phi1, 0)
	
	if ok { return new_delta, new_phi }; return nil, nil
}

func Read(n, m int) ([][]int, [][]string) {
	a, b := make([][]int, n), make([][]string, n)
	h1, h2, time = make([]int, n), make([]int, 0, n), 0
	for i := 0; i < n; i++ { a[i], b[i], h1[i] = make([]int, m), make([]string, m), -1 }
	for i := 0; i < n; i++ { for j := 0; j < m; j++ { input.Scanf("%d", &a[i][j]) } }
	for i := 0; i < n; i++ { for j := 0; j < m; j++ { input.Scanf("%s", &b[i][j]) } }
	return a, b
}

func main() {
	var n2, m2, q_start2, mark_delta, mark_phi int = 0,0,0,0,0
//Auto 1
	input.Scanf("%d\n%d\n%d\n", &n2, &m2, &q_start2)
	delta2, phi2 := Read(n2, m2)
	minimized_delta2, minimized_phi2 := AufenkampHohn(delta2, phi2, q_start2)
//Auto 2
	var n3, m3, q_start3 int
	input.Scanf("%d\n%d\n%d\n", &n3, &m3, &q_start3)
	delta3, phi3 := Read(n3, m3)
	minimized_delta3, minimized_phi3 := AufenkampHohn(delta3, phi3, q_start3)
//Compare
	if len(minimized_delta2) != len(minimized_delta3 || len(minimized_delta2[0]) != len(minimized_delta3[0]) { fmt.Printf("NOT EQUAL\n"); os.Exit(1) }
	for i := 0; i < len(minimized_delta2); i++ { for j := 0; j < len(minimized_delta2[0]); j++ { if minimized_delta2[i][j] == minimized_delta3[i][j] { mark_delta++ } } }
	for i := 0; i < len(minimized_phi2); i++ { for j := 0; j < len(minimized_phi2[0]); j++ { if minimized_phi2[i][j] == minimized_phi3[i][j] { mark_phi++ } } }
	if mark_delta == len(minimized_delta2)*len(minimized_delta2[0]) && mark_phi == len(minimized_phi2)*len(minimized_phi2[0]) { fmt.Printf("EQUAL\n") } else { fmt.Printf("NOT EQUAL\n") }
}

