package main
import "fmt"

type state struct {
        num_of_hits int
}

func main() {
	var N, K, q_start, hit int
	fmt.Scanf("%d\n", &K)
	alphabet := make([]string, K)
	for i := 0; i < K; i++ { fmt.Scanf("%s\n", &alphabet[i]) }
	fmt.Scanf("%d\n", &N)
	transit_matrix, final_state, matrix := make([][]int, N), make([]int, N), make([][]state, N)
	for i := 0; i < N; i++ { matrix[i] = make([]state, N); for j := 0; j < N; j++ { matrix[i][j].num_of_hits = 0 } }
	for i := 0; i < N; i++ { transit_matrix[i] = make([]int, K) }
	for i := 0; i < N; i++ { for j := 0; j < K; j++ { fmt.Scanf("%d", &transit_matrix[i][j]) }; fmt.Scanf("\n") }
	for i := 0; i < N; i++ { fmt.Scanf("%d\n", &final_state[i]) }
	fmt.Scanf("%d", &q_start)
	fmt.Printf("digraph {\n\trankdir = LR\n\tdummy [label = \"\", shape = none]\n")
	for i := 0; i < N; i++ { if 1 == final_state[i] { fmt.Printf("\t%d [shape = doublecircle]\n", i) } else { fmt.Printf("\t%d [shape = circle]\n", i) } }
	fmt.Printf("\tdummy -> %d\n", q_start)
	for i := 0; i < N; i++ { for j := 0; j < N; j++ { for k := 0; k < K; k++ { if j == transit_matrix[i][k] { matrix[i][j].num_of_hits++ } } } }
	for i := 0; i < N; i++ {
		hit = 0
		for j := 0; j < N; j++ {
			if 0 != matrix[i][j].num_of_hits {
				fmt.Printf("\t%d -> %d [label = \"%s", i, j, alphabet[hit])
				for k := hit + 1; k < matrix[i][j].num_of_hits; k++ { fmt.Printf(", %s", alphabet[k]) }
				fmt.Printf("\"]\n")
				hit++
			}
		}
	}
}

/*in:
2
a b
4
1 3
1 2
3 3
3 3
0
0
1
0
0
*/
/*out:
digraph { 
        rankdir = LR 
        dummy [label = "", shape = none] 
        0 [shape = circle] 
        1 [shape = circle] 
        2 [shape = doublecircle] 
        3 [shape = circle] 
        dummy -> 0 
        0 -> 1 [label = "a"] 
        0 -> 3 [label = "b"] 
        1 -> 1 [label = "a"] 
        1 -> 2 [label = "b"] 
        2 -> 3 [label = "a, b"] 
        3 -> 3 [label = "a, b"] 
}
*/

