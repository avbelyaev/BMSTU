package main
import "fmt"

type state struct {
	in, out int
	char    string
}
type elem struct {
	num_of_hits int
	chars       []string
}

func main() {
	var N, M, q_start int
	fmt.Scanf("%d\n%d\n", &N, &M)
	transit, final_state, matrix := make([]state, M), make([]int, N), make([][]elem, N)
	for i := 0; i < N; i++ { matrix[i] = make([]elem, N); for j := 0; j < N; j++ { matrix[i][j].num_of_hits = 0 } }
	for i := 0; i < M; i++ { fmt.Scanf("%d%d%s\n", &transit[i].in, &transit[i].out, &transit[i].char) }
	for i := 0; i < N; i++ { fmt.Scanf("%d\n", &final_state[i]) }
	fmt.Scanf("%d", &q_start)
	fmt.Printf("digraph {\n\trankdir = LR\n\tdummy [label = \"\", shape = none]\n")
	for i := 0; i < N; i++ { if 1 == final_state[i] { fmt.Printf("\t%d [shape = doublecircle]\n", i) } else { fmt.Printf("\t%d [shape = circle]\n", i) } }
	fmt.Printf("\tdummy -> %d\n", q_start)
	for k := 0; k < M; k++ {
		for j := 0; j < N; j++ {
			if j == transit[k].out {
				matrix[transit[k].in][j].num_of_hits++
				if transit[k].char == "lambda" { matrix[transit[k].in][j].chars = append(matrix[transit[k].in][j].chars, "λ") } else { matrix[transit[k].in][j].chars = append(matrix[transit[k].in][j].chars, transit[k].char) }
			}
		}
	}
	for i := 0; i < N; i++ {
		for j := 0; j < N && 0 != matrix[i][j].num_of_hits; j++ {
			if 0 != matrix[i][j].num_of_hits {
				fmt.Printf("\t%d -> %d [label = \"%s", i, j, matrix[i][j].chars[0])
				for k := 1; k < len(matrix[i][j].chars); k++ { fmt.Printf(", %s", matrix[i][j].chars[k]) }
				fmt.Printf("\"]\n")
			}
		}
	}
}

/*in:
5 
8 
0 1 a 
0 2 a 
1 3 b 
3 1 lambda 
3 3 a 
2 4 c 
4 4 a 
4 4 b 
0
0
0
1
1
0
*/
/*out:
digraph { 
        rankdir = LR 
        dummy [label = "", shape = none] 
        0 [shape = circle] 
        1 [shape = circle] 
        2 [shape = circle] 
        3 [shape = doublecircle] 
        4 [shape = doublecircle] 
        dummy -> 0 
        0 -> 1 [label = "a"] 
        0 -> 2 [label = "a"] 
        1 -> 3 [label = "b"] 
        2 -> 4 [label = "c"] 
        3 -> 3 [label = "a"] 
        3 -> 1 [label = "λ"] 
        4 -> 4 [label = "a, b"] 
}
*/

