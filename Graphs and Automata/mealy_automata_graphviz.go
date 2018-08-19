package main
import "fmt"

type state struct {
	in  int
	out string
}

func main() {
	var n, m, q_start int
	fmt.Scanf("%d\n%d\n%d\n", &n, &m, &q_start)
	M := make([][]state, n)
	for i := 0; i < n; i++ { M[i] = make([]state, m) }
	for i := 0; i < n; i++ { for j := 0; j < m; j++ { fmt.Scanf("%d", &M[i][j].in) }; fmt.Scanf("\n") }
	for i := 0; i < n; i++ { for j := 0; j < m; j++ { fmt.Scanf("%s", &M[i][j].out) }; fmt.Scanf("\n") }
	fmt.Printf("digraph {\n\trankdir = LR\n\tdummy [label = \"\", shape = none]\n")
	for i := 0; i < n; i++ { fmt.Printf("\t%d [shape = circle]\n", i) }
	fmt.Printf("\tdummy -> %d\n", q_start)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ { fmt.Printf("\t%d -> %d [label = \"%c(%s)\"]\n", i, M[i][j].in, 'a'+j, M[i][j].out) }
	}
	fmt.Printf("}")
}
/*
in:
4
3
0
1 3 3
1 1 2
2 2 2
1 2 3
x y y
y y x
x x x
x y y
*/
/*
out:
digraph {
    rankdir = LR
    dummy [label = "", shape = none]
    0 [shape = circle]
    1 [shape = circle]
    2 [shape = circle]
    3 [shape = circle]
    dummy -> 0
    0 -> 1 [label = "a(x)"]
    0 -> 3 [label = "b(y)"]
    0 -> 3 [label = "c(y)"]
    1 -> 1 [label = "a(y)"]
    1 -> 1 [label = "b(y)"]
    1 -> 2 [label = "c(x)"]
    2 -> 2 [label = "a(x)"]
    2 -> 2 [label = "b(x)"]
    2 -> 2 [label = "c(x)"]
    3 -> 1 [label = "a(x)"]
    3 -> 2 [label = "b(y)"]
    3 -> 3 [label = "c(y)"]
}
*/

