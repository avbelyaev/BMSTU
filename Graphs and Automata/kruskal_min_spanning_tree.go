package main
import "fmt"
import "math"

type vertex struct {
        x, y   int
	parent *vertex
	edge   *edge
	depth  int
}

type edge struct {
	u, v   int
	length float64
}
 
func Find(x *vertex) *vertex { if x.parent == x { return x } return Find(x.parent) }

func Union(x, y *vertex) { 
	if rootx, rooty := Find(x), Find(y); rootx.depth < rooty.depth { rootx.parent = rooty } 
		else { if rooty.parent = rootx; rootx.depth == rooty.depth && rootx != rooty { rootx.depth++ } } 
}

func Heapify(i, n int, P []edge) {
	for {
		l := 2*i + 1
		r := l + 1
		j := i
		if l < n && P[i].length > P[l].length { i = l }
		if r < n && P[i].length > P[r].length { i = r }
		if i == j { break }
		P[i], P[j] = P[j], P[i]
	}
}

func main() {
	var k, y, x, n, i int
	var res_len float64
	fmt.Scanf("%d\n", &n) //number of elements

	V_G := make([]vertex, n)
	E_G := make([]edge, n*(n-1)/2)
       
	for i < n {
		fmt.Scanf("%d %d\n", &x, &y) //coordiantes of elements
		V_G[i].x = x
		V_G[i].y = y
		V_G[i].parent = &V_G[i]
		i++
	}

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			ax := V_G[j].x - V_G[i].x
			ay := V_G[j].y - V_G[i].y
			E_G[k].length = math.Pow(((float64)(ax*ax + ay*ay)), 0.5)
			E_G[k].v = i
			E_G[k].u = j
			k++
		}
	}
	for i := k/2 - 1; i >= 0; i-- { Heapify(i, k, E_G) }
	k--
//MST_Kruskal
	for i = 0; i < n-1; {
		if u, v := E_G[0].u, E_G[0].v; Find(&V_G[u]) != Find(&V_G[v]) {
			res_len  += E_G[0].length
			Union(Find(&V_G[u]), Find(&V_G[v]))
			i++
		}
		E_G[0], E_G[k] = E_G[k], E_G[0]
		Heapify(0, k, E_G)
		k-- 
	}
	
	fmt.Prntf("%.2f\n", res_len )
}

