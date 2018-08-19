package main
import "fmt"
//import "github.com/skorobogatov/input"

type vertex struct {
        x      int
	mark   string
	comp   int
	edge   []*vertex
	parent *vertex
}

func main() {
	var i, n, m, v, u, comp_dfs1, comp_dfs2, t_c int
	var VisitVertex2 func(v *vertex, component int)
	var VisitVertex1 func(v *vertex) int
	var queue [](*vertex)

	fmt.Scanf("%d\n", &n)
	fmt.Scanf("%d\n", &m)
	
	list := make([]vertex, n)

	for i < m {
		fmt.Scanf("%d %d\n", &u, &v)
		list[u].edge, list[v].edge = append(list[u].edge, &list[v]), append(list[v].edge, &list[u])
		i++
	}

	Enqueue := func(Q [](*vertex), vert *vertex) []*vertex { return append(Q, vert) }
	Dequeue := func(Q [](*vertex), num int) *vertex { return Q[num] }
	
	VisitVertex1 = func(v *vertex) {
		E_G := v.edge
		v.mark = "gray"
		queue = Enqueue(queue, v)
		for i, _ := range E_G {
			u := E_G[i]
			if u.mark == "white" {
				u.parent = v
				VisitVertex1(u)
			}
		}
		v.mark = "black"
	}

	VisitVertex2 = func(v *vertex, component int) {
		v.comp = component
		E_G := v.edge
		for i, _ := range E_G {
			u := E_G[i]
			if u.comp == -1 && u.parent != v { VisitVertex2(u, component) }
		}
	}
	V_G := list
//DFS1
	for i, _ := range V_G { V_G[i].mark = "white" }
	for i, _ := range V_G {
		v := V_G[i]
		if v.mark == "white" { 
			VisitVertex1(&v)
 			comp_dfs1++ 
		}
	}
//DFS2
	for i, _ := range V_G { V_G[i].comp = -1 }
	for i, _ := range queue {
		v := Dequeue(queue, i)
		if v.comp == -1 {
			VisitVertex2(v, comp_dfs2)
			comp_dfs2++
		}
	}
	fmt.Printf("comp=%d\n", comp_dfs2-comp_dfs1)
}

