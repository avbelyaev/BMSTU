package main

import "fmt"
import "os"
//import "github.com/skorobogatov/input"

type vertex struct {
        x      int
	depth  int
	parent *vertex
}

func MakeSet(x int) *vertex {
	var t vertex
	t.x = x
	t.depth = 0
	t.parent = &t
	return &t
}

func Find(x *vertex) *vertex {
	if x.parent == x { return x } 
	return Find(x.parent)
}

func Union(x, y *vertex) {
	rootx, rooty := Find(x), Find(y)
	if rootx.depth < rooty.depth {
		rootx.parent = rooty
	} else {
		rooty.parent = rootx
		if rootx.depth == rooty.depth && rootx != rooty {
			rootx.depth++
		}
	}
}

func get_set(size, k int) [](*vertex) {
	set := make([](*vertex), size)
	for k < size {
		set[k] = MakeSet(k)
		k++
	}
	return set
}

func get_edge(set [](*vertex)) (*vertex, *vertex, bool) {
	var u, v int
	input.Scanf("%d %d\n", &u, &v)
	
	su := set[u]
	sv := set[v]
	du := Find(su)
	dv := Find(sv)
	
	if du != dv { return su, sv, false }
	return su, sv, true
}

func main() {
	var i, n, m, ks int
	input.Scanf("%d\n", &n)
	input.Scanf("%d\n", &m)
	set := get_set(n, i)

	for i < m {
		su, sv, ptr := get_edge(set)
		if Find(su) != Find(sv) && !(ptr) {
			Union(su, sv)
			n--
		}
		i++
	}
	fmt.Printf("%d\n", n)
}

