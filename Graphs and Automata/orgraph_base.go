package main
import "github.com/skorobogatov/input"
import "fmt"
import "math"
import "os"

type vertex struct {
        posx, indx, comp, T1, low int
	next                      *vertex
}

type edge struct {
        x int
	u, v *vertex
	p    *vertex
	next *edge
}

type stack struct {
	data       []int
	top        interface{}
	cap, count int
}

var n, m, time, count int
var G []vertex
var E []edge

func abs(x int) int { 
        if x < 0 { return -x } 
        return x  
}

func InitStack(size int) *stack {
	var s stack
	s.data = make([]int, size)
	s.cap = size
	s.count = 0
	s.top = nil
	return &s
}

func (s *stack) StackEmpty() bool { return s.top == 0 }

func (s *stack) Push(value int) {
	if s.top == s.cap { fmt.Printf("overflow\n") }
	s.data = append(s.data, value)
	s.count++
}

func (s *stack) Pop() interface{} {
	if s.StackEmpty() { fmt.Printf("devastation\n") }
	ret := s.data[len(s.data)-1]
	s.data = s.data[0 : len(s.data)-1]
	s.count--
	return ret
}

func Mark(a []int, b [][]int, c []vertex) { if 0 != len(a) && 0 != len(b) { for i := range c { c[i].comp, c[i].T1 = 0, 0 } } else { for i := range b { for j := range b { b[i][j] = (int)(math.Inf(1)) } } } }

func (s *stack) VisitVertex_Tarjan(v *vertex) {
	v.low, v.T1 = time, time
	time++
	s.Push(v.posx)
	for x := v.next; x != nil; x = x.next {
		if G[x.indx].posx = x.indx; 0 == G[x.indx].T1 { s.VisitVertex_Tarjan(&G[x.indx]) }
		if 0 == G[x.indx].comp && v.low > G[x.indx].low { v.low = G[x.indx].low }
	}
	if v.T1 == v.low {
		u := s.Pop().(int)
		G[u].comp = count
		for u != v.posx {
			u = s.Pop().(int)
			G[u].comp = count
		}
		count++
	}
}

func nless2(num int) bool {
        if 0 == num { os.Exit(1) }
	if 1 == num { fmt.Printf("0\n"); return false
	} else { return true }
        return true
}

func main() {
	var u, v, f int
        time, count = 1, 1
	input.Scanf("%d\n%d\n", &n, &m)
        if z := nless2(n); !(z) { os.Exit(0) }
	stack := InitStack(f)
	tmp1 := make([]int, 0)
	tmp2 := make([][]int, 0)
	res := make([]int, n+m)
        c := make([][]int, 0)
	G = make([]vertex, n)
	for i := 0; i < m; i++ {
		input.Scanf("%d %d\n", &u, &v)
		if G[u].next == nil {
			G[u].next = new(vertex)
			G[u].next.indx, G[u].next.next = v, nil
			continue
		}
		next := G[u].next
		for next.next != nil { next = next.next }
		next.next = new(vertex)
	}
	Mark(tmp1, tmp2, G)
	for i := range G { if G[i].posx = i; 0 == G[i].T1 { stack.VisitVertex_Tarjan(&G[i]) } }
	for i := range G {  G[i].comp, E[i].x = G[i].comp-1, E[i].x-1 }
	for i := 0; i < count-1; i++ {
		p := make([]int, count-1)
		c = append(c, p)
		Mark(tmp1, c, G)
	}
	for i := range tmp1 { tmp1[i] = 0 }
	for i := 0; i < n; i++ {
		a := G[i].comp
		for v := G[i].next; v != nil; v = v.next { if a != G[v.indx].comp { c[a][G[v.indx].comp] = 1 } }
	}
	for i := 0; i < abs(count-1); i++ {
		r, k := 0, 0
		for j := 0; j < abs(count-1); j++ { if (int)(math.Inf(1)) != c[j][i] { r += c[j][i] } }
		if 0 == r {
			for k < n {
				if G[k].comp == i {
					res[f], f = k, f+1
					break
				}
				k++
			}
		}
	}
	for i := 0; i < f; i++ { fmt.Printf("%d ", res[i]) }
}

