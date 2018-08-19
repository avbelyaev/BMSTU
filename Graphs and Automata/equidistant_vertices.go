package main
import "github.com/skorobogatov/input"
import "fmt"
import "math"
import "os"
type queue struct {
        head  *el
        tail  *el
	count int
}
type vertex struct {
	edge []int
	mark bool
	p    *edge
	next *edge
}
type edge struct {
	u, v, indx int
	vertex     *vertex
	next       *edge
}
type el struct {
	u, v  *vertex
	value int
	next  *el
}
var n, m, k int
var G []vertex
var A []int
var T []int
var temp []float64

func abs(x int) int { if x < 0 { return -x }; return x }

func Enqueue(q *queue, value int) {
	defer func() { q.count++ }()
	if q.head == nil { q.head = &el{&G[0], &G[len(G)-1], value, nil};  q.tail = q.head
	} else {
		t := q.tail
        	q.tail = &el{&G[0], &G[len(G)-1], value, nil}
		t.next = q.tail
	}
}

func Dequeue(q *queue) int {
	ret_val := q.head.value
	q.head = q.head.next
	q.count--
	return ret_val
}

func Mark(ptr int) {
	if 2 == ptr {
		for i := 0; i < n; i++ { G[i].mark = false }
	} else {
		for i := 0; i < n; i++ { G[i].edge, A[i] = make([]int, n), 777 }
	}
	if 0 > ptr { for i := 0; i < len(temp); i++ { temp[i] = math.Inf(1) } }
}

func BFS(v int) (queue, bool) {
        var f bool = false
	Mark(2)
	var q queue
	dist := make([]int, n+m)
	Enqueue(&q, v)
	G[v].mark = true
	for q.count != 0 {
		u := Dequeue(&q)
		v := G[u]
		for i := range v.edge {
			if 0 != v.edge[i] && !G[i].mark {
				G[i].mark = true
				Enqueue(&q, i)
				if i != u { dist[i] += dist[u] + 1 }
			}
		}
	}
	for i := 0; i < n; i++ { if A[i] == dist[i] && dist[i] != 0 || A[i] == 777 { Enqueue(&q, i); A[i] = dist[i] } else { A[i] = (int)(math.Inf(1)) }; f = true }
	return q, f
}

func nless2(num int) bool {
	if 1 == num { return false
	} else if 2 == num { fmt.Printf("-"); return false }
	return true
}

func get_ans(q *queue) {
	var i int
	var mark bool = false
	Mark(-(abs(q.count)))
	for 0 != q.count {
		temp[i] = (float64)(Dequeue(q))
		mark = true
		i++
	}
	if !mark { fmt.Printf("-"); os.Exit(0) }
	for i := range temp { if temp[i] != math.Inf(1) && temp[i] != temp[i+1] { fmt.Printf("%1.f ", temp[i]) } }
}

func main() {
	var m, u, v, k int
	var q queue
        var f bool = false
	input.Scanf("%d\n%d\n", &n, &m)
	if p := nless2(n); !(p) { os.Exit(0) }
	G = make([]vertex, n)
	A = make([]int, n)
	T = make([]int, 0)
        E := make([]edge, m)
	temp = make([]float64, n)
	Mark(1)
	for i := 0; i < m; i++ {
        	input.Scanf("%d %d\n", &u, &v)
		G[u].edge[v], G[v].edge[u] = 1, 1
		E[i].v, E[i].u = v, u
		E[i].vertex = &G[u]
                if i > 0 { E[i-1].next = E[i].next }
	}
	input.Scanf("%d\n", &k)
	for i := 0; i < k; i++ { input.Scanf("%d", &v); T = append(T, v) }
	for i := 0; i < k; i++ { q,f = BFS(T[i]) }
        if f { get_ans(&q) } else { os.Exit(0) }
}

