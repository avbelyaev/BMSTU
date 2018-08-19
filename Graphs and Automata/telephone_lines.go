package main
import "github.com/skorobogatov/input"
import "fmt"
import "os"

type queue struct {
        heap       []*vertex
	cap, count int
}

type vertex struct {
	index, key int
	value      *vertex
	p          *edge
}

type edge struct {
	u, v, len int
	p         *edge
}
var n int
var q queue

func InitPriorityQueue(size int) {
	q.heap = make([]*vertex, size)
	q.cap, q.count = size, 0
}

func Insert(ptr *vertex) {
	i := q.count
	if i == q.cap { fmt.Printf("overflow\n") }
	q.count++
	q.heap[i] = ptr
	for i > 0 && q.heap[(i-1)/2].key > q.heap[i].key {
		q.heap[(i-1)/2], q.heap[i] = q.heap[i], q.heap[(i-1)/2]
		q.heap[i].index = i
		i = (i - 1) / 2
	}
	q.heap[i].index = i
}

func ExtractMin() *vertex {
	var Heapify func(i int, k int, P []*vertex)
	Heapify = func(i int, k int, P []*vertex) {
		for {
			l := 2*i + 1
			r := l + 1
			j := i
			if l < k && P[i].key > P[l].key { i = l }
			if r < k && P[i].key > P[r].key { i = r }
			if i == j { break }
			P[i], P[j] = P[j], P[i]
			P[i].index = i
			P[j].index = j
		}
	}
	if 0 == q.count { fmt.Printf("Queue is empty\n") }
	ptr := q.heap[0]
	q.count--
	if 0 < q.count {
		q.heap[0] = q.heap[q.count]
		q.heap[0].index = 0
		Heapify(0, q.count, q.heap)
	}
	return ptr
}

func DecreaseKey(i int, k int) {
	q.heap[i].key = k
	for i > 0 && q.heap[(i-1)/2].key > k {
		q.heap[(i-1)/2], q.heap[i] = q.heap[i], q.heap[(i-1)/2]
		q.heap[i].index = i
		i = (i - 1) / 2
	}
	q.heap[i].index = i
}

func nless2(num int) bool {
	var a1, b1, c1, a2, b2, c2 int
	if 0 == num {
		fmt.Println(c1)
		return false
	} else if 1 == num {
		input.Scanf("%d %d %d\n", &a1, &b1, &c1)
		fmt.Printf("%d\n", c1)
		return false
	}
	if 2 == num {
		input.Scanf("%d %d %d\n", &a1, &b1, &c1)
		input.Scanf("%d %d %d\n", &a2, &b2, &c2)
		fmt.Printf("%d\n", c1+c2)
		return false
	}
	return true
}

func main() {
	var u, v, len, n, m, i, res int
	Mark := func(G []vertex, size int) { for i := 0; i < n; i++ { G[i].index = -1 } }
	input.Scanf("%d\n", &n)
	input.Scanf("%d\n", &m)
	if k := nless2(m); !(k) { os.Exit(0) }
	Graph := make([]vertex, n)
	for i < m {
                var x, y edge
		input.Scanf("%d %d %d\n", &v, &u, &len)
		x.u, x.v, x.len, x.p = u, v, len, Graph[v].p
		y.u, y.v, y.len, y.p = v, u, len, Graph[u].p
		Graph[u].p = &y
		Graph[v].p = &x
		i++
	}
	InitPriorityQueue(n - 1)
	Mark(Graph, n)
	vert := &Graph[0]
	for {
		vert.index = -2
		instant := vert.p
		for instant != nil {
			if -1 == Graph[instant.u].index {
				Graph[instant.u].key = instant.len
				Graph[instant.u].value = vert
				Insert(&Graph[instant.u])
			} else if -2 != Graph[instant.u].index && instant.len < Graph[instant.u].key {
				Graph[instant.u].value = vert
				DecreaseKey(Graph[instant.u].index, instant.len)
			}
			instant = instant.p
		}
		if 0 == q.count { break }
		vert = ExtractMin()
		res += vert.key
	}
	fmt.Printf("%d", res)
}

