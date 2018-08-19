package main
import "github.com/skorobogatov/input"
import "fmt"
import "os"
import "math"
type vertex struct {
        value, dist float64
        index, x, y int
	edge        []*vertex
        parent      *vertex
}
type edge struct {
	u, v, ind int
	vertex    *vertex
	weight    int
}
type queue struct {
	heap       []*vertex
	count, cap int
}

func InitPriorityQueue(size int) *queue {
	var q queue
	q.heap, q.count, q.cap = make([]*vertex, size), 0, size
	return &q
}

func (q *queue) InitAttributes(G [][]vertex, n int) {
        for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			v := &G[i][j]
			if v == &G[0][0] { v.dist = v.value } else { v.dist = math.Inf(1) }
			v.parent = nil
			q.QueueInsert(&G[i][j])
		}
	}
}

func InitGraph(size int) [][]vertex {
	graph_map := make([][]vertex, size)
	for i := 0; i < size; i++ { graph_map[i] = make([]vertex, size) }
	return graph_map
}

func (q *queue) QueueInsert(ptr *vertex) {
	i := q.count
	if i == q.cap { fmt.Printf("overflow\n") }
	q.count++
	q.heap[i] = ptr
	for i > 0 && q.heap[i].dist < q.heap[(i-1)/2].dist {
		q.heap[(i-1)/2], q.heap[i] = q.heap[i], q.heap[(i-1)/2]
		q.heap[i].index = i
		i = (i - 1) / 2
	}
	q.heap[i].index = i
}

func (q *queue) ExtractMin() *vertex {
	var Heapify func(i int, n int, P []*vertex)
	Heapify = func(i int, n int, P []*vertex) {
		for {
			l := 2*i + 1
			r := l + 1
			j := i
			if l < n && P[i].dist > P[l].dist { i = l }
			if r < n && P[i].dist > P[r].dist { i = r }
			if i == j { break }
			P[i], P[j] = P[j], P[i]
			P[i].index, P[j].index = i, j
		}
	}
	if 0 == q.count { fmt.Printf("Queue is empty\n") }
	ptr := q.heap[0]
	q.count--
	if 0 != q.count {
		q.heap[0] = q.heap[q.count]
		q.heap[0].index = 0
		Heapify(0, q.count, q.heap)
	}
	return ptr
}

func (q *queue) DecreaseKey(i int, k float64) {
	q.heap[i].dist = k
	for i > 0 && q.heap[(i-1)/2].dist > k {
		q.heap[(i-1)/2], q.heap[i] = q.heap[i], q.heap[(i-1)/2]
		q.heap[i].index = i
		i = (i - 1) / 2
	}
	q.heap[i].index = i
}

func nless2(num int) bool {
	var a, b, c, d int
	if 1 == num {
		input.Scanf("%d\n", &a)
		fmt.Printf("%d", a)
		return false
	}
	if 2 == num {
		input.Scanf("%d %d\n%d %d\n", &a, &b, &c, &d)
		if c > b { fmt.Printf("%d", a+b+d) } else { fmt.Printf("%d", a+c+d) }
		return false
	}
	return true
}

func Relax(u *vertex, v *vertex, weight float64) (changed bool) {
	changed = false
	if u.dist+weight < v.dist { changed = true }
	if changed { v.dist = u.dist + weight }
	return changed
}

func main() {
	var n int
	var x float64
	input.Scanf("%d\n", &n)
	if n > 1200 { goto next }
	if m := nless2(n); !(m) { os.Exit(0) }
next:
	graph_map := InitGraph(n)
	edge_map := make([]edge, n*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			input.Scanf("%f", &x)
			graph_map[i][j].value = x
			graph_map[i][j].x, graph_map[i][j].y = i, j
                        graph_map[i][j].parent = nil
			edge_map[i].weight = (int)(x)
			edge_map[i].ind = i
		}
	}
	queue := InitPriorityQueue(n * n)
	queue.InitAttributes(graph_map, n)
	for 0 != queue.count {
		v := queue.ExtractMin()
		v.index = -1
		x, y := v.x, v.y
		if x+1 < n && -1 != graph_map[x+1][y].index && Relax(v, &graph_map[x+1][y], graph_map[x+1][y].value) { queue.DecreaseKey(graph_map[x+1][y].index, graph_map[x+1][y].dist) }
		if y+1 < n && -1 != graph_map[x][y+1].index && Relax(v, &graph_map[x][y+1], graph_map[x][y+1].value) { queue.DecreaseKey(graph_map[x][y+1].index, graph_map[x][y+1].dist) }
		if x > 0 && -1 != graph_map[x-1][y].index && Relax(v, &graph_map[x-1][y], graph_map[x-1][y].value) { queue.DecreaseKey(graph_map[x-1][y].index, graph_map[x-1][y].dist) }
		if y > 0 && -1 != graph_map[x][y-1].index && Relax(v, &graph_map[x][y-1], graph_map[x][y-1].value) { queue.DecreaseKey(graph_map[x][y-1].index, graph_map[x][y-1].dist) }

	}
	bottom_right := (int)(graph_map[n-1][n-1].dist)
	fmt.Printf("%d", bottom_right)
}

