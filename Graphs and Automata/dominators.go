package main
import "fmt"
import "github.com/skorobogatov/input"
import "container/list"  
import "os"

type vertex struct {
	x, indx, time                      int
	dom, sdom, ancestor, label, parent *vertex
	bucket                             []*vertex
	in, out                            []int
}

type stack struct {
	data       []*vertex
	cap, count int
	top        *vertex
}

var G map[int]*vertex
var cycle map[int]bool
var time,rs int = 1, 0

type Stack struct {
     container *list.List
}
func NewStack() *Stack { return new(Stack).Init() }
func (s *Stack) Init() *Stack {
        s.container = list.New()
        return s
}

func DFS(ind int, p *vertex) {
	G[ind].time, G[ind].parent = time, p
	time++
	for _, i := range G[ind].in { if G[i].time == 0 { DFS(i, G[ind]) } }
}

func (s *Stack) Push(value interface{}) { s.container.PushBack(value) }
func (s *Stack) Pop() interface{}       { return s.container.Remove(s.container.Back()) }

func FindMin(v *vertex, n int) (min *vertex) {
        if v.ancestor == nil {
		min = v
	} else {
	        s, u := NewStack(), v
	        for u.ancestor.ancestor != nil { s.Push(u); u = u.ancestor }
	        for 0 != s.container.Len() {
	        	v = s.Pop().(*vertex)
	        	if v.ancestor.label.sdom.time < v.label.sdom.time { v.label = v.ancestor.label }
	        	v.ancestor = u.ancestor
	        }
                min = v.label
	}
	return 
}

func Dominators() {
	temp := make([]*vertex, 0)
	for _, w := range G {
		w.sdom, w.label, w.ancestor = w, w, nil
		if w.time != 0 && nil != w.parent { temp = append(temp, w) }
	}
	last := len(temp) - 1
	for i := 0; i < last; i++ {
		iMin := i
		for j := i + 1; j < len(temp); j++ { if 0 != temp[j].time && 0 != temp[iMin].time && temp[j].time > temp[iMin].time { iMin = j } }
		temp[i], temp[iMin] = temp[iMin], temp[i]
	}
	for _, w := range temp {
		for _, index := range w.out { if u := FindMin(G[index], len(G)); u.sdom.time < w.sdom.time && 0 != G[index].time { w.sdom = u.sdom } }
		w.ancestor = w.parent
		w.sdom.bucket = append(w.sdom.bucket, w)
		for _, v := range w.parent.bucket { if u := FindMin(v, len(G)); u.sdom == v.sdom { v.dom = w.parent } else { v.dom = u } }
		w.parent.bucket = w.parent.bucket[:0]
	}
	for i := 1; i < len(temp); i++ { for j := i; j > 0 && temp[j].time < temp[j-1].time; j-- { temp[j-1], temp[j] = temp[j], temp[j-1] } }
	for i := 0; i < len(temp); i++ { if w := temp[i]; w.dom != w.sdom { w.dom = w.dom.dom } }
}

func nless2 (num int) bool {
        var a int
        if num < 3 {
                if 0 == num { os.Exit(1)
                } else if 1 == num {
                    input.Scanf ("%d\n", &a)
                    fmt.Printf ("0"); os.Exit(0)
                }
        return false
        } else { return true }
        return true
}

func main() {
	var n, label, operand, val, start, prv int
	var command string
	input.Scanf("%d", &n)
	if z:= nless2(n); !(z) { os.Exit(0) }
	cycle = make(map[int]bool)
	G = make(map[int]*vertex)
	tmp := make([]int, n)
	for i := 0; i < n; i++ {
		input.Scanf("%d %s", &label, &command)
		if 0 == i { start = label }
		G[label] = new(vertex)
		G[label].indx = label
		if command == "ACTION" { val = 1 } else if command == "JUMP" { val = 2 } else { val = 3 }
		G[label].x, tmp[i] = val, label
		if 0 != prv  && 2 != G[prv].x { G[prv].in = append(G[prv].in, label) }
		prv = label
		if command == "BRANCH" || command == "JUMP" {
			input.Scanf("%d", &operand)
			G[label].in = append(G[label].in, operand)
		}
		input.Scanf("\n")
	}
	for i, k := range G { for _, j := range k.in { G[j].out = append(G[j].out, i) } }
	DFS(start, nil)
	Dominators()
	for _,i := range tmp { 
                if 0 != G[i].time  { 
                        v := G[i]
                        for _, i := range v.out {
        	                dom := G[i]
		                for dom != nil && dom != v { dom = dom.dom }
		                if dom == v { cycle[v.indx], rs = true, rs + 1 }
	                }
                } 
        }
        fmt.Println(len(cycle))
}

