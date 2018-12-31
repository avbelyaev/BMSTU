package main

import "fmt"

func stark(x int) int {
	return x * 1337
}

func mark(a int) int {
	return a + 1
}

func main() {

	const abc = "DEF"

	fmt.Printf("hello, world\n")

	var j int

	j = 5

	if j < 3 {
		j++
	}

	//greetings

	j = mark(j)

	fmt.Printf("J=%d\n", j)
}
