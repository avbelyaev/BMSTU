package main
import "fmt"
import "github.com/skorobogatov/input"

func main() {
	var x, y rune
	var i, min, current, pos_x, pos_y int
        min = 1000000
        
	str := input.Gets ()
	input.Scanf("%c %c", &x, &y)
 
        var first, last int
        
        i = 0
	for _, p := range str {
                i++
		if p == x { //нашли x
			pos_x = i 
		}
		if p == y { //нашли y
			pos_y = i 
		}
                
		if (p == x || p == y) && (0 < pos_x && 0 < pos_y) { 
                        current = pos_x - pos_y
                        /*if 0 == current { 
                                min = current
			        goto exit
	                }*/
			if 0 > current {
				current *= -1 
			}
			if current < min {
				min = current
			}
                        if 0 == (min - 1) {
                                goto exit
                        }
		}
	}
        exit:
        min -= 1
	fmt.Printf("%d", min)
}

