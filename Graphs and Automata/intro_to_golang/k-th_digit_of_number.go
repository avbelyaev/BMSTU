package main
import "fmt"
import "math"

func main() {
        var a, i, k, temp int64
	fmt.Scanf("%d", &k)

	for ; k >= a; i++ {
		temp = a
		a+=  (int64)(9*math.Pow(10, float64(i))*float64(i+1))
	}

	c := i - 1
	d := (k - temp) / i
	temp2 := (int64)(math.Pow(10, float64(c)) + float64(d))
	k = i - (k-temp)%i
	
	for ;k >= 2; k-- { temp2 /= 10 }
	fmt.Printf("%d", temp2%10)
}

