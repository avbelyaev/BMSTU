package main
import "fmt"
import "os"
 
type frc struct {
        numer int
        denom int
        //           numerator
        //fraction=  ----------
        //           denominator
}
 
//print matrix
func output_matrix(matrix [][]frc, n int) {
        for i := 0; i < n; i++ {
                for j := 0; j < n+1; j++ {
                        fmt.Printf("%d ", matrix[i][j].numer)
                }
                fmt.Printf("\n")
        }
        fmt.Printf("\n")
        return
}
 
//scan matrix
func get_matrix(matrix [][]frc, n int) {
        fmt.Printf("Insert elements:\n")
        for i := 0; i < n; i++ {    //i-th row (row), j-th symbol in a row (column)
                for j := 0; j < n+1; j++ {
                        fmt.Scanf("%d", &(matrix[i][j].numer))
                        matrix[i][j].denom = 1
                }
                fmt.Scanf("\n")
        }
        return
}
 
//swap 2 rows
func swap_rows(matrix [][]frc, str_1, str_2, n int) {
        for j := 0; j < n+1; j++ {
                matrix[str_1][j].numer, matrix[str_2][j].numer = matrix[str_2][j].numer, matrix[str_1][j].numer
                matrix[str_1][j].denom, matrix[str_2][j].denom = matrix[str_2][j].denom, matrix[str_1][j].denom
        }
        return
}
//================================================================================
func gcd(x, y int) int {
        for y != 0 { x, y = y, x%y }
        return x
}
 
func lcm(m, n int) int { return m / gcd(m, n) * n }
 
func fraction_reduce(x frc) frc { //   14/88 -> 7/44
        for 1 != gcd(x.numer, x.denom) {
                i := gcd(x.numer, x.denom)
                x.numer /= i
                x.denom /= i
        }
        return x
}
//=================================================================================

//check for negative diag elements (if exist -> mulitiply them by '-1')
func neg_diag_elements(matrix [][]frc, n, j int) {
        for i := j; i < n; i++ {
                if matrix[i][i].numer < 0 { for k := 0; k < n+1; k++ { matrix[i][k].numer *= -1 } } //multiplied by '-1' 
        }
}
 
//check if null diag elements exist
func check_null_diag(matrix [][]frc, n, i int) [][]frc {
        for i < n {
                if 0 == matrix[i][i].numer {
                        for t := i + 1; t < n; t++ {
                                if 0 != matrix[t][i].numer {
                                        swap_rows(matrix, i, t, n)
                                        break
                                }
                        }
                }
                i++
        }
        return matrix
}
 
func direct_gauss_method(matrix [][]frc, n int) {
        for j := 0; j < n; j++ {
                check_null_diag(matrix, n, j)
                //fmt.Printf("_check_null_diag_direct\n")
                output_matrix(matrix, n)
                for i := j + 1; i < n; i++ {
                       
                        neg_diag_elements(matrix, n, j)
 
                        a, b := matrix[j][j].numer, matrix[i][j].numer
                        //fmt.Printf("a=%d, b=%d, ", a, b)
                        if 0 == b { continue }
                        if 0 == a {
                                fmt.Printf("No solution :(\n")
                                os.Exit(0)
                        }

                        g := lcm(a, b)
                        k_a, k_b := g/a, g/b
                        //fmt.Printf("lcm=%d, k_a=%d, k_b=%d\n", g, k_a, k_b)
                        
			for k := 0; k < n+1; k++ { matrix[i][k].numer *= k_b } //mul
                        for k := 0; k < n+1; k++ { matrix[i][k].numer -= k_a * (matrix[j][k].numer) } //sub
 
                        output_matrix(matrix, n)
                }
        }
}
 
func reverse_gauss_method(matrix [][]frc, n int) {
        for j := n - 1; j > 0; j-- {
                output_matrix(matrix, n)
                for i := j - 1; i >= 0; i-- {
                        neg_diag_elements(matrix, n, j)
 
                        a, b := matrix[j][j].numer, matrix[i][j].numer
                        //fmt.Printf("a=%d, b=%d, ", a, b)
                        if 0 == b { continue }
                        if 0 == a {
                                fmt.Printf("No solution :(\n")
                                os.Exit(0)
                        }
 
                        g := lcm(a, b)
                        k_a, k_b := g/a, g/b
                        //fmt.Printf("lcm=%d, k_a=%d, k_b=%d\n", g, k_a, k_b)

                        for k := 0; k < n+1; k++ { matrix[i][k].numer *= k_b }
                        for k := 0; k < n+1; k++ { matrix[i][k].numer -= k_a * (matrix[j][k].numer) }
 
                        output_matrix(matrix, n)
                }
        }
}
 
func get_answers(matrix [][]frc, n int) {
        for i := 0; i < n; i++ {
                if 0 == matrix[i][i].numer {
                        fmt.Printf("No solution :(\n")
                        return
                }
        }

        for i := 0; i < n; i++ {
                var ans frc
                ans.numer, ans.denom = matrix[i][n].numer, matrix[i][i].numer
                x := fraction_reduce(ans)
                	if x.denom < 0 && x.numer > 0 { x.numer, x.denom = -x.numer, -x.denom }
                fmt.Printf("X%d=%d/%d\n", i, x.numer, x.denom)
        }
}
 
func gauss(matrix [][]frc, n int) int {
   
        direct_gauss_method(matrix, n)
        //fmt.Printf("_direct_gauss: DONE\n")
        reverse_gauss_method(matrix, n)
        //fmt.Printf("_reverse_gauss: DONE\n\n")
        get_answers(matrix, n)
 
        /* calculate one more time... */
        var mark int
        fmt.Printf("\n\nCalculate another matrix? ('1' == YES; any == NOPE):")
        fmt.Scanf("%d\n", &mark)

        if mark == 1 { return 1 } else { return 0 }
}
 
func main() {
    one_more_time:
        var n, i int
        fmt.Printf("Insert N - rank of matrix (N rows; N+1 columns):\n")
        fmt.Scanf("%d\n", &n)
 
        matrix := make([][]frc, n)
        for i = 0; i < n; i++ { matrix[i] = make([]frc, n+1) }
 
        get_matrix(matrix, n)
        mark := gauss(matrix, n)
 
        /*calculate one more time */
        if mark == 1 { goto one_more_time }
}

