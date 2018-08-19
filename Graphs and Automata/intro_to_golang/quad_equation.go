package main
import "fmt"
import "math"
import "math/cmplx"

var s_flt string = "float"
var s_cmp string = "complex"

type mytype_f float64
type mytype_c complex128

type Number interface {
        Add(x Number) Number
	Sub(x Number) Number
	Mul(x Number) Number
	Div(x Number) Number
	MulInt(k int) Number
	Sqrt() Number
	String() string
}
//==============================================
func (a mytype_f) Add(x Number) Number {
	var ret, ok = x.(mytype_f)
        var bad_ret mytype_f = 0
	if ok {	return a + ret }
	return bad_ret
}
func (a mytype_c) Add(x Number) Number {
	var ret, ok = x.(mytype_c)
        var bad_ret mytype_c = 0
	if ok {	return a + ret }
	return bad_ret
}
//==============================================
func (a mytype_f) Sub(x Number) Number {
	var ret, ok = x.(mytype_f)
        var bad_ret mytype_f = 0
	if ok {	return a - ret }
	return bad_ret
}
func (a mytype_c) Sub(x Number) Number {
	var ret, ok = x.(mytype_c)
        var bad_ret mytype_c = 0
	if ok {	return a - ret }
	return bad_ret
}
//==============================================
func (a mytype_f) Mul(x Number) Number {
	var ret, ok = x.(mytype_f)
        var bad_ret mytype_f = 0
	if ok {	return a * ret }
	return bad_ret
}
func (a mytype_c) Mul(x Number) Number {
	var ret, ok = x.(mytype_c)
        var bad_ret mytype_c = 0
	if ok {	return a * ret }
	return bad_ret
}
//==============================================
func (a mytype_f) Div(x Number) Number {
	var ret, ok = x.(mytype_f)
        var bad_ret mytype_f = 0
	if ok {	return a / ret }
	return bad_ret
}
func (a mytype_c) Div(x Number) Number {
	var ret, ok = x.(mytype_c)
        var bad_ret mytype_c = 0
	if ok { return a / ret }
	return bad_ret
}
//==============================================
func (a mytype_f) MulInt(k int) Number {
	var ret = mytype_f(k)
	return a * ret
}
func (a mytype_c) MulInt(k int) Number {
	int_part := float64(k)
	var i_part float64 = 0
	comp_construct := complex(int_part, i_part)
	t_c_c := mytype_c(comp_construct)
	return a * t_c_c
}
//==============================================
func (a mytype_f) Sqrt() Number {
	var t_a = float64(a)
	ret := math.Sqrt(t_a)
	return mytype_f(ret)
}
func (a mytype_c) Sqrt() Number {
	var t_a = complex128(a)
	ret := cmplx.Sqrt(t_a)
	return mytype_c(ret)
}
//==============================================
func (a mytype_f) String() string { return fmt.Sprintln(a) }
func (a mytype_c) String() string { return fmt.Sprintln(a) }
//==============================================

func solve(a, b, c Number) (x1, x2 Number) {
        x1 = ((b.MulInt(-1)).Add(((b.Mul(b)).Sub((a.Mul(c)).MulInt(4))).Sqrt())).Div(a.MulInt(2))
	x2 = ((b.MulInt(-1)).Sub(((b.Mul(b)).Sub((a.Mul(c)).MulInt(4))).Sqrt())).Div(a.MulInt(2))
	return x1, x2
}

func type_change(a, b, c mytype_f, ptr int) {
        var x1, x2, in_a, in_b, in_c Number

	if 1 == ptr {
		in_a, in_b, in_c = a, b, c

		x1, x2 = solve(in_a, in_b, in_c)

		var temp_x2 mytype_f = x2.(mytype_f)
		var temp_x1 mytype_f = x1.(mytype_f)
		var zero_f mytype_f = 0.00
		if (temp_x1 == -temp_x2) && (temp_x1 == zero_f) { x2 = x1 }
		fmt.Printf("%.2f %.2f", x1, x2)
	}
	if 2 == ptr {
		var i_part mytype_f = 0

		t_a_c := mytype_c(complex(a, i_part))
		t_b_c := mytype_c(complex(b, i_part))
		t_c_c := mytype_c(complex(c, i_part))

		in_a, in_b, in_c = t_a_c, t_b_c, t_c_c

		x1, x2 = solve(in_a, in_b, in_c)

		var temp_x2 mytype_c = x2.(mytype_c)
		var ret_val_real float64 = float64(real(temp_x2))
		var ret_val_imag float64 = float64(imag(temp_x2))
		var zero_c float64 = 0.00

		if ret_val_imag == zero_c {
			ret_val_imag = ret_val_imag
			//fmt.Printf("imaginable new  =%.2f\n", ret_val_imag)
			new_x2 := mytype_c(complex(ret_val_real, zero_c))
			fmt.Printf("%.2f ", x1)
			fmt.Printf("%.2f", new_x2)
			return
		}

		fmt.Printf("%.2f %.2f", x1, x2)
	}
}

func main() {
	var s_type string
	input.Scanf("%s", &s_type)
	var tc_a, tc_b, tc_c float64
	input.Scanf("%f %f %f", &tc_a, &tc_b, &tc_c)
        var a, b, c mytype_f = mytype_f(tc_a), mytype_f(tc_b), mytype_f(tc_c)
	if s_type == s_flt {
		p := 1
		type_change(a, b, c, p)
	} else {
		p := 2
		type_change(a, b, c, p)
	}
}

