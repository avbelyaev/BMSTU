package main

import "fmt"
import "strconv"

func main() {
	signs := make([]string, 0)
	chars := make([]string, 0)
	brcks := make([]string, 0)
	dict := make(map[string]string)
	
//delete expression that has already been used
	lex_cut := func(ptr bool) {
		if ptr != true {
                        brcks = brcks[:len(brcks)-1] //'('
			signs = signs[:len(signs)-1]  //'$'
			chars = chars[:len(chars)-2]   //'ab'
			brcks = brcks[:len(brcks)-1]     //')'
			return
		} else {
			chars = chars[:len(chars)-2]
			return
		}
	}
	j := 1
	var brcls, c, s, res int
	//var str string
	newj := strconv.Itoa(j)
	
	fmt.Scanf ("%d\n", &str)
	if len(str) == 1 { goto exit }
		
	for _, x := range str {
			switch {
			case x == '#':
				signs = append(signs, string(x))
				s++
				//fmt.Printf("-#-")
			case x == '$':
				signs = append(signs, string(x))
				s++
				//fmt.Printf("-$-")
			case x == '@':
				signs = append(signs, string(x))
				s++
				//fmt.Printf("-@-")
			case 'a' <= x && x <= 'z':
				chars = append(chars, string(x))
				c++
				//fmt.Printf("-a-")
			case x == '(':
				continue
			case x == ')':
				//fmt.Printf(">> ")
				brcks = append(brcks, string(x))
				brcls++
				min_expr := RPN_to_STD_notation(signs, chars, brcls)
				//check if it exists in dictionary
				_, ok := dict[min_expr]
				
				if ok {
					//if exists -> delete
					brcks = brcks[:len(brcks)-1]
					signs = signs[:len(signs)-1]
					chars := chars[:len(chars)-2]
                                        brcks = brcks[:len(brcks)-1]
					chars[len(chars)-1] = string(newj)
				} else {
					//else mark and add
					dict[min_expr] = string(brcls)
					lex_cut(ok)
					chars = append(chars, string(brcls))
					res++
				}

			}
		}
	}
	fmt.Printf("%d", res)
        exit:
	fmt.Printf("0")
}

