val fib_helper: (Int, Int => Boolean, Int, Int) => List[Int] = {
case (n, p, p1, p2) if (n < p2) => List()
case (n, p, p1, p2) if(p(p2)) => p2::fib_helper(n, p, p2, p1 + p2)
case (n, p, p1, p2) => fib_helper(n, p, p2, p1 + p2)
}

val main: (Int, Int => Boolean) => List[Int] = {
(n,p) => fib_helper (n, p, 0, 1)
}







val fib: (Int) => List[Int] = {
case (0) => List(1)
case (1) => 1 :: fib(0)
case (x) => (fib(x-1)(0) + fib(x-2)(0)) :: fib(x-1)
}

val sort_top: (List[Int], Int) => List[Int] = {
case (Nil, a) => List()
case (x :: xs, a) if(x < a) => x :: sort_top(xs, a)
case (x :: xs, a) => sort_top(xs, a)
}

val sort_pred: (List[Int], Int => Boolean) => List[Int] = {
case (Nil, a) => List()
case (x :: xs, p) if(p(x)) => x :: sort_pred(xs, p)
case (x :: xs, p) => sort_pred(xs, p)
}

val main: (Int, Int => Boolean) => List[Int] = {
(n, p) => sort_pred(sort_top(fib(n), n), p)
}




val fib: Int => Int = {
case x if(x<0) => 0
case 0 | 1 => 1
case x => fib(x-1)+fib(x-2)
}

val add: Int => List[Int] = {
case (0) => List(fib(0))
case (x) => x :: List

val chk: (Int, Int => Boolean) => List[Int] = {
case (x, p) if (p(x)) => (fib(x)(0)) :: chk(x-1, p)
case (x, p) => chk(x, p)
}

case x => (fib(x-1)(0) + fib(x-2)(0)) :: fib(x-1)
case x => ((x :: fib(x-1)) + (x :: fib(x-2))) :: fib(x-1)
8, 5, 3, 2, 1, 1

case x			=> (x + (x :: fib(x-1))) :: List(fib(x-1))