(define (my-gcd a b)
     (if (= a b)
          a
         (if (< a b)
             (my-gcd a (- b a))
             (my-gcd (- a b) b))))
 
 
(define (my-lcm a b)
  (/ (* a b) (my-gcd a b)))
 
(define (prime? n)
  (not (or (= n 1)
           (= (remainder n 2) 0)
           (= (remainder n (do ((i 2 (+ i 1)))
                               ((> i (/ n 2)) i))) 0))))

