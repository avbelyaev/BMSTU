(define (integrate f a b n)
   (define x (/ (abs (- b a)) n))
       (define (sum y k) (if (< y b)n(sum (+ y x) (+ k (abs (* x (f (+ y x)))))) k))
(sum a 0.0))

