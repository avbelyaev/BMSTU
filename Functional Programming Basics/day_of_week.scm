(define (day-of-week x y z)
  (let* (
         (a (quotient (- 14 y) 12))
         (m (- z a))
         (k (- (+ y (* 12 a)) 2)))
         (remainder (+ 7000 (+(- (+ x m (quotient m 4)) (quotient m 100)) (quotient m 400) (quotient (* 31 k) 12))) 7
        )
   )
)

