(define (newton f df x e)
  (if (< (abs (/ (f x) (df x))) e)
      x
      (newton f df (- x (/ (f x) (df x))) e)))
 
(define (golden f x0 x1 e)
  (define g (/ (- (sqrt 5) 1) 2))
  (define (big-segment g a b) (* g (abs (- a b))))
  (if (> (f (- x1 (big-segment g x1 x0))) (f (+ x0 (big-segment g x1 x0))))
      (if (< (abs (- x1 x0)) e)
          (/ (+ x0 x1) 2)
          (golden f (- x1 (big-segment g x1 x0)) x1 e))
      (if (< (abs (- x1 x0)) e)
          (/ (+ x0 x1) 2)
          (golden f x0 (+ x0 (big-segment g x1 x0)) e))))

