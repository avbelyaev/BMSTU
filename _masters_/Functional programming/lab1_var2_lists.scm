; ======= VAR 2 ========
; ------ 2. delete
(define (delete pred xs)
  (if (null? xs)
      '()
      (if (pred (car xs))
          (delete pred (cdr xs))
          (cons (car xs) (delete pred (cdr xs))))))

(display "delete(even?, [0, 1, 2, 3]):")
(delete even? '(0 1 2 3))

; ------ 3. iterate
(define (iterate f x n)
  (if (zero? n)
      '()
      (cons x (iterate f (f x) (- n 1)))))

(display "iterate(x -> 2*x, 1, 6):")
(iterate (lambda (x) (* 2 x)) 1 6)


; ------ 4. intersperse
(define (sperse e xs)
  (if (null? xs)
      '()
      (if (eq? e (car xs))
          (cons (car xs) (sperse e (cdr xs)))
          (cons (car xs) (sperse e (cons e (cdr xs)))))))

(display "intersperse('x' [1, 2, 3, 4]):")
(sperse 'x '(1 2 3))


; ------ 5. any-all
(define (all pred xs)
  (or (null? xs)
      (and (pred (car xs)) (all pred (cdr xs)))))

(define (any pred xs)
  (and (not (null? xs))
       (or (pred (car xs)) (any pred (cdr xs)))))


(define (even n)
  (zero? (remainder n 2)))

(define (odd n)
  (not (even n)))

(display "all? even? [1, 3, 5, 7]:")
(all even '(1 3 5 7))

(display "all? odd? [1, 3, 5, 7]:")
(all odd '(1 3 5 7))

(display "any? even? [1, 3, 5, 7]:")
(any even '(1 3 5 7))

(display "any? odd? [1, 3, 5, 7]:")
(any odd '(1 3 5 7))

