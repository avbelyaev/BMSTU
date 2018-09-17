; ========== VAR-1
; ------ 1. count
(define (count-helper x xs res)
  (if (null? xs)
      res
      (count-helper x (cdr xs) (if (equal? x (car xs))
                                   (+ res 1)
                                   res))))
   
(define (count x xs)
  (count-helper x xs 0))

(display "count('a, [a b c a]):")
(count 'a '(a b c a))


; ------ 2. replace
(define (replace pred proc xs)
  (if (null? xs)
      '()
      (cons (if (pred (car xs))
                (proc (car xs))
                (car xs))
            (replace pred proc (cdr xs)))))

(display "replace(even?, x -> x * 10, [0 1 2 3 4]):")
(replace even? (lambda (x) (* x 10)) '(0 1 2 3 4))


; ------- 3. replicate
(define (replicate x n)
  (if (zero? n)
      '()
      (cons x (replicate x (- n 1)))))

(display "replicate([a b], 3):")
(replicate '(a b) 3)


; ------- 4. cycle
(define (cycle-helper xs n ys)
  (if (zero? n)
      '()
      (if (null? ys)
          (cycle-helper xs (- n 1) xs)
          (cons (car ys) (cycle-helper xs n (cdr ys))))))

(define (cycle xs  n)
  (cycle-helper xs n xs))

(display "cycle([a b c], 3):")
(cycle '(a b c) 3)


; ------ 5. and/or-fold
(define (and-fold . args)
  (if (null? args)
      #t
      (if (car args)
          (apply and-fold (cdr args))
          #f)))

(display "and-fold(#t #t #f):")
(and-fold #t #t #f)
(display "and-fold(#t #t #t):")
(and-fold #t #t #t)


; ----- 6. composition
(define (compose . procs)
  (if (null? procs)
      (lambda (x) x)
      (lambda (x) ((car procs) ((apply compose (cdr procs)) x)))))

(define (f x) (* x 2))
(define (g x) (* x 3))
(define (h x) (- x))

(display "compose(f(g(h(1)))):")
((compose f g h) 1)

(display "compose(f(g(1))):")
((compose f g) 1)


; ----- 7. find-number
(define (find-number from to divisor)
  (if (> from to)
      #f
      (if (zero? (remainder from divisor))
          from
          (find-number (+ from 1) to divisor))))

(display "find-number(7, 9, 3):")
(find-number 7 9 3)

(display "find-number(3, 7, 9):")
(find-number 3 7 9)
