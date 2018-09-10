; ========== DZ-1
; ------- map
(define (m-map proc xs)
  (if (null? xs)
      '()
      (cons (proc(car(xs))
                 (m-map proc(cdr xs))))))


; ------ list-ref
(define (m-list-ref xs k)
  (if (zero? k)
      (car xs)
      (list-ref (cdr xs) (- k 1))))

(display "list-ref([1 4 8 8], 1):")
(m-list-ref '(1 4 8 8) 1)


; ------ append
(define (m-append xs x)
  (if (null? xs)
      (cons x '())
      (cons (car xs)
            (m-append (cdr xs) x))))

(display "append([1 3 3], 7):")
(m-append '(1 3 3) 7 )


; ----- filter
(define (m-filter pred xs)
  (if (null? xs)
      '()
      (if (pred (car xs))
          (cons (car xs) (m-filter pred (cdr xs)))
          (m-filter pred (cdr xs)))))
      
(display "filter(x -> 0 == x, [1 0 3])")
(m-filter zero? '(1 0 3))

(display "filter(x -> 3 <= x, [2 3 4])")
(m-filter (lambda (x) (>= x 3)) '(2 3 4))


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
(define (compose x . procs)
  (if (null? procs)
      1
      ((car procs) car (apply compose (cdr procs)))))

(define (f x) (* x 2))
(define (g x) (* x 3))
(display "compose(f(g(1))):")
;(compose 1 f g) not working atm


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

