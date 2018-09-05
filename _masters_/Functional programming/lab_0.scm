;#lang racket

(define sq
  (lambda (x) (* x x)))

(define (my-abs x) (if (< x 0) (-x) (x)))

(define (sign x)
  (cond
    ((< x 0) -1)
    ((= x 0) 0)
    (else 1)))


(define (fact n)
  (if (= n 0)
      1
      (* n (fact (- n 1)))))


; =========

(define (m-even n)
  (zero? (remainder n 2)))

(display 'even\(2\):)
(m-even 2)
(display 'even\(5\):)
(m-even 5)


; ------
(define (m-odd n)
  (not (m-even n)))

(display 'odd\(2\):)
(m-odd 2)
(display 'odd\(5\):)
(m-odd 5)



; ------
(define (m-gcd a b)
  (if (= b 0)
      a
      (m-gcd b (remainder a b))))

(display 'gcd\(5_2\):)
(m-gcd 5 2)
(display 'gcd\(10_15\):)
(m-gcd 10 15)



; ------
(define (m-lcm a b)
  (quotient (* a b) (m-gcd a b)))

(display 'lcm\(12_18\):)
(m-lcm 12 18)



; -------
(define (prm-helper x k)
  (if (= k 1)
      #t
      (if (zero? (remainder x k))
             #f
             (prm-helper x (- k 1)))))

(define (prm n)
  (prm-helper n 2))

(display 'prime\(2\):)
(prm 2)
(display 'prime\(7\):)
(prm 7)
(display 'prime\(40\):)
(prm 40)


; ==========
(define (f x y)
  (let ((a (* x x))
        (b (* y y)))
    (+ a b)))


; ==========
(define (m-list-ref xs k)
  (if (zero? k)
      (car xs)
      (list-ref (cdr xs) (- k 1))))

   
;========

(define (my-map proc xs)
  (if (null? xs)
      '()
      (cons (proc(car(xs))
                 (my-map proc(cdr xs))))))

