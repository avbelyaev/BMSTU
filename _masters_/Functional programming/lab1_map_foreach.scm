; ----- foreach -----
(define (each proc xs)
  (if (not (null? xs))
      (begin
        (proc (car xs))
        (each proc (cdr xs)))))

(display "for-each(x -> print 'x...', [1, 2, 3, 4]):")
(each (lambda (x) (display x) (display "...")) '(1 2 3 4))
(newline)


; ------- map
(define (m-map proc xs)
  (if (null? xs)
      '()
      (cons (proc (car xs)) (m-map proc (cdr xs)))))

(display "map(x -> x * 10, [1, 3, 3, 7]):")
(m-map (lambda (x) (* x 10)) '(1 3 3 7))


; ------- lib-map
(define (get-nth n xs)
  (if (zero? n)
      (car xs)
      (get-nth (- n 1) (cdr xs))))

(display "get-nth(1, [1, 2, 3, 4]):")
(get-nth 1 '(1 2 3))


(define (get-nths n lists)
  (if (null? lists)
      '()
      (cons (get-nth n (car lists)) (get-nths n (cdr lists)))))

(display "get-nths(1, [1, 2, 3], [4, 5, 6], [7, 8, 9]):")
(get-nths 1 '((1 2 3) (4 5 6) (7 8 9)))


(define (lib-map proc . lists)
  (define (mapper proc n lists)
    (if (null? lists)
        '()
        (if (= n (length (car lists)))
            '()
            (cons (apply proc (get-nths n lists)) (mapper proc (+ n 1) lists)))))
  (mapper proc 0 lists))


(display "lib-map(+, [1, 2, 3], [10, 20, 30], [100, 200, 300]):")
(lib-map + '(1 2 3) '(10 20 30) '(100 200 300))

