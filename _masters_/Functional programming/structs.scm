; ====== Structs ========
; ------- Table ---------
(define (replicate x n)
  (if (zero? n)
      '()
      (cons x (replicate x (- n 1)))))

(display "replicate([a b], 3):")
(replicate '(a b) 3)

; data is stored after this offset
(define offset 2)

(define (table rows cols . default-value)
  (list->vector (append
                 (list rows cols)
                 (replicate (if (null? default-value)
                                #f
                                (car default-value))
                            (* rows cols)))))

(display "table(rows=3, cols=4, default=0):")
(define t (table 3 4 0))
(display t)
(newline)

(define (table? t)
  (let* ((ts (vector->list t))
        (rows (car ts))
        (cols (cadr ts))
        (cells-number (+ offset (* rows cols))))
    (and (vector? t)
         (> rows 0)
         (> cols 0)
         (= cells-number (length ts)))))

(display "table? t:")
(table? t)


(define (cell-index t i j)
  (let* ((ts (vector->list t))
         (cols (cadr ts)))
    (+ (* cols i) j offset)))


(define (table-set! t i j value)
  (if (table? t)
      (vector-set! t (cell-index t i j) value)))

(display "table[2][3]=1337;") (newline)
(table-set! t 2 3 1337)
(display t)
(newline)


(define (table-ref t i j)
  (if (table? t)
      (vector-ref t (cell-index t i j))))

(display "table[2][3]:")
(table-ref t 2 3)


; ------ symbol-append ------
(define (symbol-append . symbols)
  (define (helper result-str symbols)
    (if (null? symbols)
        result-str
        (helper (string-append result-str (symbol->string (car symbols)) "-") (cdr symbols))))
  (let* ((concatenated-str (helper "" symbols))
         (removed-last-dash (substring concatenated-str 0 (- (string-length concatenated-str) 1))))
    (string->symbol removed-last-dash)))
  
(define s (symbol-append 'very 'long 'symbol))
(display s)
(symbol? s)
