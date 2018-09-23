; ================= LAB 2 =====================
; ------------- 2. Trace ---------------
(define-syntax trace-ex
  (syntax-rules ()
    ((_ expr) (begin
                  (write (quote expr))
                  (display " => ")
                  (write expr)
                  (newline)
                  expr))))

(define (zip . xss)
  (if (or (null? xss)
          (null? (trace-ex (car xss)))) ; Здесь...
      '()
      (cons (map car xss)
            (apply zip (map cdr (trace-ex xss)))))) ; ... и здесь

(zip '(1 2 3) '(one two three))


; ---------------- 4. ref -------------------
(define (ref seq index)
  (cond
    ((vector? seq) (if (and (< index (vector-length seq))
                            (>= index 0))
                       (vector-ref seq index)
                       #f))
        ((list? seq) (if (and (< index (length seq))
                              (>= index 0))
                         (list-ref seq index)
                         #f))
        ((string? seq) (if (and (< index (string-length seq))
                                (>= index 0))
                           (string-ref seq index)
                           #f))
        (else #f)))

(define suite
  (list
   (test (ref '(1 2 3) 1) 2)
   (test (ref '(1 2 3) 3) #f)
   (test (ref '(1 2 3) -1) #f)
   
   (test (ref #(1 2 3) 1) 2)
   (test (ref #(1 2 3) 3) #f)
   (test (ref #(1 2 3) -1) #f)
   
   (test (ref "123" 1) #\2)
   (test (ref "123" 3) #f)
   (test (ref "123" -1) #f)))
(run-tests suite)



(define (ins seq index value)
  (define (list-ins lst idx val)
    (let* ((tail (list-tail lst index))
           (head (reverse (list-tail (reverse lst) (- (length lst) index)))))
      (append head (list value) tail))) ; "substring(0..k) + val + substring(k+1..len)
  (cond
    ((vector? seq) (if (and (<= index (vector-length seq))
                            (>= index 0))
                       (list->vector (list-ins (vector->list seq) index value))
                       #f))
    ((list? seq) (if (and (<= index (length seq))
                          (>= index 0))
                     (list-ins seq index value)
                     #f))
    ((string? seq) (if (and (<= index (string-length seq))
                            (>= index 0)
                            (char? value))
                       (list->string (list-ins (string->list seq) index value))
                       #f))
    (else #f)))

(define suite
  (list
   (test (ins '(1 2 3 4 5) 2 0) '(1 2 0 3 4 5))
   (test (ins '(1 2 3) 3 0) '(1 2 3 0))
   (test (ins '(1 2 3) 0 0) '(0 1 2 3))
   (test (ins '(1 2 3) -1 0) #f)
   (test (ins '(1 2 3) 5 0) #f)
   
   (test (ins #(1 2 3 4 5) 2 0) #(1 2 0 3 4 5))
   (test (ins #(1 2 3) 3 0) #(1 2 3 0))
   (test (ins #(1 2 3) 0 0) #(0 1 2 3))
   (test (ins #(1 2 3) 1 #\0) #(1 #\0 2 3))
   (test (ins #(1 2 3) -1 0) #f)
   (test (ins #(1 2 3) 5 0) #f)
   
   (test (ins "123" 1 #\0) "1023")
   (test (ins "123" 0 #\0) "0123")
   (test (ins "123" 3 #\4) "1234")
   (test (ins "123" 1 0) #f)
   (test (ins "123" -1 #\0) #f)
   (test (ins "123" 5 #\4) #f)

   (test (ins '123 1 0) #f)))
(run-tests suite)

