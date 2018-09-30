; ---------- utils -----------
(define (and-fold . args)
  (if (null? args)
      #t
      (if (car args)
          (apply and-fold (cdr args))
          #f)))

(define ie (interaction-environment))

; ------- unit tests -----------
(define-syntax test
  (syntax-rules ()
    ((_ expression expected-result)
     (list (quote expression)
           expected-result))))


(define (run-test the-test)
  (let ((expression (car the-test))
        (expected-result (cadr the-test)))
    ; show before crash
    (write expression)
    (display " => ")
    (write expected-result)
    (let* ((actual-result (eval expression (interaction-environment)))
            (status (equal? expected-result actual-result)))
      (if status
          (display " ok")
          (display " FAIL"))
      (newline)
      (if (not status)
          (begin
            (display "   expected: ")
            (write (cadr the-test))
            (newline)
            (display "     actual: ")
            (write actual-result)
            (newline)))
      status)))


(define (run-tests the-tests)
  (and-fold #t (map run-test the-tests)))



; =================== Deriv =======================
(display "================= Deriv ===================")
(newline)



(define (eq-number? x num)
  (and (number? x)
       (= x num)))


(define (var? v)
  (or (symbol? v)
      (number? v)))

(define suite
  (list
   (test (var? 'a) #t)
   (test (var? 99) #t)
   (test (var? '(a)) #f)
   (test (var? '(1 2 3)) #f)))
(run-tests suite)



(define (same-var? v1 v2)
  (and (var? v1)
       (var? v2)
       (equal? v1 v2)))

(define suite
  (list
   (test (same-var? 'a 'a) #t)
   (test (same-var? 'a 'a) #t)
   (test (same-var? 11 11) #t)
   (test (same-var? 'a 'b) #f)
   (test (same-var? 11 22) #f)))
(run-tests suite)


; ---- Head -----
; first elem is ignored
; (+ (a b) c d) => (a b)
(define (head s)
  (cadr s))

(define suite
  (list
   (test (head '(+ 1 2)) 1)
   (test (head '(/ 2 3)) 2)
   (test (head '(+ 3 2 1)) 3)
   (test (head '(+ (+ a b) c)) '(+ a b))
   (test (head '(+ a b c)) 'a)
   (test (head '(/ a b c)) 'a)))
(run-tests suite)


; ---- Tail ----
; (1 2 (3 4) 5) => (2 (3 4) 5)
(define (tail s)
  (define (helper xs)
    (if (null? xs)
        '()
        (cond (car xs) (helper (cdr xs)))))
  (helper (cddr s)))

(define suite
  (list
   (test (tail '(+ 1 2)) '(2))
   (test (tail '(/ 1 2)) '(2))
   (test (tail '(+ 3 2 1)) '(2 1))
   (test (tail '(+ a b c)) '(b c))))
(run-tests suite)
  


; (+ a b) => true //as binary operation
(define (bin-op? operator exp)
  (and (list? exp)
       (= 3 (length exp))
       (eq? operator (car exp))
       (or (var? (cadr exp)) (list? (cadr exp)))
       (or (var? (caddr exp)) (list? (caddr exp)))))

(define suite
  (list
   (test (bin-op? '+ '(+ 1 2)) #t)
   (test (bin-op? '+ '(+ a b)) #t)
   (test (bin-op? '- '(- 1 2)) #t)
   (test (bin-op? '- '(- a b)) #t)
   (test (bin-op? '+ '(+ 1 2 3)) #f)
   (test (bin-op? '+ '(+ 1 2 a)) #f)
   (test (bin-op? '- '(- 1 2 3)) #f)
   ))
(run-tests suite)


; ---- (+) tail ----
; (+ a b c d) => (+ b c d)
(define (addtail vars)
  (cond ((bin-op? '+ vars) (caddr vars))
        (else (cons '+ (tail vars)))))

(define suite
  (list
   (test (addtail '(+ 1 2)) '2)
   (test (addtail '(+ a b)) 'b)
   (test (addtail '(+ (+ a b) c)) 'c)
   (test (addtail '(+ 1 2 3)) '(+ 2 3))
   (test (addtail '(+ 1 2 a)) '(+ 2 a))
   (test (addtail '(+ a b c)) '(+ b c))
   (test (addtail '(+ a b c d)) '(+ b c d))))
(run-tests suite)


; ---- (-) tail ----
; (- a b c) => (- 0 b c) //it only looks incorrect :)
(define (subtail vars)
  (cond ((bin-op? '- vars) (caddr vars))
        (else (cons '- (cons 0 (tail vars))))))
  
(define suite
  (list
   (test (subtail '(- 1 2)) '2)
   (test (subtail '(- 1 2 3)) '(- 0 2 3))
   (test (subtail '(- a b)) 'b)
   (test (subtail '(- a b c)) '(- 0 b c))
   (test (subtail '(- a b c d)) '(- 0 b c d))))
(run-tests suite)


; ---- (*) tail ----
; (* a b c) => (* b c)
(define (multail vars)
  (cond ((bin-op? '* vars) (caddr vars))
        (else (cons '* (tail vars)))))
  
(define suite
  (list
   (test (multail '(* 1 2)) '2)
   (test (multail '(* 1 2 3)) '(* 2 3))
   (test (multail '(* a b)) 'b)
   (test (multail '(* a b c)) '(* b c))
   (test (multail '(* a b c d)) '(* b c d))))
(run-tests suite)


; ---- (/) tail ----
; (/ a b c) => (/ b c) //also looks incorrect :)
(define (divtail vars)
  (cond ((bin-op? '/ vars) (caddr vars))
        (else (cons '/ (tail vars)))))



; ============= Folding rules ==============

; (x y) => (+ x y)
(define (make-add a1 a2) 
  (cond ((eq-number? a1 0) a2)
        ((eq-number? a2 0) a1)
        ((and (number? a1) (number? a2)) (+ a1 a2))
        (else `(+ ,a1 ,a2))))

(define suite
  (list
   (test (make-add 'x 0) 'x)
   (test (make-add 0 'x) 'x)
   (test (make-add 1 2) 3)
   (test (make-add 'x 'x) '(+ x x))
   (test (make-add 'x 'y) '(+ x y))
   (test (make-add '(+ x y) 'z) '(+ (+ x y) z))
   ))
(run-tests suite)


; (x y) => (- x y)
(define (make-sub a1 a2) 
    (cond ((and (number? a1) (number? a2)) (- a1 a2))
          ((eq-number? a1 0) `(- 0 ,a2))
          ((eq-number? a2 0) a1)
          ((equal? a1 a2) 0)
          (else `(- ,a1 ,a2))))

(define suite
  (list
   (test (make-sub 'x 0) 'x)
   (test (make-sub 0 'x) '(- 0 x)) ; TODO simplify
   (test (make-sub 0 1) -1)
   (test (make-sub 1 1) 0)
   (test (make-sub 'x 'x) 0)
   (test (make-sub 'x 'y) '(- x y))
   (test (make-sub '(+ x y) 'y) '(- (+ x y) y))
   ))
(run-tests suite)


; (x y) => (* x y)
(define (make-mul m1 m2)
  (cond ((or (eq-number? m1 0) (eq-number? m2 0)) 0)
        ((eq-number? m1 1) m2)
        ((eq-number? m2 1) m1)
        ((and (number? m1) (number? m2)) (* m1 m2))
        ((equal? m1 m2) `(expt ,m1 2))
        (else `(* ,m1 ,m2))))

(define suite
  (list
   (test (make-mul 'x 0) 0)
   (test (make-mul 0 'x) 0)
   (test (make-mul 0 0) 0)
   (test (make-mul 'x 1) 'x)
   (test (make-mul 1 'x) 'x)
   (test (make-mul 2 3) 6)
   (test (make-mul 'x 'y) '(* x y))
   (test (make-mul 'x 'x) '(expt x 2))
   (test (make-mul 'x '(* y z)) '(* x (* y z)))
   ))
(run-tests suite)


; (x y) => (/ x y)
(define (make-div m1 m2)
  (cond ((eq-number? m1 0) 0)
        ((eq-number? m2 0) 'zero-div) ; zero-division
        ((eq-number? m2 1) m1)
        ((and (number? m1) (number? m2)) `(/ ,m1 ,m2)) ; count not needed
        (else `(/ ,m1 ,m2))))

(define suite
  (list
   (test (make-div 1 1) 1)
   (test (make-div 1 0) 'zero-div)
   (test (make-div 14 88) '(/ 14 88))
   (test (make-div 'x 1) 'x)
   (test (make-div 0 'x) 0)
   (test (make-div 'x 'y) '(/ x y))
   (test (make-div '(+ x y) 'y) '(/ (+ x y) y))
   (test (make-div '(* x y) '(+ x 3)) '(/ (* x y) (+ x 3)))
   ))
(run-tests suite)



(define (deriv f var)
  (cond ((number? f) 0)
        ((var? f) (if (same-var? f var) 1 0))
        ((list? f) (let ((op (car f)))
                     (case op
                       ('+ (make-add (deriv (head f) var)
                                     (deriv (addtail f) var)))
                       ('- (make-sub (deriv (head f) var)
                                     (deriv (subtail f) var)))
                       ('* (make-add ; uv' + u'v
                            (make-mul (head f)
                                      (deriv (multail f) var))
                            (make-mul (deriv (head f) var)
                                      (multail f))))
                       ('/ (make-div
                            (make-sub ; u'v - uv'
                             (make-mul (deriv (head f) var)
                                       (divtail f))
                             (make-mul (head f)
                                       (deriv (divtail f) var)))
                            (make-mul (divtail f) (divtail f))))
                       (else 'unknown-operator))))
        (else 'unknown-type)))

(define suite
  (list
   (test (deriv '1 'x) '0)
   (test (deriv 'x 'x) '1)
   (test (deriv 'y 'x) '0)
   (test (deriv '(+ x x) 'x) '2)
   (test (deriv '(+ x y) 'x) '1)
   (test (deriv '(+ x y) 'y) '1)
   (test (deriv '(+ x y) 'z) '0)
   
   (test (deriv '(- x x) 'x) '0)
   (test (deriv '(- x y) 'x) '1)
   (test (deriv '(- x y) 'y) '-1)
   (test (deriv '(- x y) 'z) '0)

   (test (deriv '(* x x) 'x) '(+ x x))
   (test (deriv '(* x y) 'x) 'y)
   (test (deriv '(* x 3) 'x) '3)
   (test (deriv '(* 2 3) 'x) '0)

   (test (deriv '(/ 1 1) 'x) '0)
   (test (deriv '(/ x 1) 'x) '1)
   (test (deriv '(/ x y) 'x) '(/ y (expt y 2))) ; TODO simplify

   (test (deriv '(+ x x x) 'x) '3)
   (test (deriv '(+ x y z) 'z) '1)
   (test (deriv '(+ x z y) 'z) '1)
   (test (deriv '(+ 1 2 z) 'z) '1)

   (test (deriv '(+ (+ x y) z) 'z) '1)
   (test (deriv '(+ (+ x y) (+ x y)) 'x) '2)
   (test (deriv '(+ (- x y) (+ x y)) 'x) '2)
   (test (deriv '(+ (- x y) (+ x y)) 'y) '0)
   (test (deriv '(+ (- x y) (- x y)) 'x) '2)
   (test (deriv '(+ (- x y) (- x y)) 'y) '-2)
   (test (deriv '(- (- x y) (- x y)) 'y) '0)
   (test (deriv '(+ (+ x x x) (+ x x y)) 'x) '5)

   (test (deriv '(* (* x y) (+ x 3)) 'x)
         '(+ (* x y) (* y (+ x 3))))

   ; d/dx: xy / (x+3)
   (test (deriv '(/ (* x y) (+ x 3)) 'x)
         '(/ (- (* y (+ x 3)) (* x y)) (expt (+ x 3) 2)))

   ; d/dx: (x+y) / 5xy
   (test (deriv '(/ (+ x y) (* 5 x y)) 'x)
         '(/ (- (* 5 x y) (* (+ x y) (* 5 y))) (expt (* 5 x y) 2)))
   ))
(run-tests suite)

