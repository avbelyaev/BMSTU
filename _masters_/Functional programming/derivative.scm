(define (and-fold . args)
  (if (null? args)
      #t
      (if (car args)
          (apply and-fold (cdr args))
          #f)))

(define-syntax trace-ex
  (syntax-rules ()
    ((_ expr) (begin
                  (write (quote expr))
                  (display " => ")
                  (write expr)
                  (newline)
                  expr))))

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
      ;(and (list? v)
       ;    (= 1 (length v)))))

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
  


(define (bin-op? op-given exp)
  ;(newline)
  ;(display "bin-op: ")
  ;(write exp)
  ;(newline)
  (and (list? exp)
       (= 3 (length exp))
       (eq? op-given (car exp))
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


  
(define (addtail s)
  ;(newline)
  ;(display "adtail: ")
  ;(write s)
  ;(newline)
  (cond ((bin-op? '+ s) (caddr s))
        (else (cons '+ (tail s)))))

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



(define (subtail s)
  (cond ((bin-op? '- s) (caddr s))
        (else (cons '- (cons 0 (tail s))))))
  
(define suite
  (list
   (test (subtail '(- 1 2)) '2)
   (test (subtail '(- 1 2 3)) '(- 0 2 3))
   (test (subtail '(- a b)) 'b)
   (test (subtail '(- a b c)) '(- 0 b c))
   (test (subtail '(- a b c d)) '(- 0 b c d))))
(run-tests suite)



(define (multail s)
  (cond ((bin-op? '* s) (caddr s))
        (else (cons '* (tail s)))))
  
(define suite
  (list
   (test (multail '(* 1 2)) '2)
   (test (multail '(* 1 2 3)) '(* 2 3))
   (test (multail '(* a b)) 'b)
   (test (multail '(* a b c)) '(* b c))
   (test (multail '(* a b c d)) '(* b c d))))
(run-tests suite)



(define (divtail s)
  (cond ((bin-op? '/ s) (caddr s))
        (else (cons '/ (tail s)))))



(define (make-add a1 a2) ; (x y) => (+ x y)
  ;(newline)
  ;(display "make-add a1:")
  ;(write a1)
  ;(display " a2:")
  ;(write a2)
  ;(newline)
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



(define (make-sub a1 a2) ; (x y) => (- x y)
    (cond ((eq-number? a1 0) `(- 0 ,a2))
          ((eq-number? a2 0) a1)
          ((and (number? a1) (number? a2)) (- a1 a2))
          (else `(- ,a1 ,a2))))

(define suite
  (list
   (test (make-sub 'x 0) 'x)
   (test (make-sub 0 'x) '(- 0 x))
   (test (make-sub 1 2) -1)
   (test (make-sub 1 1) 0)
   (test (make-sub 'x 'x) '(- x x))
   (test (make-sub 'x 'y) '(- x y))
   (test (make-sub '(+ x y) 'y) '(- (+ x y) y))
   ))
(run-tests suite)



(define (make-mul m1 m2) ; (x y) => (* x y)
    (cond ((or (eq-number? m1 0)
               (eq-number? m2 0)) 0)
          ((eq-number? m1 1) m2)
          ((eq-number? m2 1) m1)
          ((and (number? m1) (number? m2)) (* m1 m2))
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
   (test (make-mul 'x '(* y z)) '(* x (* y z)))
   ))
(run-tests suite)



(define (make-div m1 m2) ; (x y) => (/ x y)
    (cond ((eq-number? m1 0) 0)
          ((eq-number? m2 1) m1)
          ((and (number? m1) (number? m2)) `(/ ,m1 ,m2)) ; do not count
          (else `(/ ,m1 ,m2))))

(define suite
  (list
   (test (make-div 1 1) 1)
   (test (make-div 14 88) '(/ 14 88))
   (test (make-div 'x 1) 'x)
   (test (make-div 0 'x) 0)
   (test (make-div 'x 'y) '(/ x y))
   (test (make-div '(+ x y) 'y) '(/ (+ x y) y))
   (test (make-div '(* x y) '(+ x 3)) '(/ (* x y) (+ x 3)))
   ))
(run-tests suite)



(define (deriv f var)
  ;(newline)
  ;(display "deriv ")
  ;(write f)
  ;(display " ")
  ;(write var)
  ;(newline)
  (cond ((number? f) 0)
        ((var? f) (if (same-var? f var) 1 0))
        ((list f) (let ((op (car f)))
                     (case op
                       ('+ (make-add (deriv (head f) var)
                                     (deriv (addtail f) var)))
                       ('- (make-sub (deriv (head f) var)
                                     (deriv (subtail f) var)))
                       ('* (make-add
                            (make-mul (head f)
                                      (deriv (multail f) var))
                            (make-mul (deriv (head f) var)
                                      (multail f))))
                       ('/ (make-div
                            (make-sub
                             (make-mul (deriv (head f) var)
                                       (divtail f))
                             (make-mul (head f)
                                       (deriv (divtail f) var)))
                            (make-mul (divtail f) (divtail f)))))))))

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
   (test (deriv '(- x y) 'y) '(- 0 1)) ; TODO fix
   (test (deriv '(- x y) 'z) '(- 0 0)) ; TODO fix

   (test (deriv '(* x x) 'x) '(+ x x))
   (test (deriv '(* x y) 'x) 'y)
   (test (deriv '(* x 3) 'x) '3)
   (test (deriv '(* 2 3) 'x) '0)

   (test (deriv '(/ 1 1) 'x) '(- 0 0)) ; TODO fix
   (test (deriv '(/ x 1) 'x) '1)
   (test (deriv '(/ x y) 'x) '(/ y (* y y))) ; TODO fix

   (test (deriv '(+ x x x) 'x) '3)
   (test (deriv '(+ x y z) 'z) '1)
   (test (deriv '(+ x z y) 'z) '1)
   (test (deriv '(+ 1 2 z) 'z) '1)

   (test (deriv '(+ (+ x y) z) 'z) '1)
   (test (deriv '(+ (+ x y) (+ x y)) 'x) '2)
   (test (deriv '(+ (- x y) (+ x y)) 'x) '2)
   (test (deriv '(+ (- x y) (+ x y)) 'y) '(+ (- 0 1) 1)) ; TODO fix
   (test (deriv '(+ (- x y) (- x y)) 'x) '2)
   (test (deriv '(+ (- x y) (- x y)) 'y) '(+ (- 0 1) (- 0 1))) ; TODO fix
   (test (deriv '(- (- x y) (- x y)) 'y) '(- (- 0 1) (- 0 1))) ; TODO fix
   (test (deriv '(* (* x y) (+ x 3)) 'x) '(+ (* x y) (* y (+ x 3))))
   (test (deriv '(/ (* x y) (+ x 3)) 'x) '(/ (- (* y (+ x 3)) (* x y)) (* (+ x 3) (+ x 3))))
    
   (test (deriv '(+ (+ x x x) (+ x x y)) 'x) '5)
   ))
(run-tests suite)

                    
(define (df f var)
  (eval (deriv f var) ie))

(define suite
  (list
   (test (df '1 'x) 0)
   (test (df 'y 'x) 0)
   (test (df 'x 'x) 1)
   
   (test (df '(+ x y) 'x) 1)
   (test (df '(+ x 3) 'x) 1)
   (test (df '(+ y x) 'x) 1)
   (test (df '(+ x x) 'x) 2)
   ;(test (deriv '(- x y) 'y) -1)
   ;(test (deriv '(- x y) 'z) 0)
   
   (test (df '(+ x x x) 'x) 3)
   ))
(run-tests suite)

