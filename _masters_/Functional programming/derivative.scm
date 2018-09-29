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




; ------- factorize cubes and quads --------
(display "=============== Factorize ===============")
(newline)


(define (is-var? var)
  (or (symbol? var)
      (number? var)
      (list? var)))

(define suite
  (list
   (test (is-var? 'a) #t)
   (test (is-var? 12) #t)
   (test (is-var? '(1 2 3)) #t)
   (test (is-var? '(expt a 2)) #t)))
(run-tests suite)



; ------------- expt x n --------------
(define (expt-x-n? expr n)
  (and (list? expr)
       (equal? (car expr) 'expt)
       (is-var? (cadr expr))
       (equal? n (caddr expr))))

(define suite
  (list
   (test (expt-x-n? '(expt a 2) 2) #t)
   (test (expt-x-n? '(expt b 3) 3) #t)
   (test (expt-x-n? '(expt a 3) 0) #f)))
(run-tests suite)



(define (expt-x-2? expr)
  (expt-x-n? expr 2))

(define suite
  (list
   (test (expt-x-2? '(expt a 2)) #t)
   (test (expt-x-2? '(expt b 2)) #t)
   (test (expt-x-2? '(expt a 3)) #f)))
(run-tests suite)



(define (expt-x-3? expr)
  (expt-x-n? expr 3))

(define suite
  (list
   (test (expt-x-3? '(expt a 3)) #t)
   (test (expt-x-3? '(expt b 3)) #t)
   (test (expt-x-3? '(expt a 0)) #f)))
(run-tests suite)



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



(define (make-sum a b)
  `(+ a b))

(define (make-prod a b)
  `(* a b))



(define (sum? x)
  (and (pair? x)
       (eq? (car x) '+)))

(define (prod? x)
  (and (pair? x)
       (eq? (car x) '*)))

(define suite
  (list
   (test (sum? '(+ 1 2)) #t)
   (test (sum? '(+ a b)) #t)
   (test (sum? '(+ a b c)) #t)
   (test (sum? '(* a b)) #f)

   (test (prod? '(* 1 2)) #t)
   (test (prod? '(* a b)) #t)
   (test (prod? '(* a b c)) #t)
   (test (prod? '(+ a b)) #f)))
(run-tests suite)



(define (addend s)
  (cadr s))

(define suite
  (list
   (test (addend '(+ 1 2)) 1)
   (test (addend '(+ 3 2 1)) 3)
   (test (addend '(+ a b c)) 'a)))
(run-tests suite)



(define (augend s)
  (define (helper xs)
    (if (null? xs)
        '()
        (cond (car xs) (helper (cdr xs)))))
  (if (= 3 (length s)) ; (+ 1 2) => 2
      (caddr s)
      (cons '+ (helper (cddr s)))))

(define suite
  (list
   (test (augend '(+ 1 2)) '2)
   (test (augend '(+ 1 2 3)) '(+ 2 3))
   (test (augend '(+ a b c)) '(+ b c))
   (test (augend '(+ a b c d)) '(+ b c d))))
(run-tests suite)




(define (make-add a1 a2)
    (cond ((eq-number? a1 0) a2)
          ((eq-number? a2 0) a1)
          ((and (number? a1) (number? a2)) (+ a1 a2))
          (else (list '+ a1 a2))))

(define suite
  (list
   (test (make-add 'x 0) 'x)
   (test (make-add 0 'x) 'x)
   (test (make-add 1 2) 3)
   (test (make-add 'x 'x) '(+ x x))
   (test (make-add 'x 'y) '(+ x y))
   ))
(run-tests suite)



(define (deriv f var)
  (cond ((number? f) 0)
        ((var? f) (if (same-var? f var) 1 0))
        ((list f) (let ((op (car f)))
                     (case op
                       ('+ (make-add (deriv (addend f) var)
                                     (deriv (augend f) var)))
                       ('- (make-sub f var)))))))

(define suite
  (list
   (test (deriv '1 'x) '0)
   (test (deriv 'x 'x) '1)
   (test (deriv 'y 'x) '0)
   (test (deriv '(+ x x) 'x) '2)
   (test (deriv '(+ x y) 'x) '1)
   (test (deriv '(+ x y) 'y) '1)
   (test (deriv '(+ x y) 'z) '0)

   (test (deriv '(+ x x x) 'x) '3)
   (test (deriv '(+ x y z) 'z) '1)
   (test (deriv '(+ x z y) 'z) '1)

   (test (deriv '(+ (+ x y) z) 'z) '1)
   (test (deriv '(+ (+ x y) (+ x y)) 'x) '2)
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
   
   (test (df '(+ x x x) 'x) 3)
   ))
(run-tests suite)

