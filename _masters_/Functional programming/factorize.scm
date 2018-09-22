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



; ---------------- quad ---------------------
(define (quad? op expr)
  (and (list? expr)
       (equal? op (car expr))
       (expt-x-2? (cadr expr))
       (expt-x-2? (caddr expr))))

(define suite
  (list
   (test (quad? '- '(- (expt a 2) (expt b 2))) #t)
   (test (quad? '+ '(+ (expt a 2) (expt b 2))) #t)
   (test (quad? '- '(- (expt x 2) (expt y 2))) #t)
   
   (test (quad? '- '(- (expt a 0) (expt b 2))) #f)
   (test (quad? '- '(- (expt a 2) (expt b 0))) #f)
   (test (quad? '- '(- (expt a 0) (expt b 0))) #f)
   (test (quad? '+ '(- (expt a 2) (expt b 2))) #f)
   (test (quad? '- '(- (exp a 2) (exp b 2))) #f)))
(run-tests suite)



(define (quad-sub? expr)
  (quad? '- expr))

(define suite
  (list
   (test (quad-sub? '(- (expt a 2) (expt b 2))) #t)
   (test (quad-sub? '(- (expt x 2) (expt y 2))) #t)
   
   (test (quad-sub? '(+ (expt a 2) (expt b 2))) #f)
   (test (quad-sub? '(- (expt a 3) (expt b 2))) #f)
   (test (quad-sub? '(- (expt a 2) (expt b 3))) #f)
   (test (quad-sub? '(- (exp a 2) (exp b 2))) #f)))
(run-tests suite)



(define (quad-add? expr)
  (quad? '+ expr))

(define suite
  (list
   (test (quad-add? '(+ (expt a 2) (expt b 2))) #t)
   (test (quad-add? '(+ (expt x 2) (expt y 2))) #t)
   
   (test (quad-add? '(- (expt a 2) (expt b 2))) #f)
   (test (quad-add? '(+ (expt a 3) (expt b 2))) #f)
   (test (quad-add? '(+ (expt a 2) (expt b 3))) #f)
   (test (quad-add? '(+ (exp a 2) (exp b 2))) #f)))
(run-tests suite)



; --------------- cube ------------------
(define (cube? op expr)
  (and (list? expr)
       (equal? op (car expr))
       (expt-x-3? (cadr expr))
       (expt-x-3? (caddr expr))))

(define suite
  (list
   (test (cube? '- '(- (expt a 3) (expt b 3))) #t)
   (test (cube? '+ '(+ (expt a 3) (expt b 3))) #t)
   (test (cube? '- '(- (expt x 3) (expt y 3))) #t)
   
   (test (cube? '+ '(- (expt a 3) (expt b 3))) #f)
   (test (cube? '- '(- (expt a 0) (expt b 3))) #f)
   (test (cube? '- '(- (expt a 3) (expt b 0))) #f)
   (test (cube? '- '(- (expt a 0) (expt b 0))) #f)
   (test (cube? '- '(- (exp a 3) (exp b 3))) #f)))
(run-tests suite)



(define (cube-sub? expr)
  (cube? '- expr))

(define suite
  (list
   (test (cube-sub? '(- (expt a 3) (expt b 3))) #t)
   (test (cube-sub? '(- (expt x 3) (expt y 3))) #t)
   
   (test (cube-sub? '(+ (expt a 3) (expt b 3))) #f)
   (test (cube-sub? '(- (expt a 0) (expt b 3))) #f)
   (test (cube-sub? '(- (expt a 3) (expt b 0))) #f)
   (test (cube-sub? '(- (expt a 0) (expt b 0))) #f)
   (test (cube-sub? '(- (exp a 3) (exp b 3))) #f)))
(run-tests suite)



(define (cube-add? expr)
  (cube? '+ expr))

(define suite
  (list
   (test (cube-add? '(+ (expt a 3) (expt b 3))) #t)
   (test (cube-add? '(+ (expt x 3) (expt y 3))) #t)
   
   (test (cube-add? '(- (expt a 3) (expt b 3))) #f)
   (test (cube-add? '(+ (expt a 0) (expt b 3))) #f)
   (test (cube-add? '(+ (expt a 3) (expt b 0))) #f)
   (test (cube-add? '(+ (expt a 0) (expt b 0))) #f)
   (test (cube-add? '(+ (expt a 0) (expt b 0))) #f)
   (test (cube-add? '(+ (exp a 3) (exp b 3))) #f)))
(run-tests suite)



;---------- factorize ----------
(define (factorize expr)
  (let ((a (cadr (cadr expr)))
        (b (cadr (caddr expr))))
    (cond ((quad-sub? expr) `(* (- ,a ,b) (+ ,a ,b)))
          ((quad-add? expr) `(+ (* ,a ,a) (* ,b ,b)))
          ((cube-sub? expr) `(* (- ,a ,b ) (+ (* ,a ,a) (* ,a ,b) (* ,b ,b))))
          ((cube-add? expr) `(* (+ ,a ,b ) (+ (- (* ,a ,a) (* ,a ,b)) (* ,b ,b)))))))

(define suite
  (list
   (test (eval (factorize '(- (expt 1 2) (expt 1 2))) (interaction-environment)) 0)
   (test (eval (factorize '(- (expt -1 2) (expt -1 2))) (interaction-environment)) 0)
   (test (eval (factorize '(+ (expt 1 2) (expt 1 2))) (interaction-environment)) 2)
   
   (test (eval (factorize '(- (expt 1 3) (expt 1 3))) (interaction-environment)) 0)
   (test (eval (factorize '(+ (expt 1 3) (expt 1 3))) (interaction-environment)) 2)
   (test (eval (factorize '(+ (expt -1 3) (expt -1 3))) (interaction-environment)) -2)
   
   (test (eval (factorize '(- (expt 5 2) (expt 4 2))) (interaction-environment)) 9)
   (test (eval (factorize '(+ (expt 5 2) (expt 4 2))) (interaction-environment)) 41)
   
   (test (eval (factorize '(- (expt 3 3) (expt 2 3))) (interaction-environment)) 19)
   (test (eval (factorize '(+ (expt 3 3) (expt 2 3))) (interaction-environment)) 35)))
(run-tests suite)
