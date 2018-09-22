(define (and-fold . args)
  (if (null? args)
      #t
      (if (car args)
          (apply and-fold (cdr args))
          #f)))

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


(define examples
  (list
   (test 1 1)
   (test "abc" "abc")
   (test (abs -2) 2)
   (test (+ 1 2 3) 5)))


(run-tests examples)
