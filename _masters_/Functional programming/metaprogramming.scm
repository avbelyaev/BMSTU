(define (define-struct cite)
  (let ((predicate `(define (,(string->symbol (string-append (symbol->string (car cite)) "?")) arg)
                 (eq? ,(car cite) arg))) ; strinbg .comare to symbol
        (anti-pred `(define (,(string->symbol (string-append (symbol->string (car cite)) "!")) arg)
                  (not (eq? ,(car cite) arg)))))

    (display "eval1")
    (eval predicate (interaction-environment))
    (display "eval2")
    (eval anti-pred (interaction-environment))))

(define-struct '(point "321"))

(point? "1253")
(point? "point")


(display "anti-pred")(newline)

(point! "1235")
(point! "point")

