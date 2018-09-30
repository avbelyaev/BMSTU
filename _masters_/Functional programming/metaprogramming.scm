(define ie (interaction-environment))

(define (make-name cite postfix) ; symbol-list, str
  (string->symbol (string-append (symbol->string (car cite)) postfix)))


(define (define-struct cite)
  (let ((predicate `(define (,(make-name cite "?") arg)
                 (equal? ,(symbol->string (car cite)) arg)))
        (anti-pred `(define (,(make-name cite "!") arg)
                  (not (equal? ,(symbol->string (car cite)) arg)))))
    (display "eval1")
    (newline)
    (eval predicate ie)
    
    (display "eval2")
    (newline)
    (eval anti-pred ie)))


(define-struct '(point "ignored-arg"))

(point? "1253")
(point? "point")


(display "anti-pred")(newline)

(point! "1235")
(point! "point")


