; =================== Macro ====================
(display "============ Macro ==============")
(newline)


(define-syntax when
  (syntax-rules ()
    ((_ cond-expr expr . exprs)
     (if cond-expr
         (begin
           expr . exprs)))))

(define m 3)
(when (< 0 m)
  (display "when")
  (newline)
  (set! m (- m 1)))


(define-syntax unless
  (syntax-rules ()
    ((_ cond-expr expr . exprs)
     (if (not cond-expr)
         (begin
           expr . exprs)))))

(define m 3)
(unless (> 0 m)
  (display "unless")
  (newline))



(define (each proc xs)
  (if (not (null? xs))
      (begin
        (proc (car xs))
        (each proc (cdr xs)))))

(define-syntax for
  (syntax-rules (in as)
    ((_ _ in items expr) (each expr items))
    ((_ items as _ expr) (each expr items))))

(for i in '(1 2 3)
  (lambda (i) (display i)))

(for '(4 5 6) as j
  (lambda (j) (display j)))
