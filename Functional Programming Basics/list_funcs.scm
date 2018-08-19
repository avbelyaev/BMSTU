(define (my-range a b d)
  (define (iter a b d xs)
    (if (< a b)
        (iter (+ a d) b d (append xs (list a)))
        xs))
  (iter a b d '()))
 
(define (my-flatten xs)
  (define (iter xs ys)
    (if (null? xs)
        ys
        (if (list? (car xs))
            (iter (cdr xs) (append ys (iter (car xs) '())))
            (iter (cdr xs) (append ys (list (car xs)))))))
  (iter xs '()))
 
(define (my-element? x xs)
  (define (f x xs k)
    (if (null? xs)
        (not k)
        (if (equal? (car xs) x)
            k
            (f x (cdr xs) k))))
  (f x xs #t))
 
(define (my-filter pred? xs)
  (define (iter pred? xs ys)
    (if (null? xs)
        ys
        (if (pred? (car xs))
            (iter pred? (cdr xs) (append ys (list (car xs))))
            (iter pred? (cdr xs) ys))))
  (iter pred? xs '()))
 
(define (my-fold-left op xs)
  (define (iter op xs count)
    (if (null? xs)
        count
        (iter op (cdr xs) (op count (car xs)))))
  (iter op (cdr xs) (car xs)))
 
(define (my-fold-right op xs)
  (define (iter op xs count)
    (if (null? xs)
        count
        (iter op (cdr xs) (op (car xs) count))))
  (iter op (cdr (reverse xs)) (car (reverse xs))))

