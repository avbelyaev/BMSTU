
; =============================================================================
; Python function returning x^n.
; -----------------------------------------------------------------------------
; def power(x, n):
;   res = 1
;   while n != 0:              # condLoop
;     if n % 2 == 0:           # condIf
;       print "A"
;       x, n = x*x, n/2        # statementA
;     else:
;       print "B"
;       res, n = res * x, n-1  # statementB
;   return res

; Env is environment <x, n, res>
(declare-datatypes () ((Env (mk-env (env-x Int) (env-n Int) (env-res Int)))))

(define-fun condLoop((e Env)) Bool
	(not (= (env-n e) 0)))

(define-fun condIf((e Env)) Bool
	(= (mod (env-n e) 2) 0))

(define-fun statementA((e Env)) Env
	(let ((x (env-x e)) (n (env-n e)) (res (env-res e)))
		(mk-env (* x x) (div n 2) res)))

(define-fun statementB((e Env)) Env
	(let ((x (env-x e)) (n (env-n e)) (res (env-res e)))
		(mk-env x (- n 1) (* res x))))

; =============================================================================
; Definitions shared by all models
; -----------------------------------------------------------------------------

(declare-const N Int)
(declare-const X Int)
(assert (>= N 0))
(define-const e1 Env (mk-env X N 1))

; =============================================================================
; Model 1: A B
; -----------------------------------------------------------------------------

(push)
(assert (and (condLoop e1) (condIf e1)))
(define-const e2 Env (statementA e1))
(assert (and (condLoop e2) (not (condIf e2))))
(define-const e3 Env (statementB e2))
(assert (not (condLoop e3)))

(echo "A B")
(check-sat)
(get-model)
(pop)

; =============================================================================
; Model 2: A A B A B
; -----------------------------------------------------------------------------

(push)
(assert (and (condLoop e1) (condIf e1)))
(define-const e2 Env (statementA e1))
(assert (and (condLoop e2) (condIf e2)))
(define-const e3 Env (statementA e2))
(assert (and (condLoop e3) (not (condIf e3))))
(define-const e4 Env (statementB e3))
(assert (and (condLoop e4) (condIf e4)))
(define-const e5 Env (statementA e4))
(assert (and (condLoop e5) (not (condIf e5))))
(define-const e6 Env (statementB e5))
(assert (not (condLoop e6)))

(echo "A A B A B")
(check-sat)
(get-model)
(pop)
