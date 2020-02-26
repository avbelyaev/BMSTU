; Env is environment <x, n, res>
(declare-datatypes ()
    ((Env (mk-env
        (env-pi (Array Int Int))
        (env-arr (Array Int Int))
        (env-len-arr Int)
        (env-i Int)
        (env-t Int)
        ))))

; === FUNCTIONS ===
(define-fun condLoopOverArray((e Env)) Bool
	(< (env-i e) (env-len-arr e)))

(define-fun statementA((e Env)) Env
	(let ((pi (env-pi e))
	      (arr (env-arr e))
	      (len-arr (env-len-arr e))
	      (i (env-i e))
	      (t (env-t e)))
		(mk-env pi arr len-arr i (select pi (- i 1)))))

(define-fun condLoopOverPrefix((e Env)) Bool
	(and
	    (> (env-t e) 0)
	    (not (= (select (env-arr e) (env-t e)) (select (env-arr e) (env-i e))))))

(define-fun statementB((e Env)) Env
	(let ((pi (env-pi e))
	      (arr (env-arr e))
	      (len-arr (env-len-arr e))
	      (i (env-i e))
	      (t (env-t e)))
		(mk-env pi arr len-arr i (select pi (- t 1)))))

(define-fun condIf((e Env)) Bool
	(= (select (env-arr e) (env-t e)) (select (env-arr e) (env-i e))))

(define-fun statementC((e Env)) Env
	(let ((pi (env-pi e))
	      (arr (env-arr e))
	      (len-arr (env-len-arr e))
	      (i (env-i e))
	      (t (env-t e)))
		(mk-env pi arr len-arr i (+ t 1))))

(define-fun statementD((e Env)) Env
	(let ((pi (env-pi e))
	      (arr (env-arr e))
	      (len-arr (env-len-arr e))
	      (i (env-i e))
	      (t (env-t e)))
		(mk-env (store pi i t) arr len-arr (+ i 1) t)))

; =============================================================================
; Definitions shared by all models
; -----------------------------------------------------------------------------

(declare-const PI (Array Int Int))
(declare-const ARR (Array Int Int))

; =============================================================================
; Model 1: A D
; -----------------------------------------------------------------------------

(push)
(define-const e1 Env (mk-env PI ARR 2 1 0))

; i = 1
(assert (condLoopOverArray e1))
(define-const e2 Env (statementA e1))
(assert (not (condLoopOverPrefix e2)))
(assert (not (condIf e2)))
(define-const e3 Env (statementD e2))
; i = 2
(assert (not (condLoopOverArray e3)))

(echo "")
(echo "Model 1: A D")
(check-sat)
(eval (select ARR 0))
(eval (select ARR 1))
(get-model)
(pop)


; =============================================================================
; Model 2: A C D
; -----------------------------------------------------------------------------

(push)
(define-const e1 Env (mk-env PI ARR 2 1 0))

; i = 1
(assert (condLoopOverArray e1))
(define-const e2 Env (statementA e1))
(assert (not (condLoopOverPrefix e2)))
(assert (condIf e2))
(define-const e3 Env (statementC e2))
(define-const e4 Env (statementD e3))
; i = 2
(assert (not (condLoopOverArray e4)))

(echo "")
(echo "Model 2: A C D")
(check-sat)
(eval (select ARR 0))
(eval (select ARR 1))
(get-model)
(pop)


; =============================================================================
; Model 3: A D A C D
; -----------------------------------------------------------------------------

(push)
(define-const e1 Env (mk-env PI ARR 3 1 0))

; i = 1
(assert (condLoopOverArray e1))
(define-const e2 Env (statementA e1))
(assert (not (condLoopOverPrefix e2)))
(assert (not (condIf e2)))
(define-const e3 Env (statementD e2))
; i = 2
(assert (condLoopOverArray e3))
(define-const e4 Env (statementA e3))
(assert (not (condLoopOverPrefix e4)))
(assert (condIf e4))
(define-const e5 Env (statementC e4))
(define-const e6 Env (statementD e5))
; i = 3
(assert (not (condLoopOverArray e6)))

(echo "")
(echo "Model 3: A D A C D")
(check-sat)
(eval (select ARR 0))
(eval (select ARR 1))
(eval (select ARR 2))
(get-model)
(pop)


; =============================================================================
; Model 4: A D A D A C D
; -----------------------------------------------------------------------------

(push)
(define-const e1 Env (mk-env PI ARR 4 1 0))

; i = 1
(assert (condLoopOverArray e1))
(define-const e2 Env (statementA e1))
(assert (not (condLoopOverPrefix e2)))
(assert (not (condIf e2)))
(define-const e3 Env (statementD e2))
; i = 2
(assert (condLoopOverArray e3))
(define-const e4 Env (statementA e3))
(assert (not (condLoopOverPrefix e4)))
(assert (not (condIf e4)))
(define-const e5 Env (statementD e4))
; i = 3
(assert (condLoopOverArray e5))
(define-const e6 Env (statementA e5))
(assert (not (condLoopOverPrefix e6)))
(assert (condIf e6))
(define-const e7 Env (statementC e6))
(define-const e8 Env (statementD e7))
; i = 4
(assert (not (condLoopOverArray e8)))

(echo "")
(echo "Model 4: A D A D A C D")
(check-sat)
(eval (select ARR 0))
(eval (select ARR 1))
(eval (select ARR 2))
(eval (select ARR 3))
(get-model)
(pop)
