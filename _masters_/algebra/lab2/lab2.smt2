(declare-const X Int)
(declare-const Y Int)
(assert (not (=> (= X (+ Y 1)) (not (= (+ X 1) Y)))))
(check-sat)