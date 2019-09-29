from typing import Optional

import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, Condition


def override(f):
    return f


class BranchAndBound(Simplexx):
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, condition: Condition):
        super().__init__(a, b, lambdas, condition)

    @override
    def run(self):
        solution = super().run()
        non_int = self.first_non_integer_solution(solution)

        if non_int is not None:
            # branch 1
            a1 = np.array([[6, -1],
                          [2, 5],
                          [1, 0]])
            b1 = np.array([[12],
                          [20],
                          [2]])
            b1_solution = Simplexx(a1, b1, self.lambdas, self.condition).run()
            print(b1_solution)

            # branch 2
            a2 = np.array([[6, -1],
                          [2, 5],
                          [1, 0]])
            b2 = np.array([[12],
                          [20],
                          [-3]])
            try:
                b2_solution = Simplexx(a2, b2, self.lambdas, self.condition).run()
                print(b2_solution)
            except Exception as e:
                print(f'error: {e}')


    def first_non_integer_solution(self, solution: dict) -> Optional[dict]:
        for x in solution.keys():
            if x != 'F':
                if not isinstance(solution[x], int):
                    return {x: solution[x]}
        return None


if __name__ == '__main__':
    a = np.array([[6, -1],
                  [2, 5]])
    b = np.array([[12],
                  [20]])
    lambdas = np.array([[12, -1]])

    # when
    bnb = BranchAndBound(a, b, lambdas, Condition.MAX)
    primary_solution = bnb.run()
    print(primary_solution)
