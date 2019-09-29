import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, Condition


class DualSimplexx(Simplexx):
    """
    Двойственность в линейном программировании
    """
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, condition: Condition):
        """
        in form of AX<=B.
        Not >=, not ==
        :param a: A
        :param b: B
        :param lambdas: aka C
        :param condition: min or max
        """
        # чтобы перейти к ДЗЛП, надо
        # - транспонировать A, B, C.
        # - поменять знаки >= -> <= путем домножения левой и правой части на -1
        # - инвертировать условие MAX -> MIN
        a_t = -1 * np.transpose(a)
        lambdas_t = np.transpose(b)
        b_t = -1 * np.transpose(lambdas)
        cond_inverted = DualSimplexx._invert_condition(condition)
        super().__init__(a_t, b_t, lambdas_t, cond_inverted)

    @staticmethod
    def _invert_condition(cond: Condition) -> Condition:
        if cond is Condition.MAX:
            return Condition.MIN

        elif cond is Condition.MIN:
            return Condition.MAX


def main():
    a = np.array([[3, 1, -4, -1],
                  [-2, -4, -1, 1]])
    b = np.array([[-3],
                  [-3]])
    lambdas = np.array([[-4, -18, -30, -5]])

    s = DualSimplexx(a, b, lambdas, Condition.MAX)
    solution = s.run()

    print('\nОтвет:')
    print(solution)


if __name__ == '__main__':
    main()
