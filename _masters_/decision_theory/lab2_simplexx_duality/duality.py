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
        if condition is Condition.MAX:
            a_t = -1 * np.transpose(a)
            lambdas_t = np.transpose(b)
            b_t = -1 * np.transpose(lambdas)
            cond = Condition.MIN

        elif condition is Condition.MIN:
            a_t = -1 * np.transpose(a)
            lambdas_t = -1 * np.transpose(b)
            b_t = np.transpose(lambdas)
            cond = Condition.MAX
        else:
            raise ValueError('Неверное условие!')

        super().__init__(a_t, b_t, lambdas_t, cond)


if __name__ == '__main__':
    print('см. тесты')
    pass
