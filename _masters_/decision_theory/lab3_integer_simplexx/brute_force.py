import numpy as np
import itertools

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, Condition
from _masters_.decision_theory.lab3_integer_simplexx.bnb import override


class BruteForce(Simplexx):
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, condition: Condition):
        """
        in form of AX<=B. Not >=, not ==
        :param a: A
        :param b: B
        :param lambdas: aka C
        :param condition: min or max
        """
        super().__init__(a, b, lambdas, condition)

    @override
    def run(self) -> dict:
        # количество x, присутствующих в F == кол-ву лямбд
        number_of_considered_xs = self.lambdas.shape[1]
        ranges_to_brute = []
        for x in range(number_of_considered_xs):
            range_from, range_to = self.range_to_brute(x)
            print(f'x_{x + 1}: диапазон перебора [{range_from}..{range_to}]')

            range_list = list(range(range_from, range_to+ 1))
            ranges_to_brute.append(range_list)

        # собираем декартово произведение. получаем набор x-векторов
        cartesian_products = list(itertools.product(*ranges_to_brute))

        # оставляем только те решения (x-векторы), которые удовлетворяют условию
        applicable_x_vectors = list(filter(lambda solution: self.is_applicable(solution),
                                           cartesian_products))

        # находим лучшее (max или min) значение F
        best_solution_vect = None
        best_f_value = -99999
        for x_vect in applicable_x_vectors:
            xs = np.array(x_vect)
            # вычисление F
            f_value = np.dot(self.lambdas, xs).item()

            if self.condition is Condition.MAX:
                if f_value > best_f_value:
                    best_f_value = f_value
                    best_solution_vect = x_vect
            else:
                if f_value < best_f_value:
                    best_f_value = f_value
                    best_solution_vect = x_vect

        return BruteForce.as_result(best_f_value, best_solution_vect)


    # найти диапазон, в котором будет вариьироваться 'x_3'. x_id=3
    def range_to_brute(self, x_id: int) -> (int, int):
        # нумерация x'ов - с единицы. нумерация строк - с нуля
        x_id -= 1

        # находим максимальное отношение b к x_i
        maximal_b_to_x_relation = 0
        for i in range(self.matrix_rows):
            # деление на 0
            if 0 == self.matr[i, x_id]:
                i += 1
                continue

            relation = self.b[i, 0] / self.matr[i, x_id]
            if relation > maximal_b_to_x_relation:
                maximal_b_to_x_relation = relation

        # диапазон всегда идет от нуля
        return 0, int(maximal_b_to_x_relation.item())

    # решение удовлетворяет услвоию ?
    def is_applicable(self, solution_vect: list) -> bool:
        xs = np.array(solution_vect).transpose()
        # [[A]] * [x]^T
        result = np.dot(self.matr, xs)

        for i in range(len(self.b)):
            # условие записано, как Ax <= b. теперь мы знаем 'x' и проверяем, что это так
            if result[i] > self.b[i]:
                return False
        return True

    @staticmethod
    def as_result(f_value: float, x_vect: list) -> dict:
        res = dict()
        for i in range(len(x_vect)):
            x_key = f'x_{i + 1}'
            res[x_key] = x_vect[i]

        res['F'] = f_value
        return res

    @property
    def matrix_rows(self):
        return self.matr.shape[0]

    @property
    def matrix_cols(self):
        return self.matr.shape[1]


if __name__ == '__main__':
    a = np.array([[6, -1],
                  [2, 5]])
    b = np.array([[12],
                  [20]])
    lambdas = np.array([[12, -1]])

    solution = BruteForce(a, b, lambdas, Condition.MAX).run()
    print(f'\n\nИтоговое решение: {solution}')
