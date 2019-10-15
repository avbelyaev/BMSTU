from enum import Enum
from typing import Optional

import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, Condition, \
    NoPivotalSolutionExists, NoOptimalSolutionExists


class NoIntegerSolutionExists(Exception):
    def __init__(self):
        super().__init__('Целочиселнного решения не существует!')


class Sign(Enum):
    MINUS = 0
    PLUS = 1


def override(f):
    return f


class BranchAndBound:
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, condition: Condition):
        assert condition in [Condition.MIN, Condition.MAX]
        # super().__init__(a, b, lambdas, condition)
        self.matr = a
        self.b = b
        self.lambdas = lambdas
        self.condition = condition

    def first_non_integer_solution(self, solution: dict) -> Optional[dict]:
        is_int = lambda val: isinstance(val, int) \
                             or (isinstance(val, float) and val.is_integer())

        # нам важны только те x_i, которые участвуют в функции
        valuable_xs = [f'x_{i + 1}' for i in range(self.lambdas.shape[1])]
        for x_coord in solution.keys():
            if x_coord in valuable_xs and not is_int(solution[x_coord]):
                return {x_coord: solution[x_coord]}
        return None

    @override
    def run(self):
        parent_simplex = Simplexx(self.matr, self.b, self.lambdas, self.condition)
        parent_solution = parent_simplex.run()
        print(f'Решение: {parent_solution}')

        non_int_value = self.first_non_integer_solution(parent_solution)
        if non_int_value is not None:
            print(f'Нецелочисленное решение: {non_int_value}. Разветвляем')
            branch1, branch2 = BranchAndBound.split(parent_simplex, non_int_value)

            try:
                print('\n\nРешаем ветку 1')
                child1_solution = branch1.run()
            except (NoPivotalSolutionExists, NoOptimalSolutionExists) as e:
                print(e)
                child1_solution = None

            try:
                print('\n\nРешаем ветку 2')
                child2_solution = branch2.run()
            except (NoPivotalSolutionExists, NoOptimalSolutionExists) as e:
                print(e)
                child2_solution = None

            # обе ветки решились
            if child1_solution is not None and child2_solution is not None:
                if self.condition is Condition.MAX:
                    if child1_solution['F'] > child2_solution['F']:
                        return child1_solution
                    else:
                        return child2_solution

                else:
                    if child1_solution['F'] < child2_solution['F']:
                        return child1_solution
                    else:
                        return child2_solution

            # ветка 1 не решилась, 2 - решилась
            elif child1_solution is None and child2_solution is not None:
                return child2_solution

            # ветка 1 решилась, 2 - не решилась
            elif child1_solution is not None and child2_solution is None:
                return child1_solution

            # обе ветки не решились
            else:
                raise NoIntegerSolutionExists()

        else:
            print('Решение целочисленное. Ветка закончена')
            return parent_solution

    @staticmethod
    def split(simplex: Simplexx, non_int_x: dict) -> 'BranchAndBound, BranchAndBound':
        # { x_3: 5.5 } -> (x_3, 5, 6)
        def split_x_value(non_int_kv_pair: dict) -> (str, int, int):
            key = next(iter(non_int_kv_pair.keys()))
            value = non_int_kv_pair[key]
            return key, int(value), int(value) + 1

        x_key, x_less_value, x_greater_value = split_x_value(non_int_x)

        # x <= non_int
        # e.g: x3 <= 5
        simplex1 = BranchAndBound.create_simplex_with_additional_bound(simplex, x_key, x_less_value, Sign.PLUS)
        bnb1 = BranchAndBound(simplex1.matr, simplex1.b, simplex1.lambdas, simplex1.condition)

        # x >= non_int + 1
        # e.g: x3 >= 6 <-(то же самое, но со знаком минус)-> -x3 <= -6
        simplex2 = BranchAndBound.create_simplex_with_additional_bound(simplex, x_key, x_greater_value, Sign.MINUS)
        bnb2 = BranchAndBound(simplex2.matr, simplex2.b, simplex2.lambdas, simplex2.condition)
        return bnb1, bnb2

    # добавляем к текущей задаче новое ограничение "x_3 <= 5" в виде
    # - дополнительой строки матрицы A: [0, 0, 1], где 1 отвечает за x_3
    # - и дополнительной строки столбца b: [5]
    @staticmethod
    def create_simplex_with_additional_bound(simplex: Simplexx,
                                             new_bound_key: str,
                                             new_bound_value: int,
                                             new_bound_sign: Sign) -> Simplexx:
        # находим колонку, которая соответстует нецелому "x_i". нумерация x - с единицы
        new_bound_column = int(new_bound_key.split('_')[1]) - 1

        matr_columns = simplex.matr.shape[1]
        additional_a_row = [0] * matr_columns
        if new_bound_sign is Sign.PLUS:
            additional_a_row[new_bound_column] = 1
        else:
            additional_a_row[new_bound_column] = -1
        additional_a_row = np.array(additional_a_row)

        # дпоисываем строку в низ матрицы А
        curr_a = simplex.matr
        new_a = np.vstack((curr_a, additional_a_row))

        additional_b_row = None
        if new_bound_sign is Sign.PLUS:
            additional_b_row = np.array([new_bound_value])
        else:
            additional_b_row = np.array([-1 * new_bound_value])

        # дописываем значение ограничения в низ столбца b
        curr_b = simplex.b
        new_b = np.vstack((curr_b, additional_b_row))

        return Simplexx(new_a, new_b, simplex.lambdas, simplex.condition)


if __name__ == '__main__':
    a = np.array([[6, -1],
                  [2, 5]])
    b = np.array([[12],
                  [20]])
    lambdas = np.array([[12, -1]])

    solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()
    print(f'\n\nИтоговое решение: {solution}')
