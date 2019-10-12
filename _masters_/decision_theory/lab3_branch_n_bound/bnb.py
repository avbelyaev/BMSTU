from enum import Enum
from typing import Optional

import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, Condition, \
    NoPivotalSolutionExists, NoOptimalSolutionExists


class Sign(Enum):
    MINUS = 0
    PLUS = 1


def override(f):
    return f


class Task:
    def __init__(self, task_index: int, simplex: Simplexx):
        self.k = task_index
        self.simplex = simplex
        try:
            self.solution = simplex.run()
            self.f_value = self.solution['F']
        except (NoPivotalSolutionExists, NoOptimalSolutionExists):
            self.solution = None
            self.f_value = None

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'Task[{self.k}]({self.simplex})'


class BranchAndBound:
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, condition: Condition):
        # super().__init__(a, b, lambdas, condition)
        self.matr = a
        self.b = b
        self.lambdas = lambdas
        self.condition = condition
        self.tasks = []

    def first_non_integer_solution(self, solution: dict) -> Optional[dict]:
        is_int = lambda val: isinstance(val, int) \
                             or (isinstance(val, float) and val.is_integer())
        for x_coord in solution.keys():
            if x_coord in ['x_1', 'x_2', 'x_3']:
                # FIXME берем только x \in [x1, x2, x3]. x4 не берем
                if not is_int(solution[x_coord]):
                    return {x_coord: solution[x_coord]}
        return None

    @override
    def run(self):
        # Ш.1.
        k = 0
        task = Task(k, Simplexx(self.matr, self.b, self.lambdas, self.condition))
        non_int_value = self.first_non_integer_solution(task.solution)
        print(f'Нецелочисленное решение: {non_int_value}')

        # если решение нецелочисленное, включить k = 0 в множество J = {k} номеров задач,
        # подлежащих дальнейшему ветвлению, и перейти на шаг 2
        if non_int_value is not None:
            self.tasks.insert(k, task)

            # Ш.2. Выбрать задачу для приоритетного ветвления:
            task_to_derieve = None
            if 0 == k:
                # выбрать для ветвления задачу ЗЛП-0,
                # исключить номер k из множества J = {k} и перейти на Ш.3.
                task_to_derieve = self.tasks[k]
                del self.tasks[k]

            elif 0 != k and 0 != len(self.tasks):
                # выбрать номер задачи, которму соответствует
                # максимальное значение целевой функции и перейти на Ш.3.
                task_with_max_f_value = next(sorted(self.tasks,
                                                    key=lambda tsk: tsk.f_value,
                                                    reverse=True),
                                             None)  # не должны попасть сюда. т.к. tasks[] не пустой
                task_to_derieve = task_with_max_f_value

            else:
                # перейти на Ш.7.
                pass

            # Ш.3. Осуществить ветвление задачи ЗЛП-k . Для этого выбрать нецелочисленную
            # координату и сформировать ограничения и 2 задачи
            simplex1, simplex2 = self.derieve(task_to_derieve, non_int_value)
            task1 = Task(2*k + 1, simplex1)
            task2 = Task(2*k + 2, simplex1)

            self.tasks.insert(task1.k, task1)
            self.tasks.insert(task2.k, task2)



            print('asda')

    # e.g.: x3 = 2.5
    def derieve(self, task: Task, non_int_value: dict) -> (Simplexx, Simplexx):
        # e.g: x3
        x_key = next(iter(non_int_value.keys()))

        # ЗЛП-1, x <= non_int
        # e.g: x3 <= 2
        x_less = int(next(iter(non_int_value.values())))
        simplex1 = self.create_simplex_with_additional_bound(task, x_key, x_less, Sign.PLUS)

        # ЗЛП-2, x >= non_int + 1
        # e.g: x3 >= 3 <-> -x3 <= -3 (то же самое, но со знаком минус)
        x_greater = int(next(iter(non_int_value.values()))) + 1
        simplex2 = self.create_simplex_with_additional_bound(task, x_key, x_greater, Sign.MINUS)

        return simplex1, simplex2

    # добавляем к текущей задаче новое ограничение "x_3 <= 5" в виде
    # дополнительой строки матрицы A: [0, 0, 0, 1], где 1 отвечает за x_3
    # и дополнительной строки столбца b: [5]
    def create_simplex_with_additional_bound(self, task: Task,
                                             new_bound_key: str,
                                             new_bound_value: int,
                                             new_bound_sign: Sign) -> Simplexx:
        # находим колонку, которая соответстует нецелому "x_i"
        # e.g: x_3 -> "x_" + "3" -> 3
        new_bound_column = int(new_bound_key.split('_')[1])

        matr_columns = task.simplex.matr.shape[0]
        additional_a_row = [0] * matr_columns
        if new_bound_sign is Sign.PLUS:
            # нумерация x_i идет с единицы
            additional_a_row[new_bound_column - 1] = 1
        else:
            # нумерация x_i идет с единицы
            additional_a_row[new_bound_column - 1] = -1
        additional_a_row = np.array(additional_a_row)

        # дпоисываем строку в низ матрицы А
        curr_a = task.simplex.matr
        new_a = np.vstack((curr_a, additional_a_row))

        additinal_b_row = None
        if new_bound_sign is Sign.PLUS:
            additinal_b_row = np.array([new_bound_value])
        else:
            additinal_b_row = np.array([-1 * new_bound_value])

        # дописываем значение ограничения в низ столбца b
        curr_b = task.simplex.b
        new_b = np.vstack((curr_b, additinal_b_row))

        return Simplexx(new_a, new_b, task.simplex.lambdas, task.simplex.condition)


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
