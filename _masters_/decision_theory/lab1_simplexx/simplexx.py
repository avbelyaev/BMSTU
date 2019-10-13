from enum import Enum
from typing import Optional

import numpy as np
from copy import deepcopy


class NoPivotalSolutionExists(Exception):
    def __init__(self, msg: str = 'Допустимого(опорного) решения не существует'):
        super().__init__(msg)


class NoOptimalSolutionExists(Exception):
    def __init__(self):
        super().__init__('Фукнция не ограничена! Оптимального решения не существует')


class Condition(Enum):
    MIN = 0
    MAX = 1


class Simplexx:
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, condition: Condition):
        """
        in form of AX<=B
        :param a: A
        :param b: B
        :param lambdas: aka C
        :param condition: min or max
        """
        self.matr = a
        self.b = b
        self.lambdas = lambdas
        self.tbl = None
        self.header_top = []
        self.header_left = []
        self.solutions = []
        self.condition = condition

    def create_simplex_table(self) -> np.ndarray:
        row_num = self.matr.shape[0]
        col_num = self.matr.shape[1]

        self.header_top = ['b']
        i = 1
        while i <= col_num:
            self.header_top.append(f'x_{i}')
            i += 1

        self.header_left = []
        j = 0
        while j < row_num:
            # продолжаем счет переменных
            self.header_left.append(f'x_{i + j}')
            j += 1

        tbl = self.matr

        # добавляем колонку b'шек слева
        tbl = np.hstack((self.b, tbl))

        # алгоритм ниже расчитан на поиск MAX
        # но ответ домножаем на -1
        if self.condition is Condition.MAX:
            pass
        # если надо найти MIN, то
        # - строку лямбд (aka С) домножаем на -1
        # - ответ (значение F в конце решения) оставляем как есть
        elif self.condition is Condition.MIN:
            self.lambdas = -1 * self.lambdas
        else:
            raise ValueError('Неверное условие')

        # добавляем 0 в начало строки лямбд
        pos = 0
        additional_zero_elem = 0
        lambdas = np.insert(self.lambdas[0], pos, additional_zero_elem, axis=0)

        # добавляем строку лямбд внизу
        tbl = np.vstack((tbl, lambdas))
        tbl = tbl.astype(dtype='float64')
        return tbl

    # проверим, что в столбце свободных членов все эл-ты положительные
    # иначе вернем строку с отрицательным элементом
    def find_negative_free_var_row(self) -> Optional[int]:
        i = 0
        while i < self._get_rows() - 1:  # пропускаем подвал таблицы
            if self.tbl[i, 0] < 0:
                return i
            i += 1
        return None

    # В строке ищем первый отрицательный элемент
    def find_first_negative_col(self, row: int) -> Optional[int]:
        #  пропускаем 0й столбец со свободными переменными
        i = 1
        while i < self._get_cols():
            if self.tbl[row, i] < 0:
                return i
            i += 1
        return None

    # поиск резрешающего столбца
    def find_determining_column(self) -> Optional[int]:
        lambdas_row = self.tbl[self._get_rows() - 1:]
        labdas_row_len = self._get_cols()
        i = 1  # пропускаем столбец свободных сленов
        while i < labdas_row_len:
            if float(lambdas_row[0, i]) > 0:
                # ищем первый положительный элемент
                # возвращаем разрешающий столбец
                return i
            i += 1
        return None

    # поиск резрешающей строки
    def find_determining_row(self, determining_col: int) -> Optional[int]:
        # Найдем минимальное положительное отношение элемента свободных членов
        # si0 к соответствующем эле- менту в разрешающем столбце
        min_relation = 999999
        determining_row = None
        i = 0
        while i < self._get_rows() - 1:
            if 0 == self.tbl[i, determining_col]:
                i += 1
                continue

            curr_relation = self.tbl[i, 0] / self.tbl[i, determining_col]
            if 0 < curr_relation < min_relation:
                min_relation = curr_relation
                determining_row = i
            i += 1
        # разрешающая строка
        return determining_row

    def change_basis(self, r: int, k: int):
        print(f'Замена базиса: {self.header_left[r]} <-> {self.header_top[k]}, row: {r}, col: {k}')
        # r - разр. строка
        # k - разр. столбец
        s_rk = self.tbl[r, k]
        self.tbl[r, k] = 1 / s_rk

        original_table = deepcopy(self.tbl)

        # меняем разрешающую строку, кроме разрешающего элемента
        col = 0
        while col < self._get_cols():
            if col == k:
                col += 1
                continue
            else:
                self.tbl[r, col] = original_table[r, col] / s_rk
            col += 1

        # обновляем разрешающий столбец
        row = 0
        while row < self._get_rows():
            if row == r:
                row += 1
                continue
            else:
                self.tbl[row, k] = -1 * original_table[row, k] / s_rk
            row += 1

        # обновляем все остальное
        row = 0
        while row < self._get_rows():
            col = 0
            while col < self._get_cols():
                if row == r or col == k:
                    col += 1
                    continue
                else:
                    t = original_table[row, k] * original_table[r, col] / s_rk
                    self.tbl[row, col] = original_table[row, col] - t
                col += 1
            row += 1

        # меняем иксы в колонке и столбце
        tmp = self.header_top[k]
        self.header_top[k] = self.header_left[r]
        self.header_left[r] = tmp

    def run(self) -> dict:
        self.tbl = self.create_simplex_table()

        print('Поиск опорного решения')
        while True:
            negative_free_var_row = self.find_negative_free_var_row()
            if negative_free_var_row is None:
                break

            determining_col = self.find_first_negative_col(negative_free_var_row)
            if determining_col is None:
                raise NoPivotalSolutionExists()

            # Найдем минимальное положительное отношение элмента свободных членов si0
            # к соответствующем эле- менту в разрешающем столбце
            determining_row = self.find_determining_row(determining_col)
            if determining_row is None:
                break

            self.change_basis(determining_row, determining_col)

        # Так как все элементы столбца si0 неотрицательны, имеем опорное решение
        self.add_solution()
        print('Опорное решение:')
        print(self.best_solution)

        print('Поиск оптимального решения')
        while True:
            determining_col = self.find_determining_column()
            if determining_col is None:
                break

            determining_row = self.find_determining_row(determining_col)
            if determining_row is None:
                raise NoOptimalSolutionExists()

            self.change_basis(determining_row, determining_col)
            self.add_solution()

            print('Более оптимальное решение:')
            print(self.best_solution)

        return self.best_solution

    @property
    def best_solution(self) -> dict:
        # положительность x'ов - это одно из условий задачи
        def all_xs_are_positive(sol: dict) -> bool:
            for x in sol.keys():
                if x != 'F' and sol[x] < 0:
                    return False
            return True

        # последнее решение, при котором все x положительные.
        i = len(self.solutions) - 1
        while i >= 0:
            # идем от последнего к первому решению
            if all_xs_are_positive(self.solutions[i]):
                return self.solutions[i]
            i -= 1
        raise NoPivotalSolutionExists('Нет решения с положительными "x"')

    def add_solution(self):
        variables = self.get_variables_mapping()
        f_value = self.get_target_func()
        self.solutions.append({**variables, **{'F': f_value}})

    def get_target_func(self) -> float:
        # так как все свобоные переменные = 0,
        # то ответ лежит в первой клетке подвала таблицы
        rows = self.tbl.shape[0]
        f_value = round(self.tbl[rows - 1, 0], ndigits=2)

        #  см заметики при составлении таблицы
        if self.condition is Condition.MIN:
            return f_value
        elif self.condition is Condition.MAX:
            return -1 * f_value
        else:
            raise ValueError('Неверное условие')

    def get_variables_mapping(self) -> dict:
        res = dict()
        for x_i in self.header_top:
            if 'b' == x_i:
                continue
            else:
                res[x_i] = 0
        j = 0
        for x_j in self.header_left:
            res[x_j] = round(self.tbl[j, 0], ndigits=2)
            j += 1
        return res

    def _get_rows(self) -> int:
        return self.tbl.shape[0]

    def _get_cols(self) -> int:
        return self.tbl.shape[1]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.matr}'


if __name__ == '__main__':
    print('см. тесты')
    pass
