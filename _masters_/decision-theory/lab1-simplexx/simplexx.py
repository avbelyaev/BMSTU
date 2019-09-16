import numpy as np


class Simplexx:
    def __init__(self, a: np.ndarray, b: np.ndarray, lambdas: np.ndarray, signs: np.ndarray):
        self.matr = a
        self.b = b
        self.lambdas = lambdas
        self.signs = signs.reshape(signs.shape[0], 1)
        self.tbl = None
        self.header_top = []
        self.header_left = []

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
        b = self._to_column(self.b)
        tbl = np.hstack((b, tbl))

        # добавляем строку лямбд внизу
        pos = 0
        additional_zero_elem = 0
        lambdas = np.insert(self.lambdas, pos, additional_zero_elem, axis=0)
        tbl = np.vstack((tbl, lambdas))
        return tbl

    def free_members_are_positive(self) -> bool:
        return True

    def run(self) -> (dict, float):
        self.tbl = self.create_simplex_table()

        if self.free_members_are_positive():
            # опорное решение найдено
            vars = self.get_variables_mapping()
            value = self.target_func()
            return (vars, value)

    def target_func(self) -> float:
        # так как все свобоные переменные = 0,
        # то ответ лежит в первой клетке подвала таблицы
        rows = self.tbl.shape[0]
        return self.tbl[rows - 1, 0]

    def get_variables_mapping(self) -> dict:
        res = dict()
        for x_i in self.header_top:
            if 'b' == x_i:
                continue
            else:
                res[x_i] = 0
        j = 0
        for x_j in self.header_left:
            res[x_j] = self.tbl[j, 0]
            j += 1
        return res

    def _to_column(self, xs: np.ndarray) -> np.ndarray:
        return xs.reshape(xs.shape[0], 1)


def main():
    a = np.array([[2, 1, 1],
                  [1, 2, 0],
                  [0, 0.5, 1]])
    b = np.array([4, 6, 2])
    lambdas = np.array([2, 8, 3])
    signs = np.array(['<=', '<=', '<='])

    s = Simplexx(a, b, lambdas, signs)
    vars, value = s.run()
    print(vars)
    print(f'F = {value}')


if __name__ == '__main__':
    main()
