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

    def create_simplex_table(self):
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

        b = self._to_column(self.b)
        tbl = np.hstack((b, tbl))

        pos = 0
        additional_zero_elem = 0
        lambdas = np.insert(self.lambdas, pos, additional_zero_elem, axis=0)
        tbl = np.vstack((tbl, lambdas))
        self.tbl = tbl

    def run(self):
        self.create_simplex_table()
        print('asda')

    def _to_column(self, xs: np.ndarray) -> np.ndarray:
        return xs.reshape(xs.shape[0], 1)

    # def run(matr: np.matrix, b: np.matrix):
    #     table = create_simplex_table(matr, b)
    #
    #     matr = matr.astype(np.float).tolist()
    #     b = b.transpose().astype(np.float).tolist()[0]
    #     all_free_vars_are_positive = True
    #     for b_i in b:
    #         if b_i < 0:
    #             all_free_vars_are_positive = False
    #             break
    #
    #     if all_free_vars_are_positive:
    #         result = step_2_opt(matr, b)
    #     else:
    #         step_1_solve()
    #         result = step_2_opt()
    #
    #     return result


def main():
    a = np.array([[2, 1, 1],
                  [1, 2, 0],
                  [0, 0.5, 1]])
    b = np.array([4, 6, 2])
    lambdas = np.array([2, 8, 3])
    signs = np.array(['<=', '<=', '<='])

    s = Simplexx(a, b, lambdas, signs)
    result = s.run()
    print(result)


if __name__ == '__main__':
    main()
