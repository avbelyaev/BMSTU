#!/usr/bin/env python3

BLOSUM62_FILENAME = 'blosum62.txt'
CONSTANT_GAP_PENALTY_SCORE = -3


class SubstMatrix:
    def __init__(self, filename: str):
        self.proteins = []
        self._m = {}
        self._load_matrix(filename)

    def _load_matrix(self, matrix_filename: str):
        with open(matrix_filename) as matrix_file:
            matrix = matrix_file.read()
        lines = matrix.strip().split('\n')

        header = lines.pop(0)
        columns = header.split()
        matrix = {}

        for row in lines:
            entries = row.split()
            row_name = entries.pop(0)
            matrix[row_name] = {}

            for column_name in columns:
                matrix[row_name][column_name] = entries.pop(0)

        self._m = matrix
        self.proteins = columns

    def lookup(self, a: str, b: str) -> int:
        return int(self._m[a][b])


def calc_scoring_matrix(s1: str, s2: str, subst: SubstMatrix, penalty: int) -> ([], int, int):
    rows = len(s1) + 1
    cols = len(s2) + 1

    matrix = [[0] * cols for _ in range(rows)]

    # первые строка и столбец = нули
    for i in range(rows):
        matrix[i][0] = 0

    for i in range(cols):
        matrix[0][i] = 0

    max_score = -1
    max_i = -1
    max_j = -1
    for i in range(1, rows):
        for j in range(1, cols):
            similarity_score = subst.lookup(s1[i - 1], s2[j - 1])

            match = matrix[i - 1][j - 1] + similarity_score
            delete = matrix[i - 1][j] + penalty
            insert = matrix[i][j - 1] + penalty

            matrix[i][j] = max(0, match, insert, delete)
            # запоминаем где было максимальное совпадение в таблице
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_i = i
                max_j = j

    return matrix, max_i, max_j


def smith_waterman(a: str, b: str, subst_matrix: SubstMatrix, penalty: int):
    m, max_i, max_j = calc_scoring_matrix(a, b, subst_matrix, penalty)

    res_a = ''
    res_b = ''

    # т.к. выравнивание _локальное_, начинаем сразу с позиции,
    # где было максмимальное совпадение. пытаемся расширить совпадение
    i = max_i
    j = max_j
    # идем назад (traceback) по максимальному пути.
    # таким образом, восстанавлиаем цепочку от конца к началу
    while i > 0 and j > 0 and m[i][j] > 0:

        score = m[i][j]
        score_diag = m[i - 1][j - 1]
        score_up = m[i][j - 1]
        score_left = m[i - 1][j]

        s = subst_matrix.lookup(a[i - 1], b[j - 1])

        # если идем по диагонали - все хорошо
        if score == (score_diag + s):
            res_a = a[i - 1] + res_a
            res_b = b[j - 1] + res_b
            i -= 1
            j -= 1

        # если идем влево, занчит в цепочке разрыв => пенальти
        elif score == (score_left + penalty):
            res_a = a[i - 1] + res_a
            res_b = '-' + res_b
            i -= 1

        # если идем вверх, занчит в цепочке разрыв => пенальти
        else:
            res_a = '-' + res_a
            res_b = b[j - 1] + res_b
            j -= 1

    return res_a, res_b


def main():
    # example from wiki
    s1 = 'TGTTACGG'
    s2 = 'GGTTGACTA'
    matrix = SubstMatrix(BLOSUM62_FILENAME)
    penalty = CONSTANT_GAP_PENALTY_SCORE

    res1, res2 = smith_waterman(s1, s2, matrix, penalty)
    assert res1 == 'GTT-AC'
    assert res2 == 'GTTGAC'

    print(f'{s1}\n{s2}')
    print('-------------')
    print(f'{res1}\n{res2}')


if __name__ == '__main__':
    main()
