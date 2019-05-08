BLOSUM62_FILENAME = 'blosum62.txt'

CONSTANT_GAP_PENALTY_SCORE = 5


class SubstMatrix:
    def __init__(self, filename: str):
        self._load_matrix(filename)

    def _load_matrix(self, matrix_filename):
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

        self.m = matrix

    def lookup(self, a: str, b: str) -> int:
        return int(self.m[a][b])


def calc_scoring_matrix(s1: str, s2: str, subst_matrix: SubstMatrix, penalty: int) -> []:
    height = len(s1) + 1
    width = len(s2) + 1

    # alloc matrix
    matrix = [[0] * width for _ in range(height)]

    for i in range(height):
        matrix[i][0] = 0

    for i in range(width):
        matrix[0][i] = 0

    max_score = 0
    max_i = 0
    max_j = 0
    for i in range(1, height):
        for j in range(1, width):
            similarity_score = subst_matrix.lookup(s1[i - 1], s2[j - 1])

            match = matrix[i - 1][j - 1] + similarity_score
            delete = matrix[i - 1][j] + penalty
            insert = matrix[i][j - 1] + penalty

            matrix[i][j] = max(0, match, insert, delete)
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_i = i
                max_j = j

    return matrix, max_i, max_j


def smith_waterman(a: str, b: str, subst_matrix: SubstMatrix, penalty: int):
    m, max_i, max_j = calc_scoring_matrix(a, b, subst_matrix, penalty)

    res_a = ''
    res_b = ''

    i = len(a)
    while i > max_i:
        res_a = a[i - 1] + res_a
        res_b = '-' + res_b
        i -= 1

    j = len(b)
    while j > max_j:
        res_a = '-' + res_a
        res_b = b[j - 1] + res_b
        j -= 1

    i = max_i
    j = max_j

    while i > 0 and j > 0 and m[i][j] > 0:

        score = m[i][j]
        score_diag = m[i - 1][j - 1]
        score_up = m[i][j - 1]
        score_left = m[i - 1][j]

        s = subst_matrix.lookup(a[i - 1], b[j - 1])

        if score == (score_diag + s):
            res_a = a[i - 1] + res_a
            res_b = b[j - 1] + res_b
            i -= 1
            j -= 1

        elif score == (score_left + penalty):
            res_a = a[i - 1] + res_a
            res_b = '-' + res_b
            i -= 1

        else:
            res_a = '-' + res_a
            res_b = b[j - 1] + res_b
            j -= 1

    while i > 0:
        res_a = a[i - 1] + res_a
        res_b = '-' + res_b
        i -= 1

    while j > 0:
        res_a = '-' + res_a
        res_b = b[j - 1] + res_b
        j -= 1

    return res_a, res_b


def main():
    s1 = 'TGTTACGG'
    s2 = 'GGTTGACTA'
    matrix = SubstMatrix(BLOSUM62_FILENAME)
    penalty = CONSTANT_GAP_PENALTY_SCORE

    res1, res2 = smith_waterman(s1, s2, matrix, penalty)

    print(res1)
    print(res2)


if __name__ == '__main__':
    main()
