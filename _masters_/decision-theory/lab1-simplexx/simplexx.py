import numpy as np

AXIS_OY = 1


LAMBDAS = [2, 8, 3]

def add_fictional_vars(a: np.matrix, signs: np.matrix, m_rows: int) -> np.matrix:
    column_width = 1
    for i in range(m_rows):
        if signs[i, 0] == '>=':
            fict_var_column = np.zeros((m_rows, column_width), dtype=int)
            fict_var_column[i][0] = -1
            a = np.c_[a, fict_var_column]

        elif signs[i, 0] == '<=':
            fict_var_column = np.zeros((m_rows, column_width), dtype=int)
            fict_var_column[i][0] = 1
            a = np.c_[a, fict_var_column]

        else:
            # equation? => do nothing
            pass
    return a


def step_1_solve():
    pass


def step_2_opt(matr: np.matrix, b: np.matrix):
    # supporting =
    pass


def create_simplex_table(matr: np.matrix, b: np.matrix):
    rows = matr.shape[0]
    cols = matr.shape[1]

    header_row = [' ', 's']
    i = 1
    while i <= cols:
        header_row.append(f'x_{i}')
        i += 1

    header_col = []
    j = 0
    while j < rows:
        # продолжаем счет переменных
        header_col.append(f'x_{i + j}')
        j += 1
    header_col = np.array(header_col)
    header_col = header_col.reshape(header_col.shape[0], 1)

    tbl = np.zeros((rows, cols))
    tbl = np.hstack((b, tbl))
    tbl = np.hstack((header_col, tbl))
    tbl = np.vstack((header_row, tbl))

    tbl = np.squeeze(np.asarray(tbl))
    print(tbl)


def simlex(matr: np.matrix, b: np.matrix):
    table = create_simplex_table(matr, b)

    matr = matr.astype(np.float).tolist()
    b = b.transpose().astype(np.float).tolist()[0]
    all_free_vars_are_positive = True
    for b_i in b:
        if b_i < 0:
            all_free_vars_are_positive = False
            break

    if all_free_vars_are_positive:
        result = step_2_opt(matr, b)
    else:
        step_1_solve()
        result = step_2_opt()

    return result


def main():
    matr = np.matrix([[2, 1, 1, '<=', 4],
                      [1, 2, 0, '<=', 6],
                      [0, 0.5, 1, '<=', 2]])

    b = matr[:, -1]
    signs = matr[:, -2]
    a = matr[:, :-2]
    m_rows = a.shape[0]
    n_cols = a.shape[1]

    print('--- Before ---')
    print(matr)

    print('--- fictional vars added ---')
    # fictionals = add_fictional_vars(a, signs, m_rows)
    # print(fictionals, b)

    result = simlex(a, b)
    print(result)


if __name__ == '__main__':
    main()
