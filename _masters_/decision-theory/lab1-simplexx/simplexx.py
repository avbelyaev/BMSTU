import numpy as np

AXIS_OY = 1


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


def step_2_opt(matr: list, b: list):
    pass


def simlex(matr: np.matrix, b: np.matrix):
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
    c = [2, 8, 3]
    b = matr[:, -1]
    signs = matr[:, -2]
    a = matr[:, :-2]
    m_rows = a.shape[0]
    n_cols = a.shape[1]

    print('--- Before ---')
    print(matr)

    print('--- fictional vars added ---')
    fictionals = add_fictional_vars(a, signs, m_rows)
    print(fictionals, b)

    result = simlex(fictionals, b)
    print(result)


if __name__ == '__main__':
    main()
