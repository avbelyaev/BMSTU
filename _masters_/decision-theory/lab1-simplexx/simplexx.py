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


def main():
    matr = np.matrix([[1, 1, '=', 10],
                      [-2, 3, '<=', -5],
                      [7, -4, '>=', 6]])
    b = matr[:, -1]
    signs = matr[:, -2]
    a = matr[:, :-2]
    m_rows = a.shape[0]
    n_cols = a.shape[1]

    print('--- Before ---')
    print(matr)

    print('--- fictional vars added ---')
    reshaped = add_fictional_vars(a, signs, m_rows)
    print(reshaped)


if __name__ == '__main__':
    main()
