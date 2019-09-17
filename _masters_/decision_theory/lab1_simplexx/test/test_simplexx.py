#! /usr/bin/env python3.6
import unittest

# relative to venv ? venv now at /BMSTU/venv
import numpy as np
import numpy.testing

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx


class TestSimplexxMethods(unittest.TestCase):
    def test_change_basis(self):
        a = np.array([[1, -2],
                      [-2, 1],
                      [1, 1]])
        b = np.array([[2],
                      [-2],
                      [5]])
        lambdas = np.array([-1, 1])
        simplex = Simplexx(a, b, lambdas, None)

        simplex.tbl = np.array([[-2, 1, -2],
                                [-2, -2, 1],
                                [5, 1, 1],
                                [0, -1, 1]])
        swap_row = 1
        swap_col = 2
        simplex.change_basis(swap_row, swap_col)

        actual_tbl = simplex.tbl
        actual_header_top = simplex.header_top
        actual_header_left = simplex.header_left

        expected_tbl = np.array([[1, 0.5, -1.5],
                                 [1, -0.5, -0.5],
                                 [4, 0.5, 1.5],
                                 -1, 0.5, -0.5])
        expected_header_top = ['b', 'x_4', 'x_2']
        expected_header_left = ['x_3', 'x_1', 'x_5']

        np.testing.assert_array_equal(actual_tbl, expected_tbl)
        self.assertEqual(actual_header_left, expected_header_left)
        self.assertEqual(actual_header_top, expected_header_top)

        if __name__ == '__main__':
            unittest.main()
