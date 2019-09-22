#! /usr/bin/env python3.6
import unittest

# relative to venv ? venv now at /BMSTU/venv
import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx


class TestSimplexxMethods(unittest.TestCase):
    def test_change_basis_1(self):
        # given
        simplex = Simplexx(None, None, None, None)
        simplex.tbl = np.array([[2, 1, -2],
                                [-2, -2, 1],
                                [5, 1, 1],
                                [0, 1, -1]], dtype='float64')
        simplex.header_top = ['b', 'x_1', 'x_2']
        simplex.header_left = ['x_3', 'x_4', 'x_5']

        # when
        determining_row = 1
        determining_col = 1
        simplex.change_basis(determining_row, determining_col)

        # then
        expected_tbl = np.array([[1, 0.5, -1.5],
                                 [1, -0.5, -0.5],
                                 [4, 0.5, 1.5],
                                 [-1, 0.5, -0.5]])
        expected_header_top = ['b', 'x_4', 'x_2']
        expected_header_left = ['x_3', 'x_1', 'x_5']

        # and
        self.assertEqual(expected_header_left, simplex.header_left)
        self.assertEqual(expected_header_top, simplex.header_top)
        np.testing.assert_array_equal(expected_tbl, simplex.tbl)


if __name__ == '__main__':
    unittest.main()
