#! /usr/bin/env python3.6
import unittest

# relative to venv ? venv now at /BMSTU/venv
import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx


class TestSimplexxMethods(unittest.TestCase):
    def test_simplex(self):
        # given
        a = np.array([[1, -2],
                      [-2, 1],
                      [1, 1]])
        b = np.array([[2],
                      [-2],
                      [5]])
        signs = np.array([['='],
                          ['='],
                          ['=']])
        lambdas = np.array([1, -1])  # TODO ищем минимум, хотя в задании указан максимум

        # when
        s = Simplexx(a, b, lambdas, signs)
        solutions = s.run()

        # then
        expected_best_solution = ({'x_3': 0, 'x_5': 0, 'x_4': 5.0, 'x_1': 4.0, 'x_2': 1.0}, {'F': -3.0})
        self.assertEqual(expected_best_solution, solutions[-1])

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
        expected_tbl = np.array([[1, 1/2, -3/2],
                                 [1, -1/2, -1/2],
                                 [4, 1/2, 3/2],
                                 [-1, 1/2, -1/2]])
        expected_header_top = ['b', 'x_4', 'x_2']
        expected_header_left = ['x_3', 'x_1', 'x_5']

        # and
        self.assertEqual(expected_header_left, simplex.header_left)
        self.assertEqual(expected_header_top, simplex.header_top)
        np.testing.assert_array_equal(expected_tbl, simplex.tbl)

    def test_change_basis_2(self):
        # given
        simplex = Simplexx(None, None, None, None)
        simplex.tbl = np.array([[1, 1/2, -3/2],
                                [1, -1/2, -1/2],
                                [4, 1/2, 3/2],
                                [-1, 1/2, -1/2]], dtype='float64')
        simplex.header_top = ['b', 'x_4', 'x_2']
        simplex.header_left = ['x_3', 'x_1', 'x_5']

        # when
        determining_row = 0
        determining_col = 1
        simplex.change_basis(determining_row, determining_col)

        # then
        expected_tbl = np.array([[2, 2, -3],
                                 [2, 1, -2],
                                 [3, -1, 3],
                                 [-2, -1, 1]])
        expected_header_top = ['b', 'x_3', 'x_2']
        expected_header_left = ['x_4', 'x_1', 'x_5']

        # and
        self.assertEqual(expected_header_left, simplex.header_left)
        self.assertEqual(expected_header_top, simplex.header_top)
        np.testing.assert_array_equal(expected_tbl, simplex.tbl)

    def test_change_basis_3(self):
        # given
        simplex = Simplexx(None, None, None, None)
        simplex.tbl = np.array([[2, 2, -3],
                                [2, 1, -2],
                                [3, -1, 3],
                                [-2, -1, 1]], dtype='float64')
        simplex.header_top = ['b', 'x_3', 'x_2']
        simplex.header_left = ['x_4', 'x_1', 'x_5']

        # when
        determining_row = 2
        determining_col = 2
        simplex.change_basis(determining_row, determining_col)

        # then
        expected_tbl = np.array([[5, 1, 1],
                                 [4, 1/3, 2/3],
                                 [1, -1/3, 1/3],
                                 [-3, -2/3, -1/3]])
        expected_header_top = ['b', 'x_3', 'x_5']
        expected_header_left = ['x_4', 'x_1', 'x_2']

        # and
        self.assertEqual(expected_header_left, simplex.header_left)
        self.assertEqual(expected_header_top, simplex.header_top)
        # comparing floats
        np.testing.assert_allclose(expected_tbl, simplex.tbl)


if __name__ == '__main__':
    unittest.main()
