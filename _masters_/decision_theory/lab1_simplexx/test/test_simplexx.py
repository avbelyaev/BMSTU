#! /usr/bin/env python3.6
import unittest

# relative to venv ? venv now at /BMSTU/venv
import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, NoOptimalSolutionExists, NoPivotalSolutionExists, \
    Condition


class TestSimplexxMethods(unittest.TestCase):
    def test_var_3_MAX_MIN(self):
        a = np.array([[2, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 1]])
        b = np.array([[4],
                      [6],
                      [2]])
        lambdas = np.array([[2, 8, 3]])

        # when
        print('=== MAX ===')
        solution = Simplexx(a, b, lambdas, Condition.MAX).run()

        # then
        actual_f_value = solution['F']
        self.assertEqual(25.5, actual_f_value)

        # when
        print('=== MIN ===')
        solution = Simplexx(a, b, lambdas, Condition.MIN).run()

        # then
        actual_f_value = solution['F']
        self.assertEqual(0, actual_f_value)

    def test_example_1_from_book(self):
        a = np.array([[1, -2],
                      [-2, 1],
                      [1, 1]])
        b = np.array([[2],
                      [-2],
                      [5]])
        lambdas = np.array([[-1, 1]])

        # when
        solution = Simplexx(a, b, lambdas, Condition.MIN).run()

        # then
        actual_f_value = solution['F']
        self.assertEqual(-3, actual_f_value)

    def test_example_2_from_book(self):
        a = np.array([[3, 1, -4, -1],
                      [-2, -4, -1, 1]])
        b = np.array([[-3],
                      [-3]])
        lambdas = np.array([[-4, -18, -30, -5]])

        # when
        solution = Simplexx(a, b, lambdas, Condition.MAX).run()

        # then
        actual_f_value = solution['F']
        self.assertEqual(-36, actual_f_value)

    def test_simplex_1(self):
        # given
        a = np.array([[1, -2],
                      [-2, 1],
                      [1, 1]])
        b = np.array([[2],
                      [-2],
                      [5]])
        lambdas = np.array([[1, -1]])

        # when
        solutions = Simplexx(a, b, lambdas, Condition.MAX).run()

        # then
        expected_best_solution = ({
            'x_1': 4.0,
            'x_2': 1.0,
            'x_3': 0,
            'x_4': 5.0,
            'x_5': 0,
            'F': 3.0
        })
        self.assertEqual(expected_best_solution, solutions)

    # неограниченное решение
    def test_unbounded_solution(self):
        # given
        a = np.array([[1, -1],
                      [1, 0]])
        b = np.array([[10],
                      [20]])
        lambdas = np.array([[1, 2]])

        # expect
        self.assertRaises(NoOptimalSolutionExists, Simplexx(a, b, lambdas, Condition.MAX).run)

    # нет допустимого решения
    def test_no_allowed_solution(self):
        # given
        a = np.array([[2, 1],
                      [-3, -4]])
        b = np.array([[2],
                      [-12]])
        lambdas = np.array([[3, 2]])

        # expect
        self.assertRaises(NoPivotalSolutionExists, Simplexx(a, b, lambdas, Condition.MAX).run)

    def test_change_basis_1(self):
        # given
        simplex = Simplexx(None, None, None, Condition.MAX)
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
        simplex = Simplexx(None, None, None, Condition.MAX)
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
        simplex = Simplexx(None, None, None, Condition.MAX)
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
