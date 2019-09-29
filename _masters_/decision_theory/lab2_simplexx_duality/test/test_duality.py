#! /usr/bin/env python3.6
import unittest

import numpy as np

# relative to venv ? venv now at /BMSTU/venv
from _masters_.decision_theory.lab1_simplexx.simplexx import Condition, NoPivotalSolutionExists, \
    NoOptimalSolutionExists, Simplexx
from _masters_.decision_theory.lab2_simplexx_duality.duality import DualSimplexx


class TestDualSimplexxMethods(unittest.TestCase):
    def test_var_3_MAX_MIN(self):
        a = np.array([[2, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 1]])
        b = np.array([[4],
                      [6],
                      [2]])
        lambdas = np.array([[2, 8, 3]])

        # when
        print('===  Прямая ===')
        primary_solution = Simplexx(a, b, lambdas, Condition.MAX).run()

        # when
        print('\n\n\n===  Двойственная ===')
        dual_solution = DualSimplexx(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 25.5
        self.assertEqual(expected_f_value, primary_solution['F'])
        self.assertEqual(expected_f_value, dual_solution['F'])

    def test_example_1_from_book(self):
        a = np.array([[1, -2],
                      [-2, 1],
                      [1, 1]])
        b = np.array([[2],
                      [-2],
                      [5]])
        lambdas = np.array([[-1, 1]])

        # when
        print('===  Прямая ===')
        primary_solution = Simplexx(a, b, lambdas, Condition.MIN).run()

        # when
        print('\n\n\n===  Двойственная ===')
        dual_solution = DualSimplexx(a, b, lambdas, Condition.MIN).run()

        # then
        expected_f_value = -3
        self.assertEqual(expected_f_value, primary_solution['F'])
        self.assertEqual(expected_f_value, dual_solution['F'])

    def test_example_2_from_book(self):
        a = np.array([[3, 1, -4, -1],
                      [-2, -4, -1, 1]])
        b = np.array([[-3],
                      [-3]])
        lambdas = np.array([[-4, -18, -30, -5]])

        # when
        print('===  Прямая ===')
        primary_solution = Simplexx(a, b, lambdas, Condition.MAX).run()

        # when
        print('\n\n\n===  Двойственная ===')
        dual_solution = DualSimplexx(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = -36
        self.assertEqual(expected_f_value, primary_solution['F'])
        self.assertEqual(expected_f_value, dual_solution['F'])

    # неограниченное решение
    def test_unbounded_solution(self):
        # given
        a = np.array([[1, -1],
                      [1, 0]])
        b = np.array([[10],
                      [20]])
        lambdas = np.array([[1, 2]])

        # expect
        self.assertRaises(NoOptimalSolutionExists, DualSimplexx(a, b, lambdas, Condition.MAX).run)

    # нет допустимого решения
    def test_no_allowed_solution(self):
        # given
        a = np.array([[2, 1],
                      [-3, -4]])
        b = np.array([[2],
                      [-12]])
        lambdas = np.array([[3, 2]])

        # expect
        self.assertRaises(NoPivotalSolutionExists, DualSimplexx(a, b, lambdas, Condition.MAX).run)


if __name__ == '__main__':
    unittest.main()
