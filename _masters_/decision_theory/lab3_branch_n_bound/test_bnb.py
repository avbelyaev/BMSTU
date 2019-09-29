#! /usr/bin/env python3.6
import unittest

import numpy as np

# relative to venv ? venv now at /BMSTU/venv
from _masters_.decision_theory.lab1_simplexx.simplexx import Condition, NoPivotalSolutionExists, \
    NoOptimalSolutionExists
from _masters_.decision_theory.lab3_branch_n_bound.bnb import BranchAndBound


class TestBranchAndBound(unittest.TestCase):
    def test_var_3(self):
        a = np.array([[2, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 1]])
        b = np.array([[4],
                      [6],
                      [2]])
        lambdas = np.array([[2, 8, 3]])

        # when
        print('===  Прямая ===')
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 25.5
        self.assertEqual(expected_f_value, solution['F'])

    def test_var_10(self):
        a = np.array([[4, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 1]])
        b = np.array([[4],
                      [3],
                      [2]])
        lambdas = np.array([[7, 5, 3]])

        # when
        print('===  Прямая ===')
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 13
        self.assertEqual(expected_f_value, solution['F'])

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
        solution = BranchAndBound(a, b, lambdas, Condition.MIN).run()

        # then
        expected_f_value = -3
        self.assertEqual(expected_f_value, solution['F'])

    # пример из методички. стр 40
    def test_example_2_from_book(self):
        a = np.array([[3, 1, -4, -1],
                      [-2, -4, -1, 1]])
        b = np.array([[-3],
                      [-3]])
        lambdas = np.array([[-4, -18, -30, -5]])

        # when
        print('===  Прямая ===')
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = -36
        self.assertEqual(expected_f_value, solution['F'])

    def test_unbounded_solution(self):
        # given
        a = np.array([[1, -1],
                      [1, 0]])
        b = np.array([[10],
                      [20]])
        lambdas = np.array([[1, 2]])

        print('===  Прямая ===')
        self.assertRaises(NoOptimalSolutionExists, BranchAndBound(a, b, lambdas, Condition.MAX).run)

    def test_no_allowed_solution(self):
        # given
        a = np.array([[2, 1],
                      [-3, -4]])
        b = np.array([[2],
                      [-12]])
        lambdas = np.array([[3, 2]])

        print('===  Прямая ===')
        self.assertRaises(NoPivotalSolutionExists, BranchAndBound(a, b, lambdas, Condition.MAX).run)


if __name__ == '__main__':
    unittest.main()
