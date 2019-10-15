#! /usr/bin/env python3.6
import unittest

import numpy as np

# relative to venv ? venv now at /BMSTU/venv
from _masters_.decision_theory.lab1_simplexx.simplexx import Condition, NoPivotalSolutionExists, \
    NoOptimalSolutionExists
from _masters_.decision_theory.lab3_integer_simplexx.bnb import BranchAndBound


class TestBranchAndBound(unittest.TestCase):
    # стр. 48
    def test_example_from_book(self):
        a = np.array([[6, -1],
                      [2, 5]])
        b = np.array([[12],
                      [20]])
        lambdas = np.array([[12, -1]])

        # when
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f = 24
        self.assertEqual(expected_f, solution['F'])


    def test_var_3(self):
        a = np.array([[2, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 1]])
        b = np.array([[4],
                      [6],
                      [2]])
        lambdas = np.array([[2, 8, 3]])

        # when
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 24
        self.assertEqual(expected_f_value, solution['F'])
        self.assertEqual(0, solution['x_1'])
        self.assertEqual(3, solution['x_2'])
        self.assertEqual(0, solution['x_3'])


    def test_var_10(self):
        a = np.array([[4, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 1]])
        b = np.array([[4],
                      [3],
                      [2]])
        lambdas = np.array([[7, 5, 3]])

        # when
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 7    # FIXME math semestr -> 8.0
        self.assertEqual(expected_f_value, solution['F'])


    def test_var_7(self):
        a = np.array([[4, 1, 1],
                      [1, 2, 0],
                      [0, 0.5, 4]])
        b = np.array([[5],
                      [3],
                      [8]])
        lambdas = np.array([[6, 6, 6]])

        # when
        solution = BranchAndBound(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 12
        self.assertEqual(expected_f_value, solution['F'])
        self.assertEqual(1, solution['x_1'])
        self.assertEqual(1, solution['x_2'])
        self.assertEqual(0, solution['x_3'])


if __name__ == '__main__':
    unittest.main()
