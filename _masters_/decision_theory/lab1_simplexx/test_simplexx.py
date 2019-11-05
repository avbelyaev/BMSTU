#! /usr/bin/env python3.6
import unittest

# relative to venv ? venv now at /BMSTU/venv
import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, NoOptimalSolutionExists, NoPivotalSolutionExists, \
    Condition


class TestSimplexxMethods(unittest.TestCase):
    def test_branch_1(self):
        a = np.array([[6, -1],
                      [2, 5],
                      [1, 0]])
        b = np.array([[12],
                      [20],
                      [2]])
        lambdas = np.array([[12, -1]])

        # when
        solution = Simplexx(a, b, lambdas, Condition.MAX).run()

        # then
        actual_f_value = solution['F']
        self.assertEqual(24, actual_f_value)

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

    # пример из методички, стр 31
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

    # пример из методички. стр 40
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

    def test_variable_mapping(self):
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
        primary_solution = Simplexx(a, b, lambdas, Condition.MAX).run()

        # then
        expected_f_value = 13
        self.assertEqual(expected_f_value, primary_solution['F'])

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


if __name__ == '__main__':
    unittest.main()
