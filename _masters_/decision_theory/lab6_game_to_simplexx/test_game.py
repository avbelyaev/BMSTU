import unittest

import numpy as np

from _masters_.decision_theory.lab6_game_to_simplexx.game import MatrixGame


class TestMatrixGame(unittest.TestCase):
    # стр 78. методички
    def test_sample_from_book(self):
        m = np.array([[1, 3, 9, 6],
                      [2, 6, 2, 3],
                      [7, 2, 6, 5]])
        game = MatrixGame(m)

        # when
        strat_a = game.play_as_a()
        # then
        expected_w_value = 0.246  # 14/57
        self.assertEqual(expected_w_value, strat_a['F'])

        # when
        strat_b = game.play_as_b()
        # then
        self.assertEqual(expected_w_value, strat_b['F'])


    def test_var_3(self):
        m = np.array([[1, 11, 6, 15, 10],
                      [3, 3, 15, 10, 7],
                      [4, 7, 16, 0, 10],
                      [16, 18, 7, 7, 5]])
        game = MatrixGame(m)

        # when
        expected_w_value = 0.126
        strategy_a = game.play_as_a()
        # then
        self.assertEqual(expected_w_value, strategy_a['F'])

        # when
        strategy_b = game.play_as_b()
        # then
        self.assertEqual(expected_w_value, strategy_b['F'])



if __name__ == '__main__':
    unittest.main()
