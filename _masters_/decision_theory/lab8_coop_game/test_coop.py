import unittest

from _masters_.decision_theory.lab8_coop_game.coop import CooperativeGame


class TestCooperativeGame(unittest.TestCase):
    # вариант из методички, стр 129
    def test_var_from_book(self):
        profits = {
            frozenset([]): 0,
            frozenset([1]): 1,
            frozenset([2]): 1,
            frozenset([3]): 1,
            frozenset([1, 2]): 3,
            frozenset([1, 3]): 3,
            frozenset([2, 3]): 3,
            frozenset([1, 2, 3]): 4,
        }
        players = [1, 2, 3]

        game = CooperativeGame(profits, players)

        self.assertTrue(game.is_super_additive())
        self.assertTrue(game.is_convex())

        for player_id in players:
            share = game.solve_sheply_for_player(player_id)
            print('Player', player_id, share)

            self.assertAlmostEqual(1.333, share)  # 4/3


if __name__ == '__main__':
    unittest.main()
