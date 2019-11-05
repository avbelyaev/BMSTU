from typing import Optional

import numpy as np

from _masters_.decision_theory.lab1_simplexx.simplexx import Simplexx, Condition, NoPivotalSolutionExists, \
    NoOptimalSolutionExists
from _masters_.decision_theory.lab2_simplexx_duality.duality import DualSimplexx


class MatrixGame:
    def __init__(self, strategy_matrix: np.ndarray):
        # приводим к форме Ax<=b
        self.matrix = np.negative(strategy_matrix.transpose())
        cols = self.matrix.shape[1]
        rows = self.matrix.shape[0]
        self.b = MatrixGame.create_b_column(rows)
        self.lambdas = np.array([[1] * cols])
        self.solution_a, self.solution_b = None, None

    def play_as_a(self) -> Optional[dict]:
        print('Стратегия игрока A')
        try:
            self.solution_a = Simplexx(self.matrix, self.b, self.lambdas, Condition.MIN).run()
            return self.solution_a
        except (NoPivotalSolutionExists, NoOptimalSolutionExists) as e:
            print(e)
            return None

    def play_as_b(self) -> Optional[dict]:
        print('Стратегия игрока B')
        try:
            self.solution_b = DualSimplexx(self.matrix, self.b, self.lambdas, Condition.MIN).run()
            return self.solution_b
        except (NoPivotalSolutionExists, NoOptimalSolutionExists) as e:
            print(e)
            return None

    @staticmethod
    def create_b_column(height: int):
        b = []
        for i in range(height):
            b.append([-1])
        return np.array(b)

    @staticmethod
    def as_player_strategy(solution: dict) -> dict:
        w_value = solution['F']
        strategy = dict()
        for x in solution.keys():
            if x != 'F':
                strategy[x] = round(solution[x] / w_value, ndigits=3)
        return strategy


def main():
    m = np.array([[1, 11, 6, 15, 10],
                  [3, 3, 15, 10, 7],
                  [4, 7, 16, 0, 10],
                  [16, 18, 7, 7, 5]])
    game = MatrixGame(m)

    solution_a = game.play_as_a()
    print(solution_a)

    solution_b = game.play_as_b()
    print(solution_b)

    print('\n\n---')
    print('Оптимальная смешанная стратегия игрока A')
    print(MatrixGame.as_player_strategy(solution_a))
    print('Оптимальная смешанная стратегия игрока B')
    print(MatrixGame.as_player_strategy(solution_b))


if __name__ == '__main__':
    main()
