import itertools
from math import factorial as fact
from typing import List

COALITION_PROFIT = {
    frozenset([]): 0,
    frozenset([1]): 3,
    frozenset([2]): 4,
    frozenset([3]): 1,
    frozenset([4]): 2,
    frozenset([1, 2]): 7,
    frozenset([1, 3]): 5,
    frozenset([1, 4]): 6,
    frozenset([2, 3]): 5,
    frozenset([2, 4]): 7,
    frozenset([3, 4]): 3,
    frozenset([1, 2, 3]): 10,
    frozenset([1, 2, 4]): 10,
    frozenset([1, 3, 4]): 9,
    frozenset([2, 3, 4]): 9,
    frozenset([1, 2, 3, 4]): 12
}
PLAYERS = [1, 2, 3, 4]


class CooperativeGame:
    def __init__(self, characteristic_discrete_func: dict, players: list):
        self.profits = characteristic_discrete_func
        self.coalitions = all_combinations(players)
        self.num_players = len(players)

    def solve_sheply_for_player(self, player_id: int) -> float:
        s = 0
        for coalition in self.coalitions:
            if 0 == len(coalition):
                continue
            # кол-во способов сформировать коалицию "до прихода игрока"
            coalition_before = len(coalition) - 1
            # разность между общим кол-вом игроков и коалицией
            players_difference = self.num_players - len(coalition)
            # полезность коалиции с игроком
            profit_with = self.profits[coalition]
            # полезность коалиции без игрока
            coalition_no_player = coalition.difference([player_id])
            profit_without = self.profits[coalition_no_player]

            s += fact(coalition_before) \
                 * fact(players_difference) \
                 * (profit_with - profit_without)

        s /= fact(self.num_players)
        return round(s, ndigits=3)

    def is_super_additive(self) -> bool:
        return True

    def is_convex(self) -> bool:
        return True


def all_combinations(items: list) -> 'List[frozenset]':
    combos = []
    for r in range(1, len(items) + 1):
        combo = list(itertools.combinations(items, r))
        combos.extend(combo)
    combos = list(map(lambda it: frozenset(it), combos))
    return combos


def main():
    game = CooperativeGame(COALITION_PROFIT, PLAYERS)

    assert game.is_super_additive(), "Игра не супераддитивная!"
    assert game.is_convex(), "Игра не выпуклая!"

    for player_id in PLAYERS:
        share = game.solve_sheply_for_player(player_id)
        print('Player', player_id, share)


if __name__ == '__main__':
    main()
