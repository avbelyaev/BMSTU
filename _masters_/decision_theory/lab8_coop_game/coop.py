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

COALITION_PROFIT_AFTER_CHANGE = {
    frozenset([]): 0,
    frozenset([1]): 3,
    frozenset([2]): 1,
    frozenset([3]): 2,
    frozenset([4]): 4,
    frozenset([1, 2]): 4,
    frozenset([1, 3]): 5,
    frozenset([1, 4]): 8,
    frozenset([2, 3]): 3,
    frozenset([2, 4]): 5,
    frozenset([3, 4]): 6,
    frozenset([1, 2, 3]): 7,
    frozenset([1, 2, 4]): 10,
    frozenset([1, 3, 4]): 11,
    frozenset([2, 3, 4]): 8,
    frozenset([1, 2, 3, 4]): 13
}
PLAYERS = [1, 2, 3, 4]


class CooperativeGame:
    def __init__(self, characteristic_discrete_func: dict, players: list):
        self.profit = characteristic_discrete_func
        self.coalitions = all_combinations(players)
        self.players = players
        self.num_players = len(players)
        self.total_coalition_profit = self.profit[frozenset(self.players)]

    def solve_sheply_for_player(self, player_id: int) -> float:
        s = 0
        for coalition in self.coalitions:
            # кол-во способов сформировать коалицию "до прихода игрока"
            coalition_before = len(coalition) - 1
            # разность между общим кол-вом игроков и коалицией
            players_difference = self.num_players - len(coalition)
            # полезность коалиции с игроком
            profit_with = self.profit[coalition]
            # полезность коалиции без игрока
            coalition_no_player = coalition - {player_id}
            profit_without = self.profit[coalition_no_player]

            s += fact(coalition_before) \
                 * fact(players_difference) \
                 * (profit_with - profit_without)

        s /= fact(self.num_players)
        return round(s, ndigits=3)

    # супераддитивность
    def is_super_additive(self) -> bool:
        for i in self.coalitions:
            for j in self.coalitions:
                if {} == i.intersection(j):
                    profit_ij = self.profit[i.union(j)]
                    profit_i = self.profit[i]
                    profit_j = self.profit[j]
                    # сумма выгод по отдельности не больше выгоды при объединении
                    if (profit_i + profit_j) > profit_ij:
                        print(f'выгода {i}: {profit_i}')
                        print(f'выгода {j}: {profit_j}')
                        print(f'выгода ({i} & {j}): {profit_ij}')
                        print(f'{profit_i} + {profit_j} > {profit_ij}')
                        return False

        individual_profits = 0
        for player_id in self.players:
            individual_profits += self.profit[frozenset([player_id])]
        # выгода тотальной коалиции больше отдельных выгод
        print('Выгода тотальной коалиции:', self.total_coalition_profit)
        if self.total_coalition_profit < individual_profits:
            return False
        return True

    # выпуклость
    def is_convex(self) -> bool:
        for i in self.coalitions:
            for j in self.coalitions:
                union_ij = i.union(j)
                intersection_ij = i.intersection(j)
                profit_i = self.profit[i]
                profit_j = self.profit[j]
                if (self.profit[union_ij] + self.profit[intersection_ij]) < (profit_i + profit_j):
                    print(f'Коалиция 1: {i}, коалиция 2: {j}')
                    print(f'Объединение: {union_ij}, пересечение: {intersection_ij}')
                    print(f'({self.profit[union_ij]} + {self.profit[intersection_ij]}) '
                          f'< ({profit_i} + {profit_j})')
                    return False
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
    print('Игра суперадитивная:', game.is_super_additive())
    print('Игра выпуклая:', game.is_convex())

    print('\n--- меняем условие ---\n')
    game = CooperativeGame(COALITION_PROFIT_AFTER_CHANGE, PLAYERS)
    print('Игра суперадитивная:', game.is_super_additive())
    print('Игра выпуклая:', game.is_convex())

    shares = []
    print('Индивидуальная рационализация:')
    for player_id in PLAYERS:
        share = game.solve_sheply_for_player(player_id)
        print(f'Player {player_id}:\n'
              f'  Выгода: {share} > индивидуальной выгоды: {COALITION_PROFIT[frozenset([player_id])]}')
        shares.append(share)

    print(f'Групповая рационализация:\n'
          f'Сумма долей: {sum(shares)}, Выгода тотальной коалиции: {game.total_coalition_profit}')


if __name__ == '__main__':
    main()
