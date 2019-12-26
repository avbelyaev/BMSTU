import random

import numpy as np

N_AGENTS = 3
CONVERGENCE_EPSILON = 1e-6


class AntagonisticGame:
    def __init__(self, matrix: np.ndarray):
        self.trust_matrix = matrix

    def converge(self, opinions: np.ndarray) -> (np.ndarray, int):
        def has_converged(opinion_vector: np.ndarray) -> bool:
            # сойдется, когда число мнений с отклонением EPS будет равно числу агентов
            return len(np.where(opinion_vector < CONVERGENCE_EPSILON)[0]) == N_AGENTS

        i = 0
        curr = np.dot(self.trust_matrix, opinions)
        while True:
            prev = curr.copy()
            # текущий вектор мнений
            curr = np.dot(self.trust_matrix, prev)
            # вектор разностей между текущими мнениями и мнениями на прошлом шаге
            delta = np.abs(curr - prev)
            if has_converged(delta):
                return np.round(curr, decimals=1), i
            i += 1


def generate_matrix(n: int) -> np.ndarray:
    m = np.zeros((n, n))
    for i in range(n - 1):
        for j in range(n - 1):
            m[i, j] = round(random.random() / 5, 3)

    for i in range(n - 1):
        m[i, n - 1] = round(1 - np.sum(m[i, :n - 1]), 3)

    for j in range(n):
        m[n - 1, j] = round(1 - np.sum(m[:n - 1, j]), 3)

    while np.sign(np.min(m)) < 0:
        for i in range(n - 1):
            for j in range(n - 1):
                m[i, j] = round(random.random() / (n / 2), 3)

        for i in range(n - 1):
            m[i, n - 1] = round(1 - np.sum(m[i, :n - 1]), 3)

        for j in range(n):
            m[n - 1, j] = round(1 - np.sum(m[:n - 1, j]), 3)
    return m


def influence(opinions: np.ndarray) -> np.ndarray:
    # формируем агентов и их мнения для 1го игрока
    agents_1 = np.random.randint(low=0, high=N_AGENTS // 2,
                                 size=random.randint(1, N_AGENTS // 2))
    agents_1_opinion = random.randint(0, 100)
    # то же самое для 2го игрока
    agents_2 = np.random.randint(low=N_AGENTS // 2, high=N_AGENTS,
                                 size=random.randint(1, N_AGENTS // 2))
    agents_2_opinion = random.randint(-100, 0)

    print(f'Агенты 1: {agents_1}')
    print(f'Агенты 2: {agents_2}')

    influenced_opinions = []
    for i in range(N_AGENTS):
        if i in agents_1:
            influenced_opinions.append(agents_1_opinion)
        elif i in agents_2:
            influenced_opinions.append(agents_2_opinion)
        else:
            # агент не зависит ни от какого игрока =>его мнение не меняется
            influenced_opinions.append(opinions[i])
    return np.array(influenced_opinions)


def main():
    m = generate_matrix(N_AGENTS)
    m = np.array([[0.403, 0.541, 0.056],
                  [0.561, 0.199, 0.24],
                  [0.036, 0.26, 0.704]])
    print(m)
    print('Матрица доверия:')

    i = 0
    while i < 3:
        m = np.dot(m ,m )
        i += 1

    print(m)

    game = AntagonisticGame(m)

    independent_opinions = np.random.randint(low=1, high=20, size=N_AGENTS)
    independent_opinions = np.array([19, 6, 7])
    resulting_opinions, iterations = game.converge(independent_opinions)
    print('Без влияния')
    print(f'Независимые мнения: {independent_opinions}')
    print(f'Результирующие мнения [{iterations}]: {resulting_opinions}')

    print('C влиянием')
    influenced_opinions = influence(independent_opinions)
    resulting_opinions, iterations = game.converge(influenced_opinions)
    print(f'Сформированные мнения: {influenced_opinions}')
    print(f'Результирующие мнения [{iterations}]:{resulting_opinions}')


if __name__ == '__main__':
    main()
