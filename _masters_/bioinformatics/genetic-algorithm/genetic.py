from random import randint

INITIAL_POPULATION = [
    '11110000', '00001111'
]


class Individual:
    def __init__(self, string: str):
        self.str = string
        self.fitness = 0

    def compute_fitness(self):
        self.fitness = self.str.count('1')

    def mutate(self):
        mutate_from = randint(0, len(self.str) - 1)
        mutate_to = randint(0, len(self.str) - 1)
        if mutate_from > mutate_to:
            mutate_from, mutate_to = mutate_to, mutate_from

        mutated_str = self.str
        i = mutate_from
        while i < mutate_to:
            mutated_char = '1' if '0' == self.str[i] else '0'
            mutated_str = setCharAt(mutated_str, mutated_char, i)
            i += 1

        print(f'mutation: {self.str} -> {mutated_str}')
        self.str = mutated_str

    def __str__(self):
        return f'{self.str}'

    def __repr__(self):
        return self.__str__()


def selection(individs: list):
    print('selection')
    individs.sort(key=lambda ind: ind.fitness, reverse=True)


def setCharAt(string: str, new_char: str, pos: int):
    chars = list(string)
    chars[pos] = new_char
    return ''.join(chars)


def crossover(parent1: Individual, parent2: Individual) -> (Individual, Individual):
    parent1_str = parent1.str
    parent2_str = parent2.str

    assert len(parent1_str) == len(parent2_str)
    crossover_point = randint(0, len(parent1_str) - 1)

    child1_str = parent1_str
    child2_str = parent2_str
    i = 0
    while i < crossover_point:
        child1_str = setCharAt(child1_str, parent2_str[i], i)
        child2_str = setCharAt(child2_str, parent1_str[i], i)
        i += 1

    print(f'crossover[{crossover_point}]')
    print(f'parnts: {parent1_str} + {parent2_str}\nchilds: {child2_str} , {child1_str}')
    return Individual(child1_str), Individual(child2_str)


def mutation(individs: list):
    to_be_mutated = randint(1, len(individs) - 1)
    i = 0
    while i < to_be_mutated:
        individs[i].mutate()
        i += 1


def has_converged(individs: list) -> bool:
    return individs[0].fitness >= 6


def main():
    individuals = [Individual(s) for s in INITIAL_POPULATION]
    [ind.compute_fitness() for ind in individuals]

    i = 0
    while True:
        print(f'--- gen {i} ---')

        # sort individuals by fitness score
        selection(individuals)

        # check algorithm convergence
        if has_converged(individuals):
            break

        # get most fit and do crossover
        parent1, parent2 = individuals[0], individuals[1]
        (child1, child2) = crossover(parent1, parent2)

        # add offspring to population
        individuals.extend([child1, child2])

        # mutate some individuals from population
        mutation(individuals)

        # count fitness
        [ind.compute_fitness() for ind in individuals]
        i += 1

    print('=============')
    print(individuals)


if __name__ == '__main__':
    main()
