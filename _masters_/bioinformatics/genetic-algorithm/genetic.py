from random import randint

# генетический алгоритм.
# дан стеарата натрия
# пытаемся получить вещество максимально близкое к лаурил сульфат натрия по HLB
# функция фитнесса - HLB (гидро-липофильный баланс)


# @formatter:off
# стеарат натрия: C_18 H_35 Na O_2
SODIUM_STEARATE = [
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C(=O)[O-].[Na+]'
]
# лаурил сульфат натрия: C_12 H_25 SO_4 Na
SODIUM_LAURYL_SULFATE = [
    'C','C','C','C','C','C','C','C','C','C','C','C','O','C','C','OS(=O)([O-])=O.[Na+]'
]

# при мутации элементы меняются на одни из следующих
HEAD = ['OS(=O)([O-])=O.[Na+]', 'C(=O)[O-].[K+]', 'C(=O)[O-].[Na+]', 'C(=O)[O-]']
BODY = ['N', 'C', 'O']
TAIL = ['C']

# https://en.wikipedia.org/wiki/Hydrophilic-lipophilic_balance
RATING = {
    'OS(=O)([O-])=O.[Na+]': 38.7,   # SO4Na
    'C(=O)[O-].[K+]':       21.1,   # COOK
    'C(=O)[O-].[Na+]':      19.1,   # COONa
    'N':                    9.4,    # N
    'C(=O)[O-]':            2.1,    # COOH
    'O':                    1.9,    # OH
    'C':                    0.5     # C
}
# @formatter:on


INITIAL_POPULATION = [SODIUM_STEARATE]
CONVERGENCE_DELTA = 0.1


class Individual:
    """
    [C][C] ... [C][O][O][N] ... [SO4]
    хвост  ...      тело    ... голова
    """

    def __init__(self, molecules: list):
        self.elems = molecules
        self.hlb_score = 0
        self.compute_hlb_score()

    def compute_hlb_score(self):
        self.hlb_score = sum([RATING[elem] for elem in self.elems]) + 7

    def mutate(self):
        should_mutate_tail = lambda position: 0 == position
        should_mutate_head = lambda position: (len(self.elems) - 1) == position

        # pick range of molecules to mutate
        mutate_from, mutate_to = pick_random_range(self.elems)

        print(f'mutation: {self} -> ', end='')
        for i in range(mutate_from, mutate_to + 1):
            if should_mutate_head(i):
                self.__mutate_head()
            elif should_mutate_tail(i):
                self.__mutate_tail()
            else:
                self.__mutate_body_at(i)
        print(self)

    def __mutate_head(self):
        self.elems[len(self.elems) - 1] = pick_random_element(HEAD)

    def __mutate_tail(self):
        self.elems[0] = pick_random_element(TAIL)

    def __mutate_body_at(self, pos: int):
        self.elems[pos] = pick_random_element(BODY)

    def __str__(self):
        return f'[{self.hlb_score}] {"".join(self.elems)}'

    def __repr__(self):
        return self.__str__()


GOAL_GLB = Individual(SODIUM_LAURYL_SULFATE).hlb_score


def pick_random_element(elements: list):
    random = randint(0, len(elements) - 1)
    return elements[random]


def pick_random_range(elements: list):
    start = randint(0, len(elements) - 1)
    end = randint(0, len(elements) - 1)
    if start > end:
        return end, start
    return start, end


def selection(individs: list):
    print('selection')
    individs.sort(key=lambda ind: abs(GOAL_GLB - ind.hlb_score))


def crossover(parent1: Individual, parent2: Individual) -> (Individual, Individual):
    assert len(parent1.elems) == len(parent2.elems)
    crossover_point = randint(0, len(parent1.elems) - 1)

    child1_elems, child2_elems = parent1.elems.copy(), parent2.elems.copy()
    for i in range(crossover_point + 1):
        child1_elems[i], child2_elems[i] = child2_elems[i], child1_elems[i]

    child1, child2 = Individual(child1_elems), Individual(child2_elems)

    print(f'crossover')
    print(f'parnts: {parent1} + {parent2}\nchilds: {child1} , {child2}')
    return child1, child2


def mutation(individs: list):
    # pick range of individuals to mutate
    _, mutate_to = pick_random_range(individs)
    print(f'mutating [0..{mutate_to}]')
    for i in range(mutate_to):
        individs[i].mutate()


def main():
    individuals = [Individual(s) for s in INITIAL_POPULATION]

    mutated_stearate = Individual(SODIUM_STEARATE)
    mutated_stearate.mutate()
    individuals.append(mutated_stearate)

    i = 0
    while True:
        print(f'--- gen {i} ---')

        # count fitness
        [ind.compute_hlb_score() for ind in individuals]

        # sort individuals by fitness score
        selection(individuals)

        # check convergence
        most_fit = individuals[0]
        if abs(most_fit.hlb_score - GOAL_GLB) <= CONVERGENCE_DELTA:
            break

        # get most fit parents and do crossover
        parent1, parent2 = individuals[0], individuals[1]
        (child1, child2) = crossover(parent1, parent2)

        # add offspring to population
        individuals.extend([child1, child2])

        # mutate some individuals from population
        mutation(individuals)
        i += 1

    print(f'=== Done in generation #{i} ===')
    for individ in individuals[:10]:
        print(individ)


if __name__ == '__main__':
    main()
