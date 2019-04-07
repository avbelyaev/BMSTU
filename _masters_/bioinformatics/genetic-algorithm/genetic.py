from random import randint

# генетический алгоритм.
# пытаемся получить лаурил сульфат натрия из стеарата натрия
# функция фитнесса - GLB


# @formatter:off
# стеарат натрия: C_18 H_35 Na O_2
SODIUM_STEARATE = [
    'C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C(=O)[O-].[Na+]'
]
# лаурил сульфат натрия: C_12 H_25 SO_4 Na
SODIUM_LAURYL_SULFATE = [
    'C','C','C','C','C','C','C','C','C','C','C','C','O','C','C','OS(=O)([O-])=O.[Na+]'
]

HEAD = ['OS(=O)([O-])=O.[Na+]', 'C(=O)[O-].[K+]', 'C(=O)[O-].[Na+]', 'C(=O)[O-]']
BODY = ['N', 'C', 'O']
TAIL = ['C']

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
CONVERGENCE_DELTA = 0.01


class Individual:
    def __init__(self, molecules: list):
        self.elems = molecules
        self.head = molecules[len(molecules) - 1]
        self.tail = molecules[0]
        self.glb_score = 0

    def compute_glb_score(self):
        self.glb_score = sum([RATING[elem] for elem in self.elems])

    def mutate(self):
        should_mutate_tail = lambda position: 0 == position
        should_mutate_head = lambda position: (len(self.elems) - 1) == position

        # pick range of molecules to mutate
        mutate_from, mutate_to = pick_random_range(self.elems)

        print(f'mutation: {self} -> ', end='')
        for i in range(mutate_from, mutate_to + 1):
            if should_mutate_head(i):
                self.__set_head(pick_random_element(HEAD))
            elif should_mutate_tail(i):
                self.__set_tail(pick_random_element(TAIL))
            else:
                self.__set_body_at(pick_random_element(BODY), i)
        print(self)

    def __set_head(self, new_head_element: str):
        self.head = new_head_element
        self.elems[len(self.elems) - 1] = new_head_element

    def __set_tail(self, new_tail_element: str):
        self.tail = new_tail_element
        self.elems[0] = new_tail_element

    def __set_body_at(self, new_body_element: str, pos: int):
        self.elems[pos] = new_body_element

    def __str__(self):
        elems_as_str = ''.join(self.elems)
        return f'[{self.glb_score}] {elems_as_str}'

    def __repr__(self):
        return self.__str__()


# TODO fixme
sulfate = Individual(SODIUM_LAURYL_SULFATE)
sulfate.compute_glb_score()
GOAL_GLB = sulfate.glb_score


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
    individs.sort(key=lambda ind: abs(GOAL_GLB - ind.glb_score))


def crossover(parent1: Individual, parent2: Individual) -> (Individual, Individual):
    assert len(parent1.elems) == len(parent2.elems)
    crossover_point = randint(0, len(parent1.elems) - 1)

    child1_elems, child2_elems = parent1.elems.copy(), parent2.elems.copy()
    for i in range(crossover_point + 1):
        child1_elems[i], child2_elems[i] = child2_elems[i], child1_elems[i]

    child1, child2 = Individual(child1_elems), Individual(child2_elems)

    print(f'crossover[{crossover_point}]')
    print(f'parnts: {parent1} + {parent2}\nchilds: {child1} , {child2}')
    return child1, child2


def mutation(individs: list):
    # pick range of individuals to mutate
    mutate_from, mutate_to = pick_random_range(individs)
    print(f'mutating [{mutate_from}..{mutate_to}]')
    for i in range(mutate_from, mutate_to):
        individs[i].mutate()


def has_converged(individs: list) -> bool:
    most_fit = individs[0]
    return abs(most_fit.glb_score - GOAL_GLB) < CONVERGENCE_DELTA


def main():
    individuals = [Individual(s) for s in INITIAL_POPULATION]

    stearate = Individual(SODIUM_STEARATE)
    stearate.mutate()
    individuals.append(stearate)

    i = 0
    while True:
        print(f'--- gen {i} ---')

        # count fitness
        [ind.compute_glb_score() for ind in individuals]

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
        i += 1

    print('=== Done ===')
    for individ in individuals[:10]:
        print(individ)


if __name__ == '__main__':
    main()
