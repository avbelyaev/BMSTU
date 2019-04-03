

INITIAL_POPULATION = [
    '11110000', '00001111'
]


class Individual:
    def __init__(self, str: str):
        self.str = str
        self.fitness =  0

    def set_fitness(self, fitness: int):
        self.fitness = fitness

    def __str__(self):
        return f'[{self.fitness}] {self.str}'

    def __repr__(self):
        return self.__str__()



def compute_fitness(individual: str):
    return individual.count("1")


def main():
    individuals = [Individual(s) for s in INITIAL_POPULATION]

    for individual in individuals:
        fit_score = compute_fitness(individual.str)
        individual.set_fitness(fit_score)


    while True:
        most_fit = selection()
        crossover()
        mutation()
        fitness = compute_fitness(initial_population)
    print('done')


if __name__ == '__main__':
    main()
