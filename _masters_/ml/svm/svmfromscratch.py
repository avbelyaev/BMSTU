from math import pi, sin, cos, sqrt
from copy import deepcopy
import matplotlib.pyplot as plt

from _masters_.ml.square_contours.contours import Point

DATASET_SIZE = 5
CLAZZ_1 = 'a'
CLAZZ_2 = 'b'

FULL_DEGREE = 360


class Vector:
    def __init__(self, pt_from: Point, pt_to: Point):
        self.a = deepcopy(pt_from)
        self.b = deepcopy(pt_to)

    def rotate(self, angle_degree: float):
        radians = angle_degree / 180 * pi

        # translate to 0
        self.b.x -= self.a.x
        self.b.y -= self.a.y

        # rotate
        xnew = self.b.x * cos(radians) - self.b.y * sin(radians)
        ynew = self.b.x * sin(radians) + self.b.y * cos(radians)

        # translate back
        self.b.x = xnew + self.a.x
        self.b.y = ynew + self.a.y

    @property
    def linear_function(self) -> dict:
        """in form of: a*x + b*y + c = 0"""
        # (y1 – y2)x + (x2 – x1)y + (x1y2 – x2y1) = 0
        return {
            'a': self.a.y - self.b.y,
            'b': self.b.x - self.a.x,
            'c': (self.a.x * self.b.y) - (self.b.x * self.a.y)
        }

    @staticmethod
    def dist_between_parallel(v1: 'Vector', v2: 'Vector') -> float:
        linear1 = v1.linear_function
        linear2 = v2.linear_function
        # assert a1 == a2, b1 == b2
        denominator = sqrt(linear1['a'] ** 2 + linear1['b'] ** 2)
        return abs(linear1['c'] + linear2['c']) // denominator

    @property
    def len(self) -> float:
        return dist_between_points(self.a, self.b)

    @property
    def xs(self) -> list:
        return [self.a.x, self.b.x]

    @property
    def ys(self) -> list:
        return [self.a.y, self.b.y]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.a} -> {self.b}'


class SVMFromScratch:
    def __init__(self):
        self.data = None

    def fit(self, data: dict):
        distances = []
        for pt1 in data[CLAZZ_1]:
            for pt2 in data[CLAZZ_2]:
                distances.append({
                    CLAZZ_1: pt1,
                    CLAZZ_2: pt2,
                    'dist': dist_between_points(pt1, pt2)
                })

        distances.sort(key=lambda it: it['dist'])
        closest_1, closest_2 = distances[0][CLAZZ_1], distances[0][CLAZZ_2]
        middle = find_middle_point(closest_1, closest_2)
        middle_vect = Vector(closest_1, closest_2)

        svm1 = Vector(closest_1, closest_2)
        svm2 = Vector(closest_2, closest_1)

        most_distant_support_vectors = {
            'vect1': None,
            'vect2': None,
            'dist': -9999
        }

        step = 10
        angle = 0
        while angle < FULL_DEGREE:
            angle += step
            print(f'angle: {angle}')

            svm1.rotate(step)
            svm2.rotate(step)
            curr_dist = Vector.dist_between_parallel(svm1, svm2)
            print(f' dist: {curr_dist}')

            c1_correct = SVMFromScratch.is_correctly_classified(svm1, data[CLAZZ_1])
            c2_correct = SVMFromScratch.is_correctly_classified(svm2, data[CLAZZ_2])
            print(f'c1: {c1_correct}, c2: {c2_correct}')

            if curr_dist > most_distant_support_vectors['dist'] \
                    and SVMFromScratch.is_correctly_classified(svm1, data[CLAZZ_1]) \
                    and SVMFromScratch.is_correctly_classified(svm2, data[CLAZZ_2]):
                most_distant_support_vectors = {
                    'vect1': deepcopy(svm1),
                    'vect2': deepcopy(svm2),
                    'dist': curr_dist
                }

            # clear canvas
            plt.clf()
            plt.cla()

            # draw clusters and mid point
            plt.scatter(xs(data[CLAZZ_1]), ys(data[CLAZZ_1]))
            plt.scatter(xs(data[CLAZZ_2]), ys(data[CLAZZ_2]))
            plt.scatter([middle.x], [middle.y])

            # draw support vectors
            plt.plot(middle_vect.xs, middle_vect.ys, 'g--')
            plt.plot(svm1.xs, svm1.ys, 'r--')
            plt.plot(svm2.xs, svm2.ys, 'r--')

            plt.pause(2)
            plt.draw()

        print(most_distant_support_vectors)
        plt.show()


    @staticmethod
    def is_correctly_classified(vect: Vector, points: list) -> bool:
        signs = []
        for p in points:
            # https://stackoverflow.com/a/3461533/4504720
            sign = (vect.b.x - vect.a.x) * (p.y - vect.a.y) - (vect.b.y - vect.a.y) * (p.x - vect.a.x)
            signs.append(sign)
        # проверяем что все точки по одну строну от прямой (т.е. либо все > 0, либо все < 0)
        # та точка, через которую проходит опорный вектор дает 0 => ослабляем до >= 0 (или <= 0)
        return all(sign >= 0 for sign in signs) or all(sign <= 0 for sign in signs)


def find_middle_point(p1: Point, p2: Point) -> Point:
    if p1.x < p2.x:
        mid_x = p1.x + (p2.x - p1.x) // 2
    else:
        mid_x = p2.x + (p1.x - p2.x) // 2

    if p1.y < p2.y:
        mid_y = p1.y + (p2.y - p1.y) // 2
    else:
        mid_y = p2.y + (p1.y - p2.y) // 2

    return Point(mid_x, mid_y)


def dist_between_points(p1: Point, p2: Point) -> float:
    return sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)


def xs(points: list) -> list:
    return [p.x for p in points]


def ys(points: list) -> list:
    return [p.y for p in points]


def main():
    center1 = (30, 60)
    center2 = (80, 20)
    distance = 20

    import numpy as np
    x1 = np.random.uniform(low=center1[0], high=center1[0] + distance, size=(DATASET_SIZE,))
    y1 = np.random.normal(loc=center1[1], scale=distance, size=(DATASET_SIZE,))

    x2 = np.random.uniform(center2[0], center2[0] + distance, size=(DATASET_SIZE,))
    y2 = np.random.normal(center2[1], distance, size=(DATASET_SIZE,))

    # convert from numpy into normal points
    marked_up_data = {
        CLAZZ_1: [],
        CLAZZ_2: []
    }

    for i in range(DATASET_SIZE):
        x = x1[i]
        y = y1[i]
        marked_up_data[CLAZZ_1].append(Point(x, y))

    for i in range(DATASET_SIZE):
        x = x2[i]
        y = y2[i]
        marked_up_data[CLAZZ_2].append(Point(x, y))

    svm = SVMFromScratch()
    svm.fit(marked_up_data)


if __name__ == '__main__':
    main()

    # v = Vector(Point(2,1), Point(5,7))
    # print(v.linear_function)
    # linear1 = {'a': -31.247440394778728, 'b': 9.66091408069402, 'c': 964.2556114111776}
    # linear2 = {'a': 31.24744039477874, 'b': -9.660914080694027, 'c': -2033.9914035109805}
    # # assert a1 == a2, b1 == b2
    # denominator = sqrt(linear1['a'] ** 2 + linear1['b'] ** 2)
    # print(abs(linear1['c'] + linear2['c']) / denominator)
