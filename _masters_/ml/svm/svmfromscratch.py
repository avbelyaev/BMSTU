from math import pi, sin, cos, sqrt
from copy import deepcopy
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from _masters_.ml.square_contours.contours import Point

IMG_NAME = 'dots.png'

DATASET_SIZE = 5
CLAZZ_1 = 'a'
CLAZZ_2 = 'b'


def dist_between(p1: Point, p2: Point) -> float:
    return sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)


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


def xs(points: list) -> list:
    return [p.x for p in points]


def ys(points: list) -> list:
    return [p.y for p in points]


class Vector:
    def __init__(self, pt_from: Point, pt_to: Point):
        self.a = deepcopy(pt_from)
        self.b = deepcopy(pt_to)

    def rotate(self, angle_degree: float):
        radians = angle_degree / 180 * pi

        # translate to 0
        self.b.x -= self.a.x
        self.b.y -= self.a.y

        # roatate
        print(f'rota: {datetime.now()}')
        xnew = self.b.x * cos(radians) - self.b.y * sin(radians)
        ynew = self.b.x * sin(radians) + self.b.y * cos(radians)

        # translate back
        self.b.x = xnew + self.a.x
        self.b.y = ynew + self.a.y

    @property
    def len(self):
        return dist_between(self.a, self.b)

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
                    'dist': dist_between(pt1, pt2)
                })

        distances.sort(key=lambda it: it['dist'])
        closest_1, closest_2 = distances[0][CLAZZ_1], distances[0][CLAZZ_2]
        middle = find_middle_point(closest_1, closest_2)
        middle_vect = Vector(closest_1, closest_2)

        transcluster_12 = Vector(closest_1, closest_2)
        transcluster_21 = Vector(closest_2, closest_1)


        step = 10
        angle = 0
        while angle < 360:
            angle += step
            print(angle)

            print(f'calc: {datetime.now()}')
            transcluster_12.rotate(step)
            transcluster_21.rotate(step)

            print(f' clr: {datetime.now()}')
            plt.clf()
            plt.cla()
            # draw clusters
            print(f'draw: {datetime.now()}')
            plt.scatter(xs(data[CLAZZ_1]), ys(data[CLAZZ_1]))
            plt.scatter(xs(data[CLAZZ_2]), ys(data[CLAZZ_2]))

            # draw vectors
            plt.scatter([middle.x], [middle.y])
            plt.plot(middle_vect.xs, middle_vect.ys, 'g--')
            plt.plot(transcluster_12.xs, transcluster_12.ys, 'r--')
            plt.plot(transcluster_21.xs, transcluster_21.ys, 'r--')

            print(f'rndr: {datetime.now()}')
            # plt.savefig(IMG_NAME)
            # time.sleep(0.3)
            plt.pause(0.05)
            plt.draw()
            print(f'done: {datetime.now()}')

        # plt.show()


def main():
    center1 = (30, 60)
    center2 = (80, 20)
    distance = 20

    x1 = np.random.uniform(low=center1[0], high=center1[0] + distance, size=(DATASET_SIZE,))
    y1 = np.random.normal(loc=center1[1], scale=distance, size=(DATASET_SIZE,))

    x2 = np.random.uniform(center2[0], center2[0] + distance, size=(DATASET_SIZE,))
    y2 = np.random.normal(center2[1], distance, size=(DATASET_SIZE,))

    # plt.scatter(x1, y1)
    # plt.scatter(x2, y2)
    # plt.show()

    clazz1_points = []
    for i in range(DATASET_SIZE):
        x = x1[i]
        y = y1[i]
        clazz1_points.append(Point(x, y))

    clazz2_points = []
    for i in range(DATASET_SIZE):
        x = x2[i]
        y = y2[i]
        clazz2_points.append(Point(x, y))

    marked_up_data = {
        CLAZZ_1: clazz1_points,
        CLAZZ_2: clazz2_points
    }

    svm = SVMFromScratch()
    svm.fit(marked_up_data)
    # svm.visualize()


if __name__ == '__main__':
    main()
