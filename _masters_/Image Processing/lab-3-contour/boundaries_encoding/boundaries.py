import abc
import math
from enum import Enum
import matplotlib.pyplot as plt
from skimage import io


WHITE_CONTOUR_MASK = [200, 200, 200]
RED = 0
GREEN = 1
BLUE = 2

# FILENAME = 'circle.png'
# FILENAME = 'star.png'
FILENAME = 'letter.png'
IMG = io.imread(FILENAME).tolist()

STEP = 37


class Direction(Enum):
    CLOCKWISE = 1
    COUNTER_CLOCKWISE = -1


class Point:
    def __init__(self, x: int, y: int):
        self.x = x
        self.y = y
        self.is_contour = is_contour(IMG[self.x][self.y])

    def __eq__(self, other):
        return other is not None and \
               self.x == other.x and self.y == other.y

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'[{self.x}:{self.y}] cont:{self.is_contour}'


class Vector:
    def __init__(self, pt_from: Point, pt_to: Point):
        self.pt_from = pt_from
        self.pt_to = pt_to

    # радиус-вектор представляется точкой
    def to_radius_vector(self) -> Point:
        delta_x = self.pt_to.x - self.pt_from.x
        delta_y = self.pt_to.y - self.pt_from.y
        return Point(delta_x, delta_y)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.pt_from} -> {self.pt_to}'


class Boundary:
    def encode(self) -> str:
        raise NotImplementedError


class ThreeAttrBoundary(Boundary):
    """
    Элемент границы, закодированной по трем признакам:
    (длина вектора, угол, направление поворота к следующему вектору)
    """

    def __init__(self, v_curr: Vector, v_next: Vector):
        self.len = vector_len(v_curr)
        self.angle = self._normalize_angle(v_curr, v_next)
        self.direction = Direction.COUNTER_CLOCKWISE if self.angle > 0 \
            else Direction.CLOCKWISE

    # угол [0 360) преобразовать в угол [0 180) со знаком
    # например, 300" == -60"
    def _normalize_angle(self, v1: Vector, v2: Vector) -> float:
        big_angle = directionwise_angle(v1, v2)
        return big_angle if big_angle <= 180 else big_angle - 360

    def encode(self) -> str:
        return f'[{self.len:.2f} \tangle:{self.angle:.2f} {self.direction}]'


# угол между радиус-векторами в диапазоне [0, 360)
def directionwise_angle(v1: Vector, v2: Vector) -> float:
    p1, p2 = v1.to_radius_vector(), v2.to_radius_vector()

    v1_theta = math.atan2(p1.y, p1.x)
    v2_theta = math.atan2(p2.y, p2.x)
    r = (v2_theta - v1_theta) * (180.0 / math.pi)
    if r < 0:
        r += 360.0
    return r


def vector_len(v: Vector) -> float:
    return math.sqrt((v.pt_to.x - v.pt_from.x) ** 2 + (v.pt_to.y - v.pt_from.y) ** 2)


# pixel = [R, G, B]
def is_contour(pixel) -> bool:
    def is_of_color(color: list) -> bool:
        return pixel[RED] >= color[RED] and \
               pixel[GREEN] >= color[GREEN] and \
               pixel[BLUE] >= color[BLUE]

    return is_of_color(WHITE_CONTOUR_MASK)


def extract_contour_pts(img) -> list:
    # находим первую попавшуюся точку контура
    def find_start_pt(img) -> Point:
        rows = len(img)
        cols = len(img[0])
        for i in range(rows):
            for j in range(cols):
                if is_contour(img[i][j]):
                    return Point(i, j)

    start_pt = find_start_pt(img)
    curr_pt = start_pt
    prev_pt = None

    contour = []
    while True:
        contour.append(curr_pt)
        x, y = curr_pt.x, curr_pt.y

        nearest_points = [
            Point(x, y + 1),
            Point(x + 1, y + 1),
            Point(x + 1, y),
            Point(x + 1, y - 1),
            Point(x, y - 1),
            Point(x - 1, y - 1),
            Point(x - 1, y),
            Point(x - 1, y + 1)
        ]
        next_pt = None
        # рассматриваем 8 соседей
        # среди соседей будут несколько контурных точек
        for p in nearest_points:
            # нам надо найти еще НЕ посещенную точку
            if p.is_contour and p not in contour:
                next_pt = p
                break

        # не нашли непосещенной точки => выходим
        if next_pt is None:
            break

        curr_pt = next_pt

    print(f'contour size: {len(contour)}')
    print(f'start pt: {start_pt}')
    print(f'end   pt: {contour[len(contour) - 1]}')
    return contour


def approximate_contour(contour: list, step) -> list:
    vectors = []
    ctr_len = len(contour)
    for i in range(0, ctr_len - step, step):
        p0 = contour[i]
        p1 = contour[i + step]
        vectors.append(Vector(p0, p1))

    # соединим последнюю точку последнего вектора с первой точкой первого вектора
    # vectors[-1].pt_to.x = vectors[0].pt_from.x
    # vectors[-1].pt_to.y = vectors[0].pt_from.y
    return vectors


def three_attr_encoding(vectors: list) -> list:
    boundaries = []
    vec_len = len(vectors)
    for i in range(vec_len):
        v1 = vectors[i]
        v2 = vectors[(i + 1) % vec_len]
        boundaries.append(ThreeAttrBoundary(v1, v2))

    return boundaries


def main():
    points = extract_contour_pts(IMG)
    vectors = approximate_contour(points, STEP)

    for v in vectors:
        plt.plot([v.pt_from.x, v.pt_to.x], [v.pt_from.y, v.pt_to.y])
    plt.show()

    enc1 = three_attr_encoding(vectors)
    print('Three attribute encoding')
    [print(e.encode()) for e in enc1]


if __name__ == '__main__':
    main()
