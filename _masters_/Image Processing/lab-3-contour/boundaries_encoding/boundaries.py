import math
from skimage import io

WHITE_CONTOUR_MASK = [200, 200, 200]
RED = 0
GREEN = 1
BLUE = 2

FILENAME = 'circle.png'
IMG = io.imread(FILENAME).tolist()

STEP = 5


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

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.pt_from} -> {self.pt_to}'


class ThreeAttributeEncoding:
    """
    Кодирование по трем признакам
    """

    def __init__(self, v_curr: Vector, v_next: Vector):
        self.len = None
        self.direction = None
        self.angle = None


def counterclockwise_angle(p1, p2):
    v1_theta = math.atan2(p1[1], p1[0])
    v2_theta = math.atan2(p2[1], p2[0])
    r = (v2_theta - v1_theta) * (180.0 / math.pi)
    if r < 0:
        r += 360.0
    return r


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
        # среди соседей будут 2 контурных точки - prev и next
        for p in nearest_points:
            # нам надо найти НЕ предыдущую точку
            if p.is_contour and p != prev_pt:
                next_pt = p
                break
        else:
            raise ValueError(f'could not find next point for {curr_pt}')

        prev_pt = curr_pt
        curr_pt = next_pt

        # контур замкнулся
        if curr_pt == start_pt:
            break

    print(f'contour size: {len(contour)}')
    print(f'start pt: {start_pt}')
    print(f'end   pt: {contour[len(contour) - 1]}')
    return contour


def approximate_contour(contour: list, step) -> list:
    vectors = []
    ctr_len = len(contour)
    for i in range(0, ctr_len, step):
        p0 = contour[i]
        p1 = contour[(i + step) % ctr_len]
        vectors.append(Vector(p0, p1))
    return vectors


def main():
    points = extract_contour_pts(IMG)
    vectors = approximate_contour(points, STEP)
    print('asdas')


if __name__ == '__main__':
    main()
