from skimage import io
import math

# берем первые X строк изображения
# идет по 2м изображениям (экспертному и машинному) построчно
# как только находит в строке пиксель соответствующего цвета (красного или зеленого)
# замеряет евклидово расстояние между такими пикселями первого и второго изображения
# строит таблицу расхождений

DEVIATION_MIN_DIST = 0.1
DEVIATION_MAX_DIST = 0.3

FILENAME_EXPERT = 'expert_31.jpg'
FILENAME_MACHINE = 'machine_31.jpg'

EXPERT_POINT = [0, 255, 0]  # green
MACHINE_POINT = [255, 0, 0]  # red

ROWS_TO_COMPARE = 30


class Point:
    def __init__(self, x: int, y: int):
        self.x = x
        self.y = y

    def euclidean_dist(self, another) -> float:
        return math.sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'[{self.x}:{self.y}]'


# отклонение
class Deviation:
    def __init__(self, expectation: Point, reality: Point):
        self.exp = expectation
        self.real = reality
        self.dist = expectation.euclidean_dist(reality)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'Exp: {self.exp} Real: {self.real} dist:{self.dist}'


def main():
    # экспертное изображение
    img_expert = io.imread(FILENAME_EXPERT)
    img_expert = img_expert.tolist()
    imglen_expert = len(img_expert)

    # полученное автоматическим путем изображение
    img_machine = io.imread(FILENAME_MACHINE)
    img_machine = img_machine.tolist()
    imglen_machine = len(img_machine)

    assert imglen_expert == imglen_machine

    deviations = []
    for i in range(ROWS_TO_COMPARE):
        # ищем точку, поставленную экспертом
        expert_point = None
        for j in range(imglen_expert):
            if EXPERT_POINT == img_expert[i][j]:
                expert_point = Point(i, j)

        # ищем точку, поставленную алгоритмом
        machine_point = None
        for j in range(imglen_machine):
            if MACHINE_POINT == img_machine[i][j]:
                machine_point = Point(i, j)

        if (expert_point is None) or (machine_point is None):
            raise AttributeError(f'Could not find contour: '
                                 f'expert:{expert_point}, machine:{machine_point}')

        deviations.append(Deviation(expert_point, machine_point))

    [print(d) for d in deviations]

    # находимо отклоенения меньше THRESHOLD
    deviations = list(filter(lambda d: DEVIATION_MIN_DIST < d.dist < DEVIATION_MAX_DIST, deviations))

    print(f'points with deviation in range: ({DEVIATION_MIN_DIST}, {DEVIATION_MAX_DIST})')
    [print(d) for d in deviations]


if __name__ == '__main__':
    main()
