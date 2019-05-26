from skimage import io
import math
from prettytable import PrettyTable

# 1. размечаем в фотошопе границы текстуры в верхнем левом углу картинки
# 2. запускаем task3Machine.m -> сохраняем картинку

# алгоритм:
# берем первые X строк изображения
# идет по 2м изображениям (экспертному и машинному) построчно
# как только находит в строке пиксель соответствующего цвета (красного или зеленого)
# замеряет евклидово расстояние между такими пикселями первого и второго изображения
# строит таблицу расхождений

DEVIATION_MIN_DIST = 3.5
DEVIATION_MAX_DIST = 4.5

FILENAME_EXPERT = 'expert_31.png'
FILENAME_MACHINE = 'machine_31.png'

EXPERT_POINT_COLOR_THRESHOLD = [0, 210, 0]  # green
MACHINE_POINT_COLOR_THRESHOLD = [210, 0, 0]  # red

ROWS_TO_COMPARE = 30

RED = 0
GREEN = 1
BLUE = 2


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
        return f'Exp: {self.exp} Real: {self.real} dist: {self.dist}'


# цвет пикселя pixel[r,g,b] подходит под маску color[r,g,b]
# альфа-канал не учитываем
def is_of_color(color: list, pixel) -> bool:
    return pixel[RED] >= color[RED] and \
           pixel[GREEN] >= color[GREEN] and \
           pixel[BLUE] >= color[BLUE]


def main():
    # полученное автоматическим путем изображение
    img_machine = io.imread(FILENAME_MACHINE)
    img_machine = img_machine.tolist()
    imglen_machine = len(img_machine)

    # экспертное изображение
    img_expert = io.imread(FILENAME_EXPERT)
    img_expert = img_expert.tolist()
    imglen_expert = len(img_expert)

    deviations = []
    i = 0
    while i < ROWS_TO_COMPARE:
        # ищем точку, поставленную экспертом
        expert_point = None
        for j in range(imglen_expert):
            if is_of_color(EXPERT_POINT_COLOR_THRESHOLD, img_expert[i][j]):
                expert_point = Point(i, j)
                break

        # ищем точку, поставленную алгоритмом
        machine_point = None
        for j in range(imglen_machine):
            if is_of_color(MACHINE_POINT_COLOR_THRESHOLD, img_machine[i][j]):
                machine_point = Point(i, j)
                break

        if expert_point and machine_point:
            deviations.append(Deviation(expert_point, machine_point))
            i += 1

        else:
            print(f'skipping row {i}. expert: {expert_point}, machine: {machine_point}')
            i += 1

    [print(d) for d in deviations]

    # находим отклоенения попадающие под THRESHOLD
    deviations_thr = list(filter(lambda d: DEVIATION_MIN_DIST < d.dist < DEVIATION_MAX_DIST, deviations))

    print(f'\npoints with deviation in range: ({DEVIATION_MIN_DIST}, {DEVIATION_MAX_DIST})')
    [print(d) for d in deviations]

    # печать таблицы
    t = PrettyTable(['Ожидание', 'Реальность', 'Расхождение'])
    for d in deviations:
        t.add_row([d.exp, d.real, d.dist])

    print(t)

    print(f'Percentage with deviation: {len(deviations_thr)/len(deviations)}')


if __name__ == '__main__':
    main()
