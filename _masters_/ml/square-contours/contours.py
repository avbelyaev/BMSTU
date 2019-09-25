import random

import cv2
import math

import numpy as np
from webcolors import name_to_rgb  # TODO replace with MPL colors

SOURCE_IMG_PATH = 'squares.jpg'
RESULT_IMG_PATH = 'res.png'

COLORS = ['white', 'green', 'purple', 'black', 'blue',
          'red', 'yellow', 'orange', 'brown', 'cyan',
          'Fuchsia', 'LawnGreen']


class Point:
    def __init__(self, x: int, y: int):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return other is not None and \
               self.x == other.x and self.y == other.y

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'[{self.x}:{self.y}]'


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


class Figure:
    def __init__(self, contour: np.ndarray = None, points: 'list[Point]' = None):
        self.color_name = COLORS[random.randint(0, len(COLORS) - 1)]
        self.contour = contour
        self.points = points
        self.peri, self.area, self.approx = None, None, None

        # если были переданы 3 точки, то предсказываем 4ю на основании других
        if points is not None and len(points) == 3:
            p1 = points[0]
            p2 = points[1]
            p3 = points[2]
            p4 = self.predict_point(p1, p2, p3)
            self.points = [p1, p2, p3, p4]

    @staticmethod
    def predict_point(p1: Point, p2: Point, p3: Point) -> Point:
        delta_x = p3.x - p2.x
        delta_y = p3.y - p2.y
        new_x = p1.x + delta_x
        new_y = p1.y + delta_y
        return Point(new_x, new_y)

    def approximate(self):
        self.area = cv2.contourArea(self.contour)
        self.peri = cv2.arcLength(self.contour, closed=True)
        self.approx = cv2.approxPolyDP(self.contour, epsilon=0.01 * self.peri, closed=True)

    def is_polyhedron(self) -> bool:
        is_not_small = self.area > 20
        return is_polyhedron(self.approx) and is_not_small

    def draw_on_current_canvas(self, img):
        print(f'drawing {self.color_name}')

        # точки рисуем как линии
        if self.points is not None:
            for i in range(len(self.points)):
                p1 = self.points[i]
                p2 = self.points[(i + 1) % len(self.points)]
                cv2.line(img, (p1.x, p1.y), (p2.x, p2.y), self.color_bgr, thickness=2)
            return img

        # контуры рисуем как контуры
        if self.contour is not None:
            draw_all_contours = -1
            cv2.drawContours(img, [self.approx], draw_all_contours, self.color_bgr, thickness=2)
            return img

    def split_into_squares(self) -> 'list[Figure]':
        print(f'> splitting {self.color_name}')

        # convert numpy points to normal points
        points = []
        for np_point in self.approx:
            points.append(Point(np_point[0, 0], np_point[0, 1]))

        squares = []
        angles = []
        for i in range(len(points)):
            p1 = points[i]
            p2 = points[(i + 1) % len(points)]
            p3 = points[(i + 2) % len(points)]

            v1 = Vector(p2, p1)
            v2 = Vector(p2, p3)

            curr_angle = directionwise_angle(v1, v2)
            angles.append(curr_angle)

            # как только набралось больше 2х углов, начинаем проверять
            if len(angles) > 2:
                prev_angle = angles[i - 1]
                prev_prev_angle = angles[i - 2]
                # соседние углы в сумме ~= 180'
                if 175 < (curr_angle + prev_angle) < 185:
                    # противоположные углы примерно похожи
                    if abs(prev_prev_angle - curr_angle) < 15:
                        figure = Figure(points=[points[i - 1], points[i], points[(i + 1) % len(points)]])
                        squares.append(figure)
                        print('>>>> splitted!')
        return squares

    @property
    def color_bgr(self) -> tuple:
        rgb = name_to_rgb(self.color_name)
        return rgb[2], rgb[1], rgb[0]

    def __str__(self):
        return f'{self.color_name}: area: {self.area}\n' \
               f'   peri: {self.peri}\n' \
               f'   points: {len(self.approx)}'


is_square = lambda points: len(points) == 4
is_polyhedron = lambda points: len(points) >= 4


def directionwise_angle(v1: Vector, v2: Vector) -> float:
    p1, p2 = v1.to_radius_vector(), v2.to_radius_vector()

    v1_theta = math.atan2(p1.y, p1.x)
    v2_theta = math.atan2(p2.y, p2.x)
    r = (v2_theta - v1_theta) * (180.0 / math.pi)
    if r < 0:
        r += 360.0
    return r


def main():
    img = cv2.imread(SOURCE_IMG_PATH, cv2.IMREAD_GRAYSCALE)
    original = cv2.imread(SOURCE_IMG_PATH)

    img = cv2.threshold(img, thresh=220, maxval=230, type=cv2.THRESH_TOZERO)[1]

    contours, hierarchy = cv2.findContours(img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    figures = []
    for c in contours:
        f = Figure(contour=c)
        f.approximate()

        if f.is_polyhedron():
            if is_square(f.approx):
                figures.append(f)
            else:
                squares = f.split_into_squares()
                figures.extend(squares)

    # draw
    for f in figures:
        original = f.draw_on_current_canvas(original)
    cv2.imwrite(f'./{RESULT_IMG_PATH}', original)


if __name__ == '__main__':
    main()
