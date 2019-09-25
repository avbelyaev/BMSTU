import random

import cv2
import math
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


class SimpleSquare:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.color_name = COLORS[random.randint(0, len(COLORS) - 1)]
        p4 = self.predict_point(p1, p2, p3)
        self.pts = [p1, p2, p3, p4]

    def draw_on_current_canvas(self, img):
        line_thickness = 2
        for i in range(len(self.pts)):
            p1 = self.pts[i]
            p2 = self.pts[(i + 1) % len(self.pts)]
            cv2.line(img, (p1.x, p1.y), (p2.x, p2.y), self.color_bgr, line_thickness)
        return img

    @property
    def color_bgr(self) -> tuple:
        rgb = name_to_rgb(self.color_name)
        return rgb[2], rgb[1], rgb[0]

    def predict_point(self, p1: Point, p2: Point, p3: Point) -> Point:
        delta_x = p3.x - p2.x
        delta_y = p3.y - p2.y
        new_x = p1.x + delta_x
        new_y = p1.y + delta_y
        return Point(new_x, new_y)


class Figure:
    def __init__(self, contour):
        self.color_name = COLORS[random.randint(0, len(COLORS) - 1)]
        self.orginal_contour = contour
        # initial approximation
        self.area = cv2.contourArea(contour)
        self.peri = cv2.arcLength(contour, closed=True)
        self.approx = cv2.approxPolyDP(contour, epsilon=0.01 * self.peri, closed=True)

    def is_polyhedron(self) -> bool:
        is_not_empty = self.area > 20
        return len(self.approx) >= 4 and is_not_empty

    def draw_on_current_canvas(self, img):
        print(f'drawing {self.color_name}')
        draw_all_contours = -1
        thickness = 2
        cv2.drawContours(img, [self.approx], draw_all_contours, self.color_bgr, thickness)
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

            # try create square
            if len(angles) > 4:
                prev_angle = angles[i - 1]
                prev_prev_angle = angles[i - 2]
                if 175 < (curr_angle + prev_angle) < 185:
                    if abs(prev_prev_angle - curr_angle) < 11:
                        print('>>>> bingo')
                        figure = SimpleSquare(points[i - 1], points[i], points[(i + 1) % len(points)])
                        squares.append(figure)
        return squares

    @property
    def color_bgr(self) -> tuple:
        rgb = name_to_rgb(self.color_name)
        return rgb[2], rgb[1], rgb[0]

    def __str__(self):
        return f'{self.color_name}: area: {self.area}\n' \
               f'   peri: {self.peri}\n' \
               f'   points: {len(self.orginal_contour)} -> {len(self.approx)}'


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

    img = cv2.threshold(img, 220, 230, cv2.THRESH_TOZERO)[1]

    contours, hierarchy = cv2.findContours(img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    figures = []
    for ctr in contours:
        f = Figure(ctr)

        if f.is_polyhedron():
            if len(f.approx) == 4:
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
