import random

import cv2
from webcolors import name_to_rgb   # TODO replace with MPL colors

SOURCE_IMG_PATH = 'squares.jpg'
RESULT_IMG_PATH = 'res.png'

COLORS = ['white', 'green', 'purple', 'black', 'blue',
          'red', 'yellow', 'orange', 'brown', 'cyan',
          'Fuchsia', 'LawnGreen']


class Figure:
    def __init__(self, contour):
        self.contour = contour
        self.peri = cv2.arcLength(contour, closed=True)
        self.area = cv2.contourArea(contour)
        self.approx = cv2.approxPolyDP(contour, epsilon=0.001 * self.peri, closed=True)
        self.color_name = COLORS[random.randint(0, len(COLORS) - 1)]

    @property
    def color_bgr(self) -> tuple:
        rgb = name_to_rgb(self.color_name)
        return rgb[2], rgb[1], rgb[0]

    def __eq__(self, other):
        return self.area - other.area < 100

    def __str__(self):
        return f'{self.color_name}: peri: {self.peri}, area: {self.area}, ' \
               f'points: {len(self.contour)} -> {len(self.approx)}'


def random_color_name():
    return 'green'


def draw_squares(image, figures: 'List[Figure]'):
    for fig in figures:
        draw_all_contours = -1
        thickness = 2
        cv2.drawContours(image, [fig.contour], draw_all_contours, fig.color_bgr, thickness)

    cv2.imwrite(f'./{RESULT_IMG_PATH}', image)


def already_exists(new_figure: Figure, figures: list):
    for fig in figures:
        if abs(fig.area - new_figure.area) < 100:
            return True
    return False


def main():
    img = cv2.imread(SOURCE_IMG_PATH)
    original = img.copy()

    # img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # img = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 11, 0)
    # img = cv2.bilateralFilter(img, 110, 150, 50)

    kernelSize = (5, 5)
    img = cv2.GaussianBlur(img, kernelSize, 0)

    thresh = 250
    # img = cv2.threshold(img, thresh, 255, cv2.THRESH_MASK)[1]
    edged = cv2.Canny(img, threshold1=40, threshold2=200)

    # img = cv2.medianBlur(img, 5)[1]

    # kernel = np.ones((3, 3), np.uint8)
    # img = cv2.erode(img, kernel, iterations=1)
    # img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)

    contours, hierarchy = cv2.findContours(edged, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    figures = []
    for ctr in contours:
        fig = Figure(ctr)

        if not already_exists(fig, figures):
            figures.append(fig)
            print(fig)

    draw_squares(original, figures)


if __name__ == '__main__':
    main()
