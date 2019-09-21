import cv2
import numpy as np

SOURCE_IMG_PATH = 'squares.jpg'
RESULT_IMG_PATH = 'res.png'


class Figure:
    def __init__(self, contour):
        self.contour = contour
        self.peri = cv2.arcLength(contour, closed=True)
        self.area = cv2.contourArea(contour)
        self.approx = cv2.approxPolyDP(contour, epsilon=0.001 * self.peri, closed=True)


def draw_squares(image, figures: 'List[Figure]'):
    squares = list(map(lambda fig: fig.contour, figures))

    draw_all_contours = -1
    color_bgr = (0, 0, 255)
    thickness = 1
    cv2.drawContours(image, squares, draw_all_contours, color_bgr, thickness)

    cv2.imwrite(f'./{RESULT_IMG_PATH}', image)


def main():
    img = cv2.imread(SOURCE_IMG_PATH)
    original = img.copy()

    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # img = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 11, 0)
    # img = cv2.bilateralFilter(img, 110, 150, 50)

    kernelSize = (3, 3)
    img = cv2.GaussianBlur(img, kernelSize, 0)

    thresh = 250
    # img = cv2.threshold(img, thresh, 255, cv2.THRESH_MASK)[1]
    edged = cv2.Canny(img, threshold1=40, threshold2=355)

    # img = cv2.medianBlur(img, 5)[1]

    # kernel = np.ones((3, 3), np.uint8)
    # img = cv2.erode(img, kernel, iterations=1)
    # img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)

    contours, hierarchy = cv2.findContours(edged, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    figures = list(map(lambda ctr: Figure(ctr), contours))
    figures.sort(key=lambda fig: fig.area, reverse=True)

    draw_squares(original, figures)


if __name__ == '__main__':
    main()
