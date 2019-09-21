import cv2
import numpy as np

SOURCE_IMG_PATH = 'squares.jpg'
RESULT_IMG_PATH = 'res.png'


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
    contours = sorted(contours, key=lambda cnt: cv2.arcLength(cnt, closed=True), reverse=True)

    squares = []
    for contour in contours:
        perimeter = cv2.arcLength(contour, closed=True)
        area = cv2.contourArea(contour)
        print(f'perimeter: {perimeter}, area: {area}')

        approximated_controur = cv2.approxPolyDP(contour, epsilon=0.001*perimeter, closed=True)

        if len(approximated_controur) >= 4 and perimeter < 500:
            squares.append(approximated_controur)

    print(len(squares))

    drawAllContours = -1
    colorBGR = (0, 0, 255)
    thickness = 1
    cv2.drawContours(original, squares, drawAllContours, colorBGR, thickness)

    cv2.imwrite(f'./{RESULT_IMG_PATH}', original)


if __name__ == '__main__':
    main()
