import cv2

SOURCE_IMG_PATH = 'squares.jpg'
RESULT_IMG_PATH = 'res.png'


class ContourExtractor:
    def __init__(self, img_path: str):
        self._img_path = img_path
        self._img = cv2.imread(img_path)
        self._rows = self._img.shape[0]
        self._cols = self._img.shape[1]
        self._shape = (self._cols, self._rows)

    def extract_contour(self):
        self._img = cv2.cvtColor(self._img, cv2.COLOR_BGR2GRAY)

        # self._img = cv2.threshold(self._img, 50, 250, cv2.THRESH_BINARY)

        # self._img = cv2.bilateralFilter(self._img, 111, 17, 17)

        edged = cv2.Canny(self._img, 80, 500)
        cnts, heirarchy = cv2.findContours(edged.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

        # for c in cnts:
        #     # approximate the contour
        #     # peri = cv2.arcLength(c, True)
        #     approx = cv2.approxPolyDP(c, 0.015, True)
        #
        #     # if our approximated contour has four points, then
        #     # we can assume that we have found our screen
        #     if len(approx) == 4:
        #         screenCnt = approx
        #         break
        cv2.drawContours(self._img, cnts, -1, (0, 255, 0), 2)
        return self

    def save(self, result_path: str):
        cv2.imwrite(f'./{result_path}', self._img)


def main():
    ContourExtractor(SOURCE_IMG_PATH) \
        .extract_contour() \
        .save(RESULT_IMG_PATH)


if __name__ == '__main__':
    main()
