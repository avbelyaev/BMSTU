import cv2
import numpy as np

from _masters_.ml.square_contours.contours import Point

SOURCE_IMG_PATH = 'blurred.jpg'
RESULT_IMG_PATH = 'res.png'

WINDOW_HEIGHT = 100
WINDOW_WIDTH = 100

FRAGMENT_HEIGHT = 80
FRAGMENT_WIDTH = 100


class Fragment:
    def __init__(self, top_left: Point, bot_right: Point, img: np.ndarray):
        self.top_left = top_left
        self.bot_right = bot_right
        self.matrix = img[top_left.x: bot_right.x, top_left.y: bot_right.y]

    def local_variance(self) -> float:
        sum = 0
        i = 0
        while i < WINDOW_WIDTH:
            j = 0
            while j < WINDOW_HEIGHT:
                pass

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.top_left}..{self.bot_right}\n'


def crop(img) -> 'list[Fragment]':
    fragments = []
    img_h, img_w = img.shape[:2]
    y = 0
    while y < img_w:
        x = 0
        while x < img_h:
            fragment = Fragment(top_left=Point(x, y),
                                bot_right=Point(x + FRAGMENT_WIDTH, y + FRAGMENT_HEIGHT),
                                img=img)
            fragments.append(fragment)
            x += FRAGMENT_HEIGHT
        y += FRAGMENT_WIDTH
    return fragments


def main():
    img = cv2.imread(SOURCE_IMG_PATH, cv2.IMREAD_GRAYSCALE)
    original = cv2.imread(SOURCE_IMG_PATH)
    # img = cv2.resize(img, )
    fragments = crop(img)

    print(fragments)

    frag = fragments[3]
    cv2.rectangle(original, (frag.top_left.x, frag.top_left.y),
                  (frag.bot_right.x, frag.bot_right.y), (0, 0, 255), thickness=2)

    # size = imglen // WINDOW_SIZE
    # deviations = [[0] * size for _ in range(size)]
    #
    # print(f'computing brightness and deviation for img[{imglen}:{imglen}]')
    cv2.imwrite(f'./{RESULT_IMG_PATH}', original)


if __name__ == '__main__':
    main()
