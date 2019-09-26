import cv2
import numpy as np

from _masters_.ml.square_contours.contours import Point
# https://vk.com/doc1164151_516201556?hash=d3ac940079dcba822a&dl=37b24c9319a7ddb58b

SOURCE_IMG_PATH = 'blurred.jpg'
RESULT_IMG_PATH = 'res.png'

WINDOW_HEIGHT = 30
WINDOW_WIDTH = 30

FRAGMENT_HEIGHT = 80
FRAGMENT_WIDTH = 100


class Fragment:
    def __init__(self, top_left: Point, bot_right: Point, img: np.ndarray):
        self.top_left = top_left
        self.bot_right = bot_right
        self.img = img[top_left.x: bot_right.x, top_left.y: bot_right.y]

    def local_variance_at(self, m: int, n: int) -> float:
        summ = 0
        i = 0
        while i < WINDOW_WIDTH:
            j = 0
            while j < WINDOW_HEIGHT:
                w = Window(top_left=Point(i, j),
                           bot_right=Point(i + WINDOW_WIDTH, j + WINDOW_HEIGHT),
                           img=self.img)
                I_mean = w.mean
                I = self.img[m + i, n + j]
                summ += (I - I_mean) ** 2
                j += 1
            i += 1
        return summ / (WINDOW_HEIGHT * WINDOW_WIDTH)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.top_left}..{self.bot_right}\n'


class Window(Fragment):
    def __init__(self, top_left: Point, bot_right: Point, img: np.ndarray):
        super().__init__(top_left, bot_right, img)

    @property
    def mean(self) -> float:
        return self.img.mean()


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
    lv = frag.local_variance_at(0, 0)

    cv2.rectangle(original, (frag.top_left.x, frag.top_left.y),
                  (frag.bot_right.x, frag.bot_right.y), (0, 0, 255), thickness=2)

    cv2.imwrite(f'./{RESULT_IMG_PATH}', original)


if __name__ == '__main__':
    main()
