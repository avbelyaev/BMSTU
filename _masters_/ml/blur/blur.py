import cv2
import numpy as np
from scipy.signal import convolve2d

from _masters_.ml.square_contours.contours import Point

# https://vk.com/doc1164151_516201556?hash=d3ac940079dcba822a&dl=37b24c9319a7ddb58b

SOURCE_IMG_PATH = 'original.png'
RESULT_TENG_IMG_PATH = 'teng-res.png'
RESULT_GLVN_IMG_PATH = 'glvn-res.png'

# размер окна, проходя которым по фрагменту, проверяем его заблюренность
WINDOW_HEIGHT = 3
WINDOW_WIDTH = 3

# фрагменты, на которые бьется изображение
FRAGMENT_HEIGHT = 30
FRAGMENT_WIDTH = 30

SOBEL_X = [[-1, 0, 1],
           [-2, 0, 2],
           [-1, 0, 1]]
SOBEL_Y = [[1, 2, 1],
           [0, 0, 0],
           [-1, -2, -1]]


class Fragment:
    def __init__(self, top_left: Point, bot_right: Point, img: np.ndarray):
        self.top_left = top_left
        self.bot_right = bot_right
        self.img = img[top_left.y: bot_right.y, top_left.x: bot_right.x]
        self.height, self.width = self.img.shape[:2]
        self.teng = None    # for memoization
        self.glvn = None    # for memoization
        self.lvs = dict()

    def tenengrad(self) -> float:
        if self.teng is not None:
            return self.teng
        G_x = convolve2d(self.img, SOBEL_X)
        G_y = convolve2d(self.img, SOBEL_Y)
        S = np.sqrt(np.square(G_x) + np.square(G_y))

        threshold = np.max(S) * 0.66
        res = np.where(S > threshold)

        self.teng = np.sum(np.square(res))
        return self.teng

    def grey_level_local_variance(self) -> float:
        if self.glvn is not None:
            return self.glvn
        lv = 0
        M = self.width
        N = self.height
        m = 0
        while m < M - WINDOW_HEIGHT:
            n = 0
            while n < N - WINDOW_WIDTH:
                lv += self._local_variance_at(n, m)
                n += 1
            m += 1

        lv_mean = lv / (M * N)

        m = 0
        glvn = 0
        while m < M - WINDOW_HEIGHT:
            n = 0
            while n < N - WINDOW_WIDTH:
                lv = self._local_variance_at(n, m)
                glvn += (lv - lv_mean) ** 2
                n += 1
            m += 1

        self.glvn = glvn / (M * N)
        return self.glvn

    def _local_variance_at(self, m: int, n: int) -> float:
        if (m, n) in self.lvs:
            return self.lvs[(m, n)]
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

        self.lvs[(m, n)] = summ / (WINDOW_HEIGHT * WINDOW_WIDTH)
        return self.lvs[(m, n)]

    def draw(self, img: np.ndarray):
        bgr_orange_color = (0, 128, 255)
        cv2.rectangle(img,
                      pt1=(self.top_left.x, self.top_left.y),
                      pt2=(self.bot_right.x, self.bot_right.y),
                      color=bgr_orange_color, thickness=1)

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
    while y < img_h:
        x = 0
        while x < img_w:
            fragment = Fragment(top_left=Point(x, y),
                                bot_right=Point(x + FRAGMENT_WIDTH, y + FRAGMENT_HEIGHT),
                                img=img)
            fragments.append(fragment)
            x += FRAGMENT_WIDTH
        y += FRAGMENT_HEIGHT
    return fragments


def main():
    # read imgs
    original1 = cv2.imread(SOURCE_IMG_PATH)
    original2 = original1.copy()
    img = cv2.imread(SOURCE_IMG_PATH, cv2.IMREAD_GRAYSCALE)

    # crop into pieces
    fragments = crop(img)

    # detect blurred pieces with tenengrad
    tenengrads = list(map(lambda f: f.tenengrad(), fragments))
    teng_thresh = max(tenengrads) * 0.66

    for fragment in fragments:
        if fragment.tenengrad() > teng_thresh:
            fragment.draw(original1)

    cv2.imwrite(f'./{RESULT_TENG_IMG_PATH}', original1)

    # detect blur with GLVN
    glvns = []
    i = 0
    for f in fragments:
        glvns.append(f.grey_level_local_variance())
        print(f'{i}/{len(fragments)}: {f.grey_level_local_variance()}')
        i += 1
    glvn_thresh = max(glvns) * 0.001

    for fragment in fragments:
        if fragment.grey_level_local_variance() < glvn_thresh:
            fragment.draw(original2)

    cv2.imwrite(f'./{RESULT_GLVN_IMG_PATH}', original2)


if __name__ == '__main__':
    main()
