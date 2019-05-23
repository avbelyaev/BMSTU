import math

import numpy as np
from skimage import io

# Выделение границ статистическим методом
# https://stanislaw.ru/rus/education/university/imagine.asp?section=imageborders


FILENAME = 'regular_002.jpg'
WINDOW_SIZE = 2

# трешхолдинг можно убрать
THRESHOLD = 0.06
CONTOUR = 1
BACKGR = 0


def deviation(x: int, y: int, img: list) -> float:
    def avg_brightness(x: int, y: int, img: list) -> float:
        s = 0
        for i in range(x, x + WINDOW_SIZE):
            for j in range(y, y + WINDOW_SIZE):
                s += img[i][j]
        return s / (WINDOW_SIZE * WINDOW_SIZE)

    avg = avg_brightness(x, y, img)
    s = 0
    for i in range(x, x + WINDOW_SIZE):
        for j in range(y, y + WINDOW_SIZE):
            s += (img[i][j] - avg) ** 2

    return math.sqrt(s / (WINDOW_SIZE * WINDOW_SIZE))


def main():
    img = io.imread(FILENAME, as_gray=True)
    img = img.tolist()
    imglen = len(img)

    # создать таблицу средних яркостей в частях изображения размером с рабочее окно
    size = imglen // WINDOW_SIZE
    deviations = [[0] * size for _ in range(size)]

    print(f'computing brightness and deviation for img[{imglen}:{imglen}]')
    i = 0
    k = 0
    while i < imglen:
        j = 0
        m = 0
        while j < imglen:
            deviations[k][m] = deviation(i, j, img)
            m += 1
            j += WINDOW_SIZE
        i += WINDOW_SIZE
        k += 1

    print('adjusting input image')
    for i in range(imglen):
        for j in range(imglen):
            img[i][j] *= deviations[i // WINDOW_SIZE][j // WINDOW_SIZE]
            img[i][j] = CONTOUR if img[i][j] > THRESHOLD else BACKGR

    img_array = np.array(img)
    io.imsave('res.png', img_array)


if __name__ == '__main__':
    main()
