import cv2

from _masters_.ml.square_contours.contours import Point

SOURCE_IMG_PATH = 'blurred.jpg'
RESULT_IMG_PATH = 'res.png'
WINDOW_SIZE = 100


class Fragment:
    def __init__(self, top_left: Point, bot_right: Point):
        self.topleft = top_left
        self.botright = bot_right

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'{self.topleft}..{self.botright}\n'


def crop(img) -> 'list[Fragment]':
    fragments = []
    img_h, img_w = img.shape[:2]
    i = 0
    while i < img_h:
        j = 0
        while j < img_w:
            fragment = Fragment(Point(i, j), Point(i + WINDOW_SIZE, j + WINDOW_SIZE))
            fragments.append(fragment)
            j += WINDOW_SIZE
        i += WINDOW_SIZE
    return fragments


def main():
    img = cv2.imread(SOURCE_IMG_PATH, cv2.IMREAD_GRAYSCALE)
    original = cv2.imread(SOURCE_IMG_PATH)
    # img = cv2.resize(img, )
    fragments = crop(img)

    print(fragments)

    frag = fragments[15]
    cv2.rectangle(original, (frag.topleft.x, frag.topleft.y),
                  (frag.botright.x, frag.botright.y), (0, 0, 255), thickness=2)

    # size = imglen // WINDOW_SIZE
    # deviations = [[0] * size for _ in range(size)]
    #
    # print(f'computing brightness and deviation for img[{imglen}:{imglen}]')
    cv2.imwrite(f'./{RESULT_IMG_PATH}', original)


if __name__ == '__main__':
    main()
