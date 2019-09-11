from enum import Enum
from typing import Tuple

import cv2

SOURCE_IMG_PATH = 'yoba.png'
RESULT_IMG_PATH = 'res.png'


class RotationPoint(Enum):
    CENTER = 0
    TOP_LEFT = 1
    TOP_RIGHT = 2
    BOT_LEFT = 3
    BOT_RIGHT = 4


class ImgEditor:
    def __init__(self, img_path: str):
        self._img_path = img_path
        self._img = cv2.imread(img_path)
        self._rows = self._img.shape[0]
        self._cols = self._img.shape[1]
        self._shape = (self._cols, self._rows)

    def rotate(self, over: RotationPoint = RotationPoint.CENTER, angle: int = 0):
        """
        rotate image by angle over rotation point
        """
        point = self._get_rotation_point(over)
        rot_matrix = cv2.getRotationMatrix2D(point, angle, 1)
        self._img = cv2.warpAffine(self._img, rot_matrix, self._shape)

    def shrink(self, x_factor: float = 1.0, y_factor: float = 1.0):
        self._img = cv2.resize(self._img, None,
                               fx=x_factor, fy=y_factor,
                               interpolation=cv2.INTER_CUBIC)

    def apply_and_save(self, result_path: str):
        """
        apply changes and save at path relative to current directory
        """
        cv2.imwrite(f'./{result_path}', self._img)

    def _get_rotation_point(self, rotation_pt: RotationPoint) -> Tuple[int, int]:
        if rotation_pt is RotationPoint.CENTER:
            return self._cols / 2, self._rows / 2
        elif rotation_pt is RotationPoint.TOP_LEFT:
            return 0, 0
        elif rotation_pt is RotationPoint.TOP_RIGHT:
            return self._cols, 0
        elif rotation_pt is RotationPoint.BOT_LEFT:
            return 0, self._rows
        elif rotation_pt is RotationPoint.BOT_RIGHT:
            return self._cols, self._rows
        else:
            raise ValueError("unexpected rotation point!")


def main():
    img = ImgEditor(SOURCE_IMG_PATH)
    img.rotate(RotationPoint.CENTER, 10)
    img.shrink(x_factor=0.7, y_factor=0.7)
    img.apply_and_save(RESULT_IMG_PATH)


if __name__ == '__main__':
    main()
