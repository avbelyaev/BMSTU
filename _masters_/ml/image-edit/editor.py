from enum import Enum
from typing import Tuple

import cv2
import numpy as np

SOURCE_IMG_PATH = 'face.png'
RESULT_IMG_PATH = 'res.png'


class RotationPoint(Enum):
    CENTER = 0
    TOP_LEFT = 1
    TOP_RIGHT = 2
    BOT_LEFT = 3
    BOT_RIGHT = 4


class ImgEditor:
    """
    Manipulate the image in the Builder-pattern-like way
    """
    def __init__(self, img_path: str):
        self._img_path = img_path
        self._img = cv2.imread(img_path)
        self._rows = self._img.shape[0]
        self._cols = self._img.shape[1]
        self._shape = (self._cols, self._rows)

    def rotate(self, point: RotationPoint = RotationPoint.CENTER, angle: int = 0) -> 'ImgEditor':
        point = self._get_rotation_point(point)
        rot_matrix = cv2.getRotationMatrix2D(point, angle, 1)
        self._img = cv2.warpAffine(self._img, rot_matrix, self._shape)
        return self

    def translate(self, x_offset: int = 0, y_offset: int = 0) -> 'ImgEditor':
        translation_matrix = np.float32([[1, 0, x_offset], [0, 1, -y_offset]])
        self._img = cv2.warpAffine(self._img, translation_matrix, self._shape)
        return self

    def shrink(self, x_factor: float = 1.0, y_factor: float = 1.0) -> 'ImgEditor':
        self._img = cv2.resize(self._img, None,
                               fx=x_factor, fy=y_factor,
                               interpolation=cv2.INTER_CUBIC)
        return self

    def mirror(self, over_x: bool = False, over_y: bool = False) -> 'ImgEditor':
        if over_y:
            self._img = cv2.flip(self._img, 1)
        if over_x:
            self._img = cv2.flip(self._img, 0)
        return self

    def save(self, result_path: str):
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


if __name__ == '__main__':
    ImgEditor(SOURCE_IMG_PATH) \
        .rotate(RotationPoint.CENTER, angle=10) \
        .rotate(RotationPoint.CENTER, angle=15) \
        .shrink(x_factor=0.9) \
        .translate(y_offset=20) \
        .mirror(over_y=True) \
        .save(RESULT_IMG_PATH)
