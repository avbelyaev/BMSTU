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
        self.img_path = img_path
        self.img = cv2.imread(img_path)
        self.rows = self.img.shape[0]
        self.cols = self.img.shape[1]
        self.shape = (self.cols, self.rows)

    def rotate(self, over: RotationPoint = RotationPoint.CENTER, angle: int = 0):
        """
        rotate image by angle over rotation point
        """
        point = self._get_rotation_point(over)
        rot_matrix = cv2.getRotationMatrix2D(point, angle, 1)
        self.img = cv2.warpAffine(self.img, rot_matrix, self.shape)

    def shrink(self, x_factor: float = 1.0, y_factor: float = 1.0):
        self.img = cv2.resize(self.img, None,
                              fx=x_factor, fy=y_factor,
                              interpolation=cv2.INTER_CUBIC)

    def apply_and_save(self, result_path: str):
        """
        apply changes and save at path relative to current directory
        """
        cv2.imwrite(f'./{result_path}', self.img)

    def _get_rotation_point(self, rotation_pt: RotationPoint) -> Tuple[int, int]:
        if rotation_pt is RotationPoint.CENTER:
            return self.cols / 2, self.rows / 2
        elif rotation_pt is RotationPoint.TOP_LEFT:
            return 0, 0
        elif rotation_pt is RotationPoint.TOP_RIGHT:
            return self.cols, 0
        elif rotation_pt is RotationPoint.BOT_LEFT:
            return 0, self.rows
        elif rotation_pt is RotationPoint.BOT_RIGHT:
            return self.cols, self.rows
        else:
            raise ValueError("unexpected rotation point!")


def main():
    editor = ImgEditor(SOURCE_IMG_PATH)
    editor.rotate(RotationPoint.CENTER, 10)
    editor.shrink(x_factor=0.7, y_factor=0.7)
    editor.apply_and_save(RESULT_IMG_PATH)


if __name__ == '__main__':
    main()
