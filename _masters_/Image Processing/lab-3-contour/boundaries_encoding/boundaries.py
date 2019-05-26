import math


class Point:
    def __init__(self, x: int, y: int):
        self.x = x
        self.y = y

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'[{self.x}:{self.y}]'


class Vector:
    def __init__(self, pt_from: Point, pt_to: Point):
        self.pt_from = pt_from
        self.pt_to = pt_to


class ThreeAttributeEncoding:
    """
    Кодирование по трем признакам
    """

    def __init__(self, v_curr: Vector, v_next: Vector):
        self.len = None
        self.direction = None
        self.angle = None


def angle_between(p1, p2):
    v1_theta = math.atan2(p1[1], p1[0])
    v2_theta = math.atan2(p2[1], p2[0])
    r = (v2_theta - v1_theta) * (180.0 / math.pi)
    if r < 0:
        r += 360.0
    return r


def main():
    pass


if __name__ == '__main__':
    main()
