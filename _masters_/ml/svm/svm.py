import numpy as np
import matplotlib.pyplot as plt

DATASET_SIZE = 100

CLAZZ_1 = -1
CLAZZ_2 = +1


class SVM:
    def __init__(self):
        pass

    def fit(self, data: np.ndarray):
        all_data = []
        all_data.extend(data[CLAZZ_1])
        all_data.extend(data[CLAZZ_2])
        print('asd')

    def predict(self, data: np.ndarray):
        pass


def main():
    center1 = (50, 60)
    center2 = (80, 20)
    distance = 20

    x1 = np.random.uniform(low=center1[0], high=center1[0] + distance, size=(DATASET_SIZE,))
    y1 = np.random.normal(loc=center1[1], scale=distance, size=(DATASET_SIZE,))

    x2 = np.random.uniform(center2[0], center2[0] + distance, size=(DATASET_SIZE,))
    y2 = np.random.normal(center2[1], distance, size=(DATASET_SIZE,))

    # plt.scatter(x1, y1)
    # plt.scatter(x2, y2)
    # plt.show()

    clazz1_points = []
    for i in range(DATASET_SIZE):
        x = x1[i]
        y = y1[i]
        clazz1_points.append([x, y])

    clazz2_points = []
    for i in range(DATASET_SIZE):
        x = x2[i]
        y = y2[i]
        clazz2_points.append([x, y])

    marked_up_data = {
        CLAZZ_1: clazz1_points,
        CLAZZ_2: clazz2_points
    }

    svm = SVM()
    svm.fit(marked_up_data)


if __name__ == '__main__':
    main()
