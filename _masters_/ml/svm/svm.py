import numpy as np
import matplotlib.pyplot as plt

DATASET_SIZE = 100


def main():
    center1 = (50, 60)
    center2 = (80, 20)
    distance = 20

    x1 = np.random.uniform(low=center1[0], high=center1[0] + distance, size=(DATASET_SIZE,))
    y1 = np.random.normal(loc=center1[1], scale=distance, size=(DATASET_SIZE,))

    x2 = np.random.uniform(center2[0], center2[0] + distance, size=(DATASET_SIZE,))
    y2 = np.random.normal(center2[1], distance, size=(DATASET_SIZE,))

    plt.scatter(x1, y1)
    plt.scatter(x2, y2)
    plt.show()


if __name__ == '__main__':
    main()
