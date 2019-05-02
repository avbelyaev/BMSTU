import matplotlib.pyplot as plot
import numpy as np
import scipy.cluster.hierarchy as sch


POINTS = [
    (0.6, 1.9),
    (1.8, 1.6),
    (2.7, 2.0),
    (3.0, 2.1),
    (3.0, 2.6),
    (3.1, 4.5),
    (3.8, 0.6),
    (4.2, 2.7)
]

CLUSTERS_NUM = 2


class Cluster:
    def __init__(self):
        self.points = []

    def add_point(self, point: tuple):
        self.points.append(point)

    def coordinates(self):
        return list(point[0] for point in self.points), list(point[1] for point in self.points)


def main():
    d = sch.distance.pdist(np.array(POINTS), metric='cosine')
    L = sch.linkage(d, method='complete')
    clusterized = sch.fcluster(L, 0.5 * d.max(), 'distance')

    clusters = []
    for i in range(CLUSTERS_NUM):
        clusters.append(Cluster())

    i = 0
    for p in clusterized:
        clusters[p - 1].add_point(POINTS[i])
        i += 1

    xs1, ys1 = clusters[0].coordinates()
    plot.scatter(xs1, ys1)

    xs2, ys2 = clusters[1].coordinates()
    plot.scatter(xs2, ys2)

    plot.gca().legend(('C1', 'C2'))
    plot.show()


if __name__ == '__main__':
    main()
